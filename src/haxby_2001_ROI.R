#!/usr/bin/env Rscript
########################################################
# this assumes that we know the "block" stimuli   ######
# which allows us to majority vote over a block   ######
# to identify the final predictor for that time   ######
########################################################
library(caret)
library(randomForest)
library(ANTsR)
library(pheatmap)
library(e1071)
########################################################
########################################################
fmriPredictorMatrix <- function( fmri, mask , motionin , selector = NA , ncompcor= 3 , nsvdcomp = 4 )
{
  mat<-timeseries2matrix( fmri, mask )
  if ( is.na( selector ) ) selector<-rep(TRUE, nrow( motionin ) )
  mat<-mat[ selector, ]
  msvd <- svd(t(motionin[selector,3:ncol(motionin) ]))
  motion <-as.matrix(msvd$v[, 1:nsvdcomp])
  cmpc<-as.matrix( compcor( mat, ncompcor = ncompcor, fastsvd=TRUE ) )
  globsig<-apply( mat , FUN=mean , MARGIN=1 )
  if ( length(c(cmpc)) == nrow(mat) ) mat<-residuals( lm( mat ~ as.matrix(motion) ) ) else mat<-residuals( lm( mat ~ as.matrix(motion)  ) )
  mydot<-F
  if ( mydot ) mat<-temporalwhiten( mat )
  return( list( mat = mat, cmpc=cmpc, globsig=globsig )  )
}
########################################################
majoritylabel <- function( groundtruth, myprediction )
  {
  y  <- myprediction
  yt <- groundtruth
  ct<-1
  mylabs<-unique(  yt  )
  myyruns<-rep(NA,length(yt))
  myyruns[1]<-ct
  for ( i in 1:(length(y)-1) )
    {
    if ( yt[i] == yt[i+1] ) { myyruns[i+1]<-ct }  else { ct<-ct+1  ;   myyruns[i+1] <- ct }
    }
  myvotedlabels<-data.frame( groundtruth=rep(NA,ct) , voted=rep(NA,ct) )
  for ( i in unique( myyruns ) )
    {
    mylabs<-y[ myyruns == i ] 
    ff<-as.data.frame(table(mylabs))
    mylab<-ff$mylabs[ which.max( ff$Freq ) ]
    print( paste( i, mylab ,yt[myyruns == i][1]) )
    myvotedlabels$groundtruth[i]<-as.character(yt[myyruns == i][1])
    myvotedlabels$voted[i]<-as.character(mylab)
    }
  return(myvotedlabels)
  }
########################################################
if ( ! exists("myrates") ) myrates<-rep(NA,12)
design<-read.table('labels.txt',header=T)
unique(design$chunks)
runstotest<-unique(design$chunks)
runstotest<-runstotest[ runstotest < 12 ] 
if ( ! file.exists("AFFINE.nii.gz") )
  {
  print("FAILURE --- you need to be within a subject's directory")
  q()
  }
fmri<-antsImageRead("AFFINE.nii.gz",4)
fmriavg<-antsImageRead("AFFINE_avg.nii.gz",3)
motionin<-read.csv('AFFINEMOCOparams.csv')
############ other masking strategies ###########
maskNuis<-antsImageClone(fmriavg)
ThresholdImage(3,fmriavg,maskNuis,"Otsu",3)
ThresholdImage(3,maskNuis,maskNuis,1,1)
mask<-antsImageRead( 'mask.nii.gz', 3 )
if ( file.exists("AALlabel.nii.gz") & FALSE  )
  {
    aalimg<-antsImageRead("AALlabel.nii.gz",3)
    aalnums<-c( 55, 56, 89:100, 113 )
#   aalnums<-sort( unique( aalimg[ mask > 0 & aalimg > 0 ] ) )
    aalimg<-maskImage( aalimg,aalimg, as.list( aalnums ) , binarize=TRUE )
    aalvoxselection <- ( aalimg <= 100 & aalimg > 0  )
    mask<-antsImageClone( aalimg )
    mask[ !aalvoxselection ]<-0
  }
############ other masking strategies end ###########
library(irlba)
ncc <- 4
segmat<-timeseries2matrix( fmri, maskNuis )
segsvd<-irlba( segmat[,sample(1:ncol(segmat))[1:(ncol(segmat)/4)]] , nu=ncc )
fmripreds<-fmriPredictorMatrix( fmri, mask, motionin, selector=NA, ncompcor = ncc )
fullmat<-residuals( lm( fmripreds$mat ~ segsvd$u ) )
ff<-svd( fullmat )
for ( wrun in runstotest )
{
  mask<-antsImageRead('mask.nii.gz',3)
  selector<-( as.numeric( design$chunks ) != wrun   )
  selector2<-( as.numeric( design$chunks ) == wrun  )
  subdesign<-subset( design, selector )
  for ( nv in seq( from=175,to=175,by=5) )
    {
##### train phase ######
    mysccanimages<- (ff$u[selector,1:nv])
    whichtostudy<-( subdesign$labels != "rest" )
    mydf         <- data.frame( factpreds = as.factor((subdesign$labels))  , imgs = mysccanimages )
    mydf <- subset( mydf , whichtostudy )
    my.rf        <- svm( factpreds ~ . , data=mydf ) # , cross = 0 )
#    tempfact <-as.factor( rank( as.numeric( mydf$factpreds ) ) )
#    levels( tempfact )<-levels( mydf$factpreds )[ levels(mydf$factpreds ) != "rest" ]
#    mydf$factpreds <- tempfact
#    my.rf        <- randomForest( factpreds ~ . , data=mydf ) 
##### test phase ######
    whichtostudy2<-( subdesign2$labels != "rest" )
    mydf2<-data.frame(  imgs = ff$u[selector2,1:nv] ) 
    mydf2 <- subset( mydf2 , whichtostudy2 )
    mypred2<-predict( my.rf , newdata = mydf2 )
    subdesign2<-subset( design, selector2 )
    sublabels<-as.factor((subdesign2$labels))
    zz<-majoritylabel( sublabels[whichtostudy2] , mypred2 )
    myrate<-100*(sum(zz$groundtruth==zz$voted)/length(zz$groundtruth))
    print(paste("CorrectClassify:",myrate,"%",getwd(),wrun+1,' npreds ',nv) )
    }
  if ( wrun == 0 ) predictedlabels<-zz else predictedlabels<-rbind( predictedlabels, zz )
  myrates[ wrun+1 ]<-myrate
} # wrun loop
############################################################################################
ratedf<-data.frame( RunNumber=c(0:11), CrossValidatedPredictionForRun=myrates )
write.csv(ratedf,'mypredictionresults.csv',row.names=F)
write.csv(predictedlabels,'mypredictedLabels.csv',row.names=F)
