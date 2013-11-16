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
  if ( is.na( selector ) ) selector<-rep(TRUE, 1:nrow( motionin ) )
  mat<-mat[ selector, ]
  msvd <- svd(t(motionin[selector,3:ncol(motionin) ]))
  motion <-as.matrix(msvd$v[, 1:nsvdcomp])
  cmpc<-as.matrix( compcor( mat, ncompcor = ncompcor, fastsvd=FALSE ) )
  globsig<-apply( mat , FUN=mean , MARGIN=1 )
  if ( length(c(cmpc)) == nrow(mat) ) mat<-residuals( lm( mat ~ as.matrix(motion)   ) ) else mat<-residuals( lm( mat ~ as.matrix(motion)  ) )
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
maskFull<-getMask(fmriavg,100,1.e9,TRUE)
############ other masking strategies ###########
maskbrain<-antsImageRead("fullMask.nii.gz",3)
maskFull<-antsImageClone(maskbrain)
ImageMath(3,maskFull,"MD",maskFull,2)
ImageMath(3,maskFull,"-",maskFull,maskbrain)
# maskFull<-getMask(fmriavg,250,1.e9,TRUE)
# ImageMath(3,maskFull,"ME",maskFull,1)
# ImageMath(3,maskFull,"GetLargestComponent",maskFull)
antsImageWrite(maskFull , "../temp.nii.gz") 
mask<-antsImageRead( 'mask.nii.gz', 3 )
if ( file.exists("AALlabel.nii.gz") & FALSE  )
  {
    aalimg<-antsImageRead("AALlabel.nii.gz",3)
    aalnums<-c( 55, 56, 89:100, 113 )
#    aalnums<-sort( unique( aalimg[ mask > 0 & aalimg > 0 ] ) )
    aalimg<-maskImage( aalimg,aalimg, as.list( aalnums ) , binarize=TRUE )
    aalvoxselection <- ( aalimg <= 100 & aalimg > 0 & maskFull > 0 )
    mask<-antsImageClone( aalimg )
    mask[ !aalvoxselection ]<-0
  }
############ other masking strategies end ###########
# for ( wrun in runstotest )
for ( wrun in 11 )
{
  mask<-antsImageRead('mask.nii.gz',3)
  selector<-( as.numeric( design$chunks ) != wrun   )
  selector2<-( as.numeric( design$chunks ) == wrun  )
  subdesign<-subset( design, selector )
  ncc <- 4
  fmripreds<-fmriPredictorMatrix( fmri, mask, motionin, selector, ncompcor = ncc )
  fmripredsFull<-fmriPredictorMatrix( fmri, maskFull, motionin, selector, ncompcor = ncc )
  mat<-residuals( lm( fmripreds$mat ~ fmripredsFull$cmpc ) )
  myclasses <- levels( subdesign$labels )
#######################
##### train phase ######
#######################
  for ( nv in seq( from=100,to=100,by=5) )
    {
    ff<-svd( mat )
    mysccanimages<-t(ff$v[,1:nv])
    mysccanpreds <- ( mat  ) %*% t( mysccanimages )
    mydf         <- data.frame( factpreds = as.factor((subdesign$labels))  , imgs = mysccanpreds )
    my.rf        <- svm( factpreds ~ . , data=mydf, probability = TRUE ) # , kernel = "linear" ) # , cross = 0 )
#######################
##### test phase ######
#######################
  fmriTest     <- fmriPredictorMatrix( fmri, mask,     motionin, selector2, ncompcor = ncc )
  fmriTestFull <- fmriPredictorMatrix( fmri, maskFull, motionin, selector2, ncompcor = ncc )
  mat2<-residuals( lm( fmriTest$mat ~ fmriTestFull$cmpc ) )
  mysccanpreds2 <- ( mat2  ) %*% t( mysccanimages )
  mydf2<-data.frame(  imgs = mysccanpreds2 ) # , corrs = cor( t( mysccanpreds2 ) ) )
  mypred2<-predict( my.rf , newdata = mydf2 )
  subdesign2<-subset( design, selector2 )
  sublabels<-as.factor((subdesign2$labels))
  zz<-majoritylabel( sublabels , mypred2 )
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
