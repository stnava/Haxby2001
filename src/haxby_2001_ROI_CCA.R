#!/usr/bin/env Rscript
#
# see http://cran.r-project.org/web/packages/bootfs/vignettes/bootfs.pdf
#
########################################################
# this assumes that we know the "block" stimuli   ######
# which allows us to majority vote over a block   ######
# to identify the final predictor for that time   ######
########################################################
library(irlba)
library(caret)
library(randomForest)
library(ANTsR)
library(pheatmap)
library(e1071)
data("aal",package='ANTsR')
########################################################
########################################################
fmriPredictorMatrix <- function( fmri, mask , motionin , selector = NA , ncompcor= 3 , nsvdcomp = 4 )
{
  mat<-timeseries2matrix( fmri, mask )
  if ( is.na( selector ) ) selector<-rep(TRUE, nrow( motionin ) )
  mat<-mat[ selector, ]
  msvd <- svd(t(motionin[selector,3:ncol(motionin) ]))
  motion <-as.matrix(msvd$v[, 1:nsvdcomp])
  cmpc<-as.matrix( compcor( mat, ncompcor = ncompcor, fastsvd=FALSE ) )
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
if ( file.exists("AALlabel.nii.gz") & TRUE  )
  {
  haxbymask<-antsImageRead("mask.nii.gz",3)
  aalimg<-antsImageRead("AALlabel.nii.gz",3)
  aalnums<-c( 55, 56, 89:92 )
  aalnums<-c( 1:92 )
  print( aal$label_name[aalnums] )
  aalimg<-maskImage( aalimg,aalimg, as.list( aalnums ) , binarize=TRUE )
  aalvoxselection <- ( aalimg <= 100 & aalimg > 0  )
  mask<-antsImageClone( aalimg )
  mask[ !aalvoxselection ]<-0
  ncc <- 4
  segmat<-timeseries2matrix( fmri, maskNuis )
  segsvd<-irlba( segmat[,sample(1:ncol(segmat))[1:(ncol(segmat)/4)]] , nu=ncc , nv = 0)
  fmripreds<-fmriPredictorMatrix( fmri, mask, motionin, selector=NA, ncompcor = ncc )
  fullmat<-residuals( lm( fmripreds$mat ~ segsvd$u ) )
# now cca
  myclasses <- levels( design$labels )
  nclasses<-length(myclasses )
  myblocks<-matrix( rep(0,(nrow(design))*nclasses ), nrow=( nrow(design)  ))
  for ( i in 1:nclasses ) myblocks[,i]<-as.numeric(  design$labels == myclasses[i] )
  for ( i in 1:ncol(myblocks) ) myblocks<-cbind(myblocks, predict(smooth.spline(myblocks[,i],df=100))$y )
  nv<-20
  whichtotest <- 2
  mselection<-( design$chunks != whichtotest   )
  sparseness <- (-1) * sum( haxbymask > 0 ) / sum( mask > 0 ) * 1.5 / nv
  mycca<-sparseDecom2( inmatrix=list(fullmat[mselection,],as.matrix(myblocks[mselection,])), inmask=list(mask,NA), perms=0, its=19, mycoption=1, sparseness=c( sparseness , 1/nv ) , nvecs=nv, smooth=1, robust=0, cthresh=c(5,0), ell1 = 0.1 , z=-1 )
  sccanmask<-eigSeg(mask,mycca$eig1)
  antsImageWrite(sccanmask,"../temp.nii.gz")
  sccanmask[ sccanmask > 0 & sccanmask < 19 ]<-1
  sccanmask[ sccanmask > 1.5 ]<-0
  sccanmat<-imageListToMatrix( mycca$eig1 , mask )
# now svd on cca output  
  }
nv<-min( c( 175, ncol( fullmat )*0.8 ) )
fmripreds<-fmriPredictorMatrix( fmri, sccanmask, motionin, selector=NA, ncompcor = ncc )
for ( wrun in runstotest )
  {
  selector<-( as.numeric( design$chunks ) != wrun   )
  selector2<-( as.numeric( design$chunks ) == wrun  )
  subdesign<-subset( design, selector )
  fullmat<-residuals( lm( fmripreds$mat[selector,] ~ segsvd$u[selector,] ) )
  mysvd<-svd( (fullmat) , nv = 0 )
  testmat<-residuals( lm( fmripreds$mat[selector2,] ~ segsvd$u[selector2,] ) )
  mytestimages<- testmat %*% t( t(mysvd$u[,1:nv]) %*% fullmat )
##### train phase ######
  mytrainimages<- (mysvd$u[,1:nv])
  whichtostudy<-( subdesign$labels != "rest" )
  mydf         <- data.frame( factpreds = as.factor((subdesign$labels))  , imgs = scale( mytrainimages ) )
  mydf <- subset( mydf , whichtostudy )
  my.rf        <- svm( factpreds ~ . , data=mydf ) # , cross = 0 )
##### test phase ######
  subdesign2<-subset( design, selector2 )
  whichtostudy2<-( subdesign2$labels != "rest" )
  mydf2<-data.frame(  imgs = scale( mytestimages ) ) # mysvd$u[selector2,1:nv] ) 
  mydf2 <- subset( mydf2 , whichtostudy2 )
  mypred2<-predict( my.rf , newdata = mydf2 )
  sublabels<-as.factor((subdesign2$labels))
  zz<-majoritylabel( sublabels[whichtostudy2] , mypred2 )
  myrate<-100*(sum(zz$groundtruth==zz$voted)/length(zz$groundtruth))
  print(paste("CorrectClassify:",myrate,"%",getwd(),wrun+1,' npreds ',nv) )
  if ( wrun == 0 ) predictedlabels<-zz else predictedlabels<-rbind( predictedlabels, zz )
  myrates[ wrun+1 ]<-myrate
} # wrun loop
print( paste( "%var used:",sum(mysvd$d[1:nv])/sum(mysvd$d) * 100 , " mean-rate ",mean(myrates)) )
############################################################################################
ratedf<-data.frame( RunNumber=c(0:11), CrossValidatedPredictionForRun=myrates )
write.csv(ratedf,'mypredictionresults.csv',row.names=F)
write.csv(predictedlabels,'mypredictedLabels.csv',row.names=F)
