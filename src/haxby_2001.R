########################################################
fmriPredictorMatrix <- function( fmri, mask , motionin , selector = NA , ncompcor= 3 , nsvdcomp = 4 )
{
  mat<-timeseries2matrix( fmri, mask )
  if ( is.na( selector ) ) selector<-rep(TRUE, 1:nrow( motionin ) )
  mat<-mat[ selector, ]
  msvd <- svd(t(motionin[selector,3:ncol(motionin) ]))
  motion <-as.matrix(msvd$v[, 1:nsvdcomp])
  cmpc<-as.matrix( compcor( mat, ncompcor = ncompcor ) )
  globsig<-apply( mat , FUN=mean , MARGIN=1 )
  if ( length(c(cmpc)) == nrow(mat) ) mat<-residuals( lm( mat ~ as.matrix(motion) + c(cmpc)  ) ) else mat<-residuals( lm( mat ~ as.matrix(motion) + cmpc + globsig ) )
  mydot<-F
  if ( mydot ) mat<-temporalwhiten( mat )
  return( mat )
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
if ( ! exists("myrates") ) myrates<-rep(NA,12)
for ( wrun in 1:1 )
{
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
doit<-T
if ( ! exists("fmri")  | doit )
  {
    fmri<-antsImageRead("AFFINE.nii.gz",4)
    fmriavg<-antsImageRead("AFFINE_avg.nii.gz",3)
    print(dim(fmri))
    mask<-antsImageRead('mask.nii.gz',3)
    dofeaturesel<-TRUE 
    if ( dofeaturesel ) mask<-getMask(fmriavg,250,1.e9,TRUE)
    design<-read.table('labels.txt',header=T)
    selector<-( as.numeric( design$chunks ) %% 2 == 0  )
    selector2<-( as.numeric( design$chunks ) %% 2 == 1  )
    selector<-( as.numeric( design$chunks ) != wrun   )
    selector2<-( as.numeric( design$chunks ) == wrun  )
    subdesign<-subset( design, selector )
    motionin<-read.csv('AFFINEMOCOparams.csv')
    ncc <- 3
    if ( dofeaturesel )  ncc <- 8 
    mat<-fmriPredictorMatrix( fmri, mask, motionin, selector, ncompcor = ncc )
  }
myclasses <- levels( subdesign$labels )
nclasses<-length(myclasses )
myblocks<-matrix( rep(0,(nrow(mat))*nclasses ), nrow=( nrow(mat)  ))
for ( i in 1:nclasses ) myblocks[,i]<-as.numeric(  subdesign$labels == myclasses[i] )
mysblocks<-myblocks
for ( i in 1:ncol(mysblocks) ) mysblocks[,i]<-predict(smooth.spline(mysblocks[,i],df=100))$y
myblocks<-cbind( myblocks, mysblocks )
nv<-100
if ( FALSE ) {
  wmat<-whiten(mat)
  ff<-sparseDecom2( inmatrix=list(wmat,as.matrix(myblocks)), inmask=list(mask,NA), perms=0, its=12, mycoption=1, sparseness=c( 0.02 , 0.1 ) , nvecs=nv, smooth=1, robust=0, cthresh=c(0,0), ell1 = 0.1 , z=-1 ) ;  mysccanimages<-imageListToMatrix( imageList=ff$eig1, mask=mask) 
  } else {
    if  ( dofeaturesel ) {
    # quick cca based data selection 
      print("begin feature selection")
      fsvecs<-9
      ff<-sparseDecom2( inmatrix=list(mat,as.matrix(myblocks)), inmask=list(mask,NA), perms=0, its=11, mycoption=1, sparseness=c( -0.002 , 1/fsvecs ) , nvecs=fsvecs, smooth=1, robust=0, cthresh=c(0,0), ell1 = 0.01 , z=-1 )
      sccamask<-antsImageClone( ff$eig1[[1]] )
      sccamask[ mask > 0 ]<-0
      for ( img in ff$eig1 ) {
        ImageMath(img@dimension,img,'abs',img)
        ImageMath(img@dimension,sccamask,'+',sccamask,img)
      }
      sccamask[ sccamask > 1.e-6 ]<-1 
      sccamask[ sccamask < 1     ]<-0
      antsImageWrite(sccamask,'sccamask.nii.gz')
      mask<-sccamask
      mat<-fmriPredictorMatrix( fmri, mask, motionin, selector )
    }
    ff<-svd( mat )
    mysccanimages<-t(ff$v[,1:nv])
    print( paste( "100:",sum(ff$d[1:100]/sum(ff$d)),"200:",sum(ff$d[1:200]/sum(ff$d))) )
  }
mysccanpreds <- ( mat  ) %*% t( mysccanimages )
mydf         <- data.frame( factpreds = as.factor((subdesign$labels))  , imgs = mysccanpreds )
my.rf        <- svm( factpreds ~ . , data=mydf, probability = TRUE  )
if ( ! exists("mat2") | doit ) 
  {
    ##### test phase ######
    mat2<-fmriPredictorMatrix( fmri, mask, motionin, selector2 )
  }
mysccanpreds2 <- ( mat2  ) %*% t( mysccanimages )
mydf2<-data.frame(  imgs = mysccanpreds2 ) # , corrs = cor( t( mysccanpreds2 ) ) )
mypred2<-predict( my.rf , newdata = mydf2 )
subdesign2<-subset( design, selector2 )
sublabels<-as.factor((subdesign2$labels))
mycomps<-rep(FALSE,length( subdesign2$labels  ))
for (  i in 1:length( subdesign2$labels  ) )
  {
    comp<-( sublabels[i] == mypred2[i] )
    mycomps[i]<-comp
  }
zz<-majoritylabel( sublabels , mypred2 )
myrate<-100*(sum(zz$groundtruth==zz$voted)/length(zz$groundtruth))
print(paste("CorrectClassify:",myrate,"%"))
myrates[ wrun+1 ]<-myrate
} # wrun loop
#######################


