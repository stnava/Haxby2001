################################################
echo you are processing $PWD
fmri=bold.nii.gz
nm=AFFINE
ref=${nm}_avg.nii.gz
if [[ ! -s ${nm}MOCOparams.csv ]] ; then 
  antsMotionCorr -d 3 -a $fmri -o $ref
  antsMotionCorr  -d 3 -o [ ${nm}, ${nm}.nii.gz,${nm}_avg.nii.gz] -m MI[${ref}, ${fmri}, 1 , 32 , Regular, 0.1  ] -t Affine[ 0.1 ] -u 1 -e 1 -s 1x0 -f 2x1 -i 30x10 -n 3  
fi
cp ${nm}_avg.nii.gz mask.nii.gz
MultiplyImages 3 mask.nii.gz 0 mask.nii.gz 
for x in mask*nii.gz ; do 
  ImageMath 3 mask.nii.gz + mask.nii.gz $x 
done
ThresholdImage 3 mask.nii.gz mask.nii.gz 1 1.e9 

/Users/stnava/data/haxby2001/Haxby2001/src/haxby_2001.R 
echo your processing of $PWD is done 

