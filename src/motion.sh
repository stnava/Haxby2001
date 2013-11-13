
fmri=bold.nii.gz
nm=AFFINE
ref=${nm}_avg.nii.gz
antsMotionCorr -d 3 -a $fmri -o $ref
antsMotionCorr  -d 3 -o [ ${nm}, ${nm}.nii.gz,${nm}_avg.nii.gz] -m MI[${ref}, ${fmri}, 1 , 32 , Regular, 0.1  ] -t Affine[ 0.1 ] -u 1 -e 1 -s 1x0 -f 2x1 -i 30x10 -n 3  

