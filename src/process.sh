mydir=$PWD
for x in subj5 ; do
 echo ${mydir}/${x}
 cd ${mydir}/${x} 
 ${mydir}/motion.sh 
done


