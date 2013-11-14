mydir=$PWD
for x in subj* ; do
 echo ${mydir}/${x}
 cd ${mydir}/${x} 
 ${mydir}/motion.sh 
done


