#!/bin/bash
file=$1 #plink filename without .bed
num=$2 #number of iterations
P=$3 #number of threads
out=$4 #output directory (no /)
K=$5 #number of populations
star=$6
bfile=`basename $file`
echo num = $num
echo p = $P
echo out = $out
echo K = $K
if [ -f $file.bed ]
then
   echo file = $file
else
   echo File $file.bed does not exist.
fi

#touch $out/$bfile.likes.tmp
#rm $out/$bfile.likes.tmp

for f in `seq $star $(($star+$num-1))`
do
    echo -n -e $f"\t" >> $out/$bfile.likes.tmp
    /kellyData/home/users/xiaodong/Software/admixture_linux-1.3.0/admixture $file.bed $K -s $f -j$P > $out/$bfile.$K.log_$f
    mv $bfile.$K.Q $out/$bfile.$K.Q_$f
    mv $bfile.$K.P $out/$bfile.$K.P_$f
    grep ^Loglikelihood  $out/$bfile.$K.log_$f | cut -f2 -d" " >> $out/$bfile.likes.tmp
    CONV=`Rscript -e "r<-read.table('$out/$bfile.likes.tmp');r<-r[order(-r[,2]),];cat(sum(r[1,2]-r[,2]<5),'\n')"`

    if [ $CONV -gt 2 ]
    then
        cp $out/$bfile.$K.Q_$f $out/$bfile.$K.Q_conv
        cp $out/$bfile.$K.P_$f $out/$bfile.$K.P_conv
        cp $out/$bfile.$K.log_$f $out/$bfile.$K.log_conv
        echo "reaach likelihood duplicate times"
        break
    fi

done
cat $out/$bfile.likes.tmp | sort -k2 -n -r > $out/$bfile.likes
