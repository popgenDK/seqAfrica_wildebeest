plink  --allow-extra-chr --bfile /home/users/long/Wildebeest/dataset/all.maf0.05.r2.0.99 --chr-set 29 --indep-pairwise 1000 100 0.8 --out all.maf0.05.r2.0.99.ld08
plink  --bfile /home/users/long/Wildebeest/dataset/all.maf0.05.r2.0.99 --chr-set 29 --extract all.maf0.05.r2.0.99.ld08.prune.in --make-bed --out all.maf0.05.r2.0.99.ld08.ldrune
for k in {2..13}
do
  bash multiConvV3.sh all.maf0.05.r2.0.99.ld08.ldrune 200 20 K$k $k
done
