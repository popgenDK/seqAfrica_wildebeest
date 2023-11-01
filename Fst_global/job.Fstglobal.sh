#!/bin/bash

cd /home/users/xiaodong/Documents/Project/Wildebeest/Fst_global

 ~/Software/plink2 --bfile ./all.maf005.r099 --within ./all.maf005.r099.strata.tsv  --fst CATPHENO method=hudson --out Wildebeest.all.HudsonFst.global --allow-extra-chr 

sed -i.bak 's/Black/Yblack/g' Wildebeest.all.HudsonFst.global.fst.summary
