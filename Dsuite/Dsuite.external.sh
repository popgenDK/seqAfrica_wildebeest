#!/bin/bash

Dsuite=/home/users/xiaodong/Software/Dsuite/Build/Dsuite
Bcftools=/home/users/xiaodong/Software/bcftools/bcftools



cd /home/users/xiaodong/Documents/Project/Wildebeest/Dsuite/external_nyassa

$Bcftools view /davidData/malthe/africa/gt/results/wildebeest/vcf/Goat_wildebeest_variable_sites_mergehartebeest-topi_nomultiallelics_noindels_10dp_2het.bcf.gz \
	   | $Dsuite Dtrios -l 61799287 stdin SETS.inNyassa.txt  -t TREE1_inNyassa.nwk 


cd /home/users/xiaodong/Documents/Project/Wildebeest/Dsuite/external

$Dsuite Fbranch TREE1_rmNyassa.nwk SETS.20220401_tree.txt > Fbranch_out.txt


python3 ~/Software/Dsuite/utils/dtools.py  Fbranch_out.txt TREE1_rmNyassa.nwk  --tree-label-size  8


