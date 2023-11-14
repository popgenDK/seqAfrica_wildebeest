#Only for Dsuite, for qpgraph, see the qpGraph folder

#!/bin/bash
Dsuite=/home/users/xiaodong/Software/Dsuite/Build/Dsuite
Bcftools=/home/users/xiaodong/Software/bcftools/bcftools
$Bcftools view /davidData/malthe/africa/gt/results/wildebeest/vcf/Goat_wildebeest_variable_sites_mergehartebeest-topi_nomultiallelics_noindels_10dp_2het.bcf.gz \
	   | $Dsuite Dtrios -l 61799287 stdin SETS.txt  -t TREE.nwk
$Dsuite Fbranch TREE.nwk SETS_tree.txt > Fbranch_out.txt
python3 /home/users/xiaodong/Software/Dsuite/utils/dtools.py  Fbranch_out.txt TREE.nwk  --tree-label-size  18 #max 20, but good for tree details
