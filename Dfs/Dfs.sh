#!/bin/bash
# P1 w.white-bearded, P2 brindled4(Etosha), P3 black
#Note that DFS is only useful if your sample size for populations P1 and P2 is at least 10 haploid genomes (5 diploids).

bcftools view -S /home/users/xiaodong/Documents/Project/Wildebeest/dfs/external/sample.black9_brindled4.list  /davidData/data/malthe/gt/results/wildebeest/vcf/Goat_wildebeest_variable_sites_mergehartebeest-topi_nomultiallelics_noindels_10dp_2het.bcf.gz | bcftools view -c 1 -o dfs.subset.vcf -O v

export PYTHONPATH=$PYTHONPATH:/home/users/xiaodong/Software/genomics_general-0.4
python3 ~/Software/genomics_general-0.4/VCF_processing/parseVCF.py  -i dfs.subset.vcf -o mydata.geno.brindled4.gz

python3  ~/Software/genomics_general-0.4/freq.py --threads 10 -g  mydata.geno.brindled4.gz -p wwbearded -p brindled4 -p black -p hartebeest --popsFile sample.black9_brindled4.popfile | gzip  > mydata.brindled4.basecounts.tsv.gz

#Note that your SFS must have the same number of samples for P1 and P2.
python3 ~/Software/genomics_general-0.4/sfs.py -i mydata.brindled4.basecounts.tsv.gz --inputType baseCounts --outgroup hartebeest --FSpops wwbearded brindled4 black --subsample 10 10 18 --pref brindled4 --suff .subsample10.sfs

