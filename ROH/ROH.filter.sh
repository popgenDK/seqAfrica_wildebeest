#!/bin/bash

vcf1="Wildebeest_wildebeest_variable_sites_nomultiallelics_noindels_10dp_3het.MAF005.vcf.gz"
bcftools=/home/krishang/software/bcftools/bcftools-1.13/bin/bcftools
plink1=/home/users/xiaodong/Software/plink_linux_x86_64_20210606/plink

### step1, global filter, min DP of 10, min het of 3, allelic balance
$bcftools +setGT --help
$bcftools  +setGT  $vcf1  -- -t q -n . --include 'FMT/GT=="het" & ((FORMAT/AD[*:0] / (FMT/DP) ) <0.25 | (FMT/AD[*:0] / (FMT/DP)) >0.75)'   >  Wildebeest_wildebeest_variable_sites_nomultiallelics_noindels_10dp_3het_maf005_ab25.vcf

# convert to plink
~/Software/plink2  -vcf ./Wildebeest_wildebeest_variable_sites_nomultiallelics_noindels_10dp_3het_maf005_ab25.vcf \
		   --allow-extra-chr  --max-alleles 2 --const-fid 0  --set-missing-var-ids @:#   --make-bed --out Wildebeest_wildebeest_variable_sites_nomultiallelics_noindels_10dp_3het_maf005_ab25

mv Wildebeest_wildebeest_variable_sites_nomultiallelics_noindels_10dp_3het_maf005_ab25.fam ildebeest_wildebeest_variable_sites_nomultiallelics_noindels_10dp_3het_maf005_ab25.fam.bak
awk '{print $2,$2,$3,$4,$5,$6}' Wildebeest_wildebeest_variable_sites_nomultiallelics_noindels_10dp_3het_maf005_ab25.fam.bak > Wildebeest_wildebeest_variable_sites_nomultiallelics_noindels_10dp_3het_maf005_ab25.fam

### step2, population(subspecies) filter,
## extract sites of het rate >0.5


#  geno 0.05, basically no datamissing per snp
for i in *list
do
    pop=${i%%\.list}
    echo $pop
    # --ibc (ported from GCTA) calculates three inbreeding coefficients for each sample, and writes a report to plink.ibc. Briefly, Fhat1 is the usual variance-standardized relationship minus 1, Fhat2 is similar to the --het estimate, and Fhat3 is based on the correlation between uniting gametes.
    ${plink1} --bfile Wildebeest_wildebeest_variable_sites_nomultiallelics_noindels_10dp_3het_maf005_ab25 --keep $i  --allow-extra-chr --maf 0.05 --geno 0.05 --het --hardy --ibc --missing \
                      --set-missing-var-ids "@:#:\$1:\$2" --make-bed --out Wildebeest_wildebeest_variable_sites_nomultiallelics_noindels_10dp_3het_maf005_ab25.${pop}_maf005_geno005

    awk '$7>=0.5 {print $2}' Wildebeest_wildebeest_variable_sites_nomultiallelics_noindels_10dp_3het_maf005_ab25.${pop}_maf005_geno005.hwe | sed '1d' > Wildebeest_wildebeest_variable_sites_nomultiallelics_noindels_10dp_3het_maf005_ab25.${pop}_maf005_geno005.OHE05.sites
    echo

done

### step 3, individual plink
mkdir geno005_pop

for i in *list # loop through pop file
do
    pop=${i%%\.list}
    echo $pop
    for j in $(awk '{print $2}'  Wildebeest_wildebeest_variable_sites_nomultiallelics_noindels_10dp_3het_maf005_ab25.${pop}_maf005_geno005.fam) # loop through sample
    do
	sample=$j
	echo $sample
	if [ -d geno005_pop/${sample}_plink ]; then
        rm -r geno005_pop/${sample}_plink
	fi
	mkdir geno005_pop/${sample}_plink
	grep -w ${sample} Wildebeest_wildebeest_variable_sites_nomultiallelics_noindels_10dp_3het_maf005_ab25.${pop}_maf005_geno005.fam  > geno005_pop/${sample}_plink/${sample}.list
	# prepare ind plink file
	$plink1 --keep  geno005_pop/${sample}_plink/${sample}.list --allow-extra-chr --bfile Wildebeest_wildebeest_variable_sites_nomultiallelics_noindels_10dp_3het_maf005_ab25.${pop}_maf005_geno005 \
              --exclude Wildebeest_wildebeest_variable_sites_nomultiallelics_noindels_10dp_3het_maf005_ab25.${pop}_maf005_geno005.OHE05.sites  --keep-allele-order \
              --make-bed --homozyg --homozyg-window-het 3 --homozyg-window-missing 20 --out geno005_pop/${sample}_plink/${sample}_OHE05 &> geno005_pop/${sample}_plink/${sample}.log
    done
	
    
done

    
 
