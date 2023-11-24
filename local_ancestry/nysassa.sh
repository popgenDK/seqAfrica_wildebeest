#!/bin/bash

#conda activate loter

folder=/home/users/xiaodong/Documents/Project/Wildebeest/loter
for i in {1..29}
do
    if [[ "$i" -eq 3 ]]
    then
        :
    else
	# bcftools view /home/zilong/project/african1k/wildebeest/output/imputed/all.phased.maf0.05.r2.0.99.bcf -S wwhitebeard.list -r HiC_scaffold_${i} -Ov -o ${folder}/data/wwhitebeard.chr${i}.vcf

	# bcftools view /home/zilong/project/african1k/wildebeest/output/imputed/all.phased.maf0.05.r2.0.99.bcf -S black.list -r HiC_scaffold_${i} -Ov -o ${folder}/data/black.chr${i}.vcf 

	bcftools view /home/zilong/project/african1k/wildebeest/output/imputed/all.phased.maf0.05.r2.0.99.bcf -S Nyassa.list -r HiC_scaffold_${i} -Ov -o ${folder}/data/nyassa.chr${i}.vcf

	python3 run_loter.py ${folder}/data/wwhitebeard.chr${i}.vcf ${folder}/data/black.chr${i}.vcf ${folder}/data/nyassa.chr${i}.vcf ${i}.out 40
    fi
done


