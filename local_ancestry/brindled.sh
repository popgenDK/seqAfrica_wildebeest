#!/bin/bash


for i in {1..29}
do
    if [[ "$i" -eq 3 ]]
    then
        :
    else
	bcftools view /home/zilong/project/african1k/wildebeest/output/imputed/all.phased.maf0.05.r2.0.99.bcf -S wwhitebeard.list -r HiC_scaffold_${i} -Ov -o wwhitebeard.chr${i}.vcf

	bcftools view /home/zilong/project/african1k/wildebeest/output/imputed/all.phased.maf0.05.r2.0.99.bcf -S black.list -r HiC_scaffold_${i} -Ov -o black.chr${i}.vcf 

	bcftools view /home/zilong/project/african1k/wildebeest/output/imputed/all.phased.maf0.05.r2.0.99.bcf -S brindled4.list -r HiC_scaffold_${i} -Ov -o brindled4.chr${i}.vcf

	python run_loter.py wwhitebeard.chr${i}.vcf black.chr${i}.vcf brindled4.chr${i}.vcf ${i}.out 40
    fi
done


