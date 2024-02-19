plink1.9/plink --keep ${subspecies.id.list} --bfile ${bfile} --allow-extra-chr --chr-set 29 --make-bed --out bfile/${dirname}/${subspecies}
plink1.9/plink --genome full --geno 0.05 --maf 0.05  --bfile bfile/${dirname}/${subspecies} --allow-extra-chr --out ./output/${dirname}/plink/${subspecies}.relaedness --ppc-gap 0 --chr-set 29
plink_2/plink2 --bfile bfile/${dirname}/${subspecies} --geno 0.05  --maf 0.05  --allow-extra-chr --make-king-table --out ./output/${dirname}/king/${subspecies}.king --chr-set 29  #--king-table-filter 0.05
plink1.9/plink --bfile bfile/${dirname}/${subspecies} --geno 0.05 --maf 0.05 --allow-extra-chr --het small-sample --out ./output/${dirname}/inbreeding_coef/${subspecies}.inbreeding_coef --chr-set 29
