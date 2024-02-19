#!/bin/bash
beagle3=/home/zilong/local/bin/beagle3.jar

VCF=/home/users/xiaodong/Documents/Project/Wildebeest/dataset/dataset1/Wildebeest_wildebeest_variable_sites_nomultiallelics_noindels_maf005.bcf.gz
OUT=/home/zilong/project/african1k/wildebeest/output
mkdir -p $OUT/vcf $OUT/imputed

runBeagle3() {
    # show verbose log
    echo "start running beagle3 by chr"
    chrs=$(bcftools index -s $VCF|cut -f1|tail -n+$1 | head -n $2)   # change chr name patterns to your own
    for chrom in $chrs;do
    {
        # convert PL tag to GL in beagle format 
        indir=$OUT/vcf
        bname=$indir/`basename $VCF`.$chrom
        # bcftools view -r $chrom -O u $VCF | bcftools +tag2tag -Ov -o ${bname}.vcf -- -r --pl-to-gl && vcftools --vcf ${bname}.vcf --out ${bname} --BEAGLE-GL --chr $chrom \
        #   && rm -f ${bname}.vcf && echo 'convert to BEAGLE-GL done'

        # run imputation
        outdir=$OUT/imputed
        out=$outdir/$chrom
        input=${bname}.BEAGLE.GL
        # java -Xss5m -Xmx20g -Djava.io.tmpdir=$outdir -jar $beagle3 like=$input omitprefix=true out=$out
        zcat $out.BEAGLE.GL.gprobs.gz | java -jar /home/zilong/local/bin/gprobs2beagle.jar 0.9 -1 | gzip -c >$out.bgl.gz && \
        zcat $out.BEAGLE.GL.gprobs.gz | awk 'NR>1{split($1,a,":");print $1,a[2],$2,$3}' >$out.bgl.sites && \
        java -jar /home/zilong/local/bin/beagle2vcf.jar $chrom $out.bgl.sites $out.bgl.gz -1 |bgzip -c >$out.vcf.gz && \
        bcftools index -f $out.vcf.gz && calc_imputed_gt_discord.py -chr $chrom $VCF $out.vcf.gz $out.sum.disc.count && ech "convert to vcf done"

        ### convert phased.gz to hap then to vcf
        beagle_phased_to_hap_sample.py $chrom $out.BEAGLE.GL.phased.gz $out.bgl.sites $out && bcftools convert --hapsample2vcf $out |bcftools annotate -I +'%CHROM:%POS' -Oz -o $out.phased.vcf.gz && bcftools index -f $out.phased.vcf.gz
        echo "$chrom done"
    } &
    done
    wait

    echo "all jobs done by chroms"


    echo "beagle3 imputation done"
}

chroms=$(bcftools index -s $VCF|cut -f1)      # not array just string sep by space
nchr=$(echo $chroms|wc -w)
n=6
# for i in $(seq 1 $n $nchr)
# do
# 	runBeagle3 $i $n
# done

echo "start concating files"
outdir=$OUT/imputed
out=$OUT/imputed/all
bcftools concat --threads 10 -Ob -o $out.bcf `for i in $chroms; do echo $outdir/$i.vcf.gz;done` && bcftools index -f $out.bcf &
bcftools concat --threads 10 `for i in $chroms; do echo $outdir/$i.phased.vcf.gz;done` | bcftools annotate -I +'%CHROM:%POS' -Ob -o $out.phased.bcf --threads 10  && bcftools index -f $out.phased.bcf &
for i in $chroms;do
    cat $outdir/$i.BEAGLE.GL.r2
done >$out.r2
echo "done concating files"
wait

# run post QC

echo "run post QC"

vcfin=$outdir/all.bcf
vcfphased=$outdir/all.phased.bcf
out=$outdir/all
stats=$out.imputed.sites.stats
bcftools +fill-tags $vcfin -- -t MAF |bcftools query -f "%ID\t%MAF\n" >$out.maf
# # merge r2 and maf together
echo -e "ID\tR2\tMAF" >$stats
paste $out.r2 <(awk '{print $2}' $out.maf) >>$stats
plot-r2-vs-maf.R $stats $stats && echo plot done
# apply filters MAF>0.05 && R2 >0.95
# remove NaN sites in R2 file
awk 'NR>1 && $3>0.05 && $2>0.99 && $2!="NaN"{print $1}' $stats >$out.maf0.05.r2.0.99.sites
bcftools view -i ID=@$out.maf0.05.r2.0.99.sites -Ob -o $out.maf0.05.r2.0.99.bcf --threads 10 $vcfin && bcftools index -f $out.maf0.05.r2.0.99.bcf
# # for phased vcf
bcftools view -i ID=@$out.maf0.05.r2.0.99.sites -Ob -o $out.phased.maf0.05.r2.0.99.bcf --threads 10 $vcfphased && bcftools index -f $out.phased.maf0.05.r2.0.99.bcf

echo "done post QC"
