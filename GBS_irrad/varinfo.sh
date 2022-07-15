#usr/bin/env bash
rm -r totalVCFhets
mkdir totalVCFhets
mkdir uniqueVCFhets
cd totalVCFs
for i in `ls *vcf`
do
  x=${i%totals*}
  y=${x##*BQSR_}
  echo $x $y
  bcftools view -i "FMT/AF<0.75 & INFO/DP>20" $i > ../totalVCFhets/$y.vcf
  bgzip ../totalVCFhets/$y.vcf
  bcftools index ../totalVCFhets/$y.vcf.gz
  #bcftools filter -R /Volumes/scratch/cancergeneticslab/SSM3/totalVCFhets/SRR2142076.vcf.gz ../totalVCFhets/$y.vcf.gz > ../uniqueVCFhets/$y.vcf
  bcftools isec -w 1 -n~10 ../totalVCFhets/$y.vcf.gz /Volumes/scratch/cancergeneticslab/SSM3/totalVCFhets/SRR2142076.vcf.gz > ../uniqueVCFhets/$y.vcf
  gunzip ../totalVCFhets/$y.vcf
done
rm ../totalVCFhets/*gz*
