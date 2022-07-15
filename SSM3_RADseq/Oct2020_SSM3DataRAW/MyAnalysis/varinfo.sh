#usr/bin/env bash
rm -r totalVCFhets
rm -r uniqueVCFhets
mkdir totalVCFhets
mkdir uniqueVCFhets
cp /Volumes/archive/cancergeneticslab/ConorM/GBS_irrad/totalVCFhets/0D5.trimmed.vcf .
bgzip 0D5.trimmed.vcf
bcftools index 0D5.trimmed.vcf.gz
cd totalVCFs
for i in `ls *vcf`
do
  x=${i%totals*}
  y=${x##*BQSR_}
  echo $x $y
  bcftools view -i "FMT/AF<0.75 & INFO/DP>20 & FMT/AF>0.30" $i > ../totalVCFhets/$y.vcf
  bgzip ../totalVCFhets/$y.vcf
  bcftools index ../totalVCFhets/$y.vcf.gz
 # bcftools filter -R /Volumes/scratch/cancergeneticslab/SSM3/totalVCFhets/SRR2142076.vcf.gz ../totalVCFhets/$y.vcf.gz > tmp.vcf
 # bgzip tmp.vcf
 # bcftools index tmp.vcf.gz
 # bcftools filter -R ../0D5.trimmed.vcf.gz tmp.vcf.gz > ../uniqueVCFhets/$y.vcf
 bcftools isec -w 1 -Ov -n~100 ../totalVCFhets/$y.vcf.gz ../0D5.trimmed.vcf.gz /Volumes/scratch/cancergeneticslab/SSM3/totalVCFhets/SRR2142076.vcf.gz > ../uniqueVCFhets/$y.vcf

 gunzip ../totalVCFhets/$y.vcf.gz
 gunzip ../uniqueVCFhets/$y.vcf.gz
  #rm tmp.vcf.gz*
done
rm ../totalVCFhets/*gz*
rm ../uniqueVCFhets/*gz*