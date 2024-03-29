#!/bin/bash


PATH2VCFs=uniqueVCFhets
rm unSNPS.txt
mkdir indelsunique
mkdir SNPsunique
touch unSNPS.txt
for i in `ls $PATH2VCFs/*vcf`
do
  x=$i
  y=${x%%.vcf*}
  z=${y##*/}
  # echo $z yo
  bgzip -c $PATH2VCFs/$z.vcf > $PATH2VCFs/$z.vcf.gz
  tabix -p vcf $PATH2VCFs/$z.vcf.gz
  gatk SelectVariants -V $PATH2VCFs/$z.vcf -O SNPsunique/$z.vcf -select-type SNP --exclude-filtered
  gatk SelectVariants -V $PATH2VCFs/$z.vcf -O indelsunique/$z.vcf -select-type INDEL --exclude-filtered
  cat SNPsunique/$z.vcf | bcftools query -f "%CHROM %POS %REF %ALT $z %INFO/DP %FORMAT \n"  >> unSNPS.txt
done
 rm -r uniqueVCFhets/*gz*

PATH2VCFs=totalVCFnoGerm/
rm totSNPS.txt
mkdir indelstotal
mkdir SNPstotal
touch totSNPS.txt
for i in `ls $PATH2VCFs/*vcf`
do
  x=$i
  y=${x%%.vcf*}
  z=${y##*/}
  # echo $z yo
  bgzip -c $PATH2VCFs/$z.vcf > $PATH2VCFs/$z.vcf.gz
  tabix -p vcf $PATH2VCFs/$z.vcf.gz
  gatk SelectVariants -V $PATH2VCFs/$z.vcf -O SNPstotal/$z.vcf -select-type SNP --exclude-filtered
  gatk SelectVariants -V $PATH2VCFs/$z.vcf -O indelstotal/$z.vcf -select-type INDEL --exclude-filtered
  cat SNPstotal/$z.vcf | bcftools query -f "%CHROM %POS %REF %ALT $z %INFO/DP %FORMAT \n"  >> totSNPS.txt
done
 rm -r totalVCFhets/*gz*
