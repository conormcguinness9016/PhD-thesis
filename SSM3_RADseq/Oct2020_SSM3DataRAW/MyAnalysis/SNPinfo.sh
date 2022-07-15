#usr/bin/env bash
rm uniqueSNPsnoheader.txt
rm AllSNPsnoheader.txt
rm SNPsbygroupnoheader.txt
rm hetSNPsnoheader.txt
mkdir SNPshets
cd SNPs/uniqueSNPs
for i in `ls`
do
  grep -v ^# $i | awk -v myvar=${i%%.vcf*} '{print $1, $2, $4, $5, myvar, $10 }' >> ../../uniqueSNPsnoheader.txt
done
cd ../../
cd SNPs/tmp
for i in `ls *gz`
do
  bcftools view -i "INFO/DP>20" $i | grep -v ^# | awk -v myvar=${i%%SNPS.vcf.gz*} '{print $1, $2, $4, $5, myvar, $10 }' >> ../../AllSNPsnoheader.txt
   bcftools view -i "FMT/AF<0.75 & FMT/AF > 0.2 & INFO/DP>20" $i | grep -v ^# | awk -v myvar=${i%%SNPS.vcf.gz*} '{print $1, $2, $4, $5, myvar, $10 }' >> ../../hetSNPsnoheader.txt
   bcftools view -i "FMT/AF<0.75 & FMT/AF > 0.2 & INFO/DP>20" $i > ../../SNPshets/${i%SNPS.vcf*}.vcf
done

cd ../SNPsbygroup

for i in `ls`
do
  grep -v ^# $i | awk -v myvar=${i%%.vcf*} '{print $1, $2, $4, $5, myvar, $10 }' >> ../../SNPsbygroupnoheader.txt
done

cd ../../
