
#usr/bin/env bash
cd totalVCFhets
rm *gz
for i in S*vcf
do
echo $i
  bgzip $i
  tabix $i.gz
done

for i in P*vcf
do
echo $i
  bgzip $i
  tabix $i.gz
done
gunzip *DT*_DT*gz
rm *NEG*
bcftools merge *gz > ../merge.vcf
gunzip *gz
rm *tbi
cd ..

 bcftools +setGT merge.vcf -Ov -o filled.vcf -- -t . -n 0
 git clone https://github.com/edgardomortiz/vcf2phylip
 python3 vcf2phylip/vcf2phylip.py -i filled.vcf -o OD5
