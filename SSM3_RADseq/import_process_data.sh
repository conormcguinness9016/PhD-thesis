#usr/bin/env bash
GERMVCF=/Volumes/scratch/cancergeneticslab/ConorM_scratch_Nov21/ref_genomes/mouse2/00-All.vcf.gz
rm -r totalVCFs
rm -r uniqueVCFtotals
rm -r totalVCFnoGerm
mkdir totalVCFs
mkdir uniqueVCFtotals
rm depthtab.txt
touch depthtab.txt
grep '^"sample"' /Volumes/archive/cancergeneticslab/ConorM/SSM3_RADseq/GBS_irrad/depthtab.txt | awk -F ':' 'BEGIN{OFS="\t";} {print $1, $2, $3, $4 }' >> depthtab.txt
grep -v '^"sample"' /Volumes/archive/cancergeneticslab/ConorM/SSM3_RADseq/GBS_irrad/depthtab.txt | awk -F ':' 'BEGIN{OFS="\t";} {gsub("\"", "", $1); print "irrad_"$1"", $2, $3, $4 }' >> depthtab.txt
grep -v '^"sample"' /Volumes/archive/cancergeneticslab/ConorM/SSM3_RADseq/GBSNov21/GBS_combined/depthtab.txt | awk -F ':' 'BEGIN{OFS="\t";} {gsub("\"", "", $1); print "SA_"$1"", $2, $3, $4 }' >> depthtab.txt
grep -v '^"sample"' /Volumes/archive/cancergeneticslab/ConorM/SSM3_RADseq/Oct2020_SSM3DataRAW/MyAnalysis/depthtab.txt | awk -F ':' 'BEGIN{OFS="\t";} {gsub("\"", "", $1); print "PIG_"$1"", $2, $3, $4 }' >> depthtab.txt
for file in /Volumes/archive/cancergeneticslab/ConorM/SSM3_RADseq/GBS_irrad/totalVCFs/*vcf
do
  x=${file%.trimmed*}
  echo x = $x
   echo ${x##*BQSR_}
  cp $file totalVCFs/irrad_${x##*BQSR_}.vcf
done
for file in /Volumes/archive/cancergeneticslab/ConorM/SSM3_RADseq/GBSNov21/GBS_combined/totalVCFs/*vcf
do
  x=${file%totals.vcf*}
  cp $file totalVCFs/SA_${x##*BQSR_}.vcf
done
# 
for file in /Volumes/archive/cancergeneticslab/ConorM/SSM3_RADseq/Oct2020_SSM3DataRAW/MyAnalysis/totalVCFs/*vcf
do
  x=${file%totals.vcf*}
  cp $file totalVCFs/PIG_${x##*BQSR_}.vcf
done
cd totalVCFs
bgzip -@ 8 irrad_0D5.vcf
bcftools index irrad_0D5.vcf.gz
bgzip -@ 8 SA_OD5.vcf
bcftools index SA_OD5.vcf.gz
mkdir ../totalVCFnoGerm/
for i in `ls *vcf*`
do
  x=${i%.vcf*}
  bgzip -@ 8 $x.vcf
  bcftools index $x.vcf.gz
  echo y is $y
  echo  x is $x
  bcftools isec -w 1 -Ov -n~10 $x.vcf.gz $GERMVCF > ../totalVCFnoGerm/$x.vcf
  bcftools isec -w 1 -Ov -n~1000 $x.vcf.gz irrad_0D5.vcf.gz SA_OD5.vcf.gz /Volumes/scratch/cancergeneticslab/SSM3/totalVCFs/SRR2142076.vcf.gz > ../uniqueVCFtotals/$x.vcf
done
cd ..
gunzip uniqueVCFtotals/*gz
gunzip totalVCFs/*gz
rm uniqueVCFtotals/*gz*
rm totalVCFs/*gz*
