#usr/bin/env bash
GERMVCF=/Volumes/scratch/cancergeneticslab/ConorM_scratch_Nov21/ref_genomes/mouse2/00-All.vcf.gz
rm -r totalVCFhets
rm -r uniqueVCFhets
rm -r totalVCFnoGerm
mkdir totalVCFhets
mkdir uniqueVCFhets
rm depthtab.txt
touch depthtab.txt
grep '^"sample"' /Volumes/archive/cancergeneticslab/ConorM/GBS_irrad/depthtab.txt | awk -F ':' 'BEGIN{OFS="\t";} {print $1, $2, $3, $4 }' >> depthtab.txt
grep -v '^"sample"' /Volumes/archive/cancergeneticslab/ConorM/GBS_irrad/depthtab.txt | awk -F ':' 'BEGIN{OFS="\t";} {gsub("\"", "", $1); print "irrad_"$1"", $2, $3, $4 }' >> depthtab.txt
grep -v '^"sample"' /Volumes/archive/cancergeneticslab/ConorM/GBSNov21/GBS_combined/depthtab.txt | awk -F ':' 'BEGIN{OFS="\t";} {gsub("\"", "", $1); print "SA_"$1"", $2, $3, $4 }' >> depthtab.txt
grep -v '^"sample"' /Volumes/archive/cancergeneticslab/ConorM/Oct2020_SSM3DataRAW/MyAnalysis/depthtab.txt | awk -F ':' 'BEGIN{OFS="\t";} {gsub("\"", "", $1); print "PIG_"$1"", $2, $3, $4 }' >> depthtab.txt
for file in /Volumes/archive/cancergeneticslab/ConorM/GBS_irrad/totalVCFhets/*vcf
do
  x=${file%trimmed*}
  cp $file totalVCFhets/irrad_${x##*hets/}vcf
done
for file in /Volumes/archive/cancergeneticslab/ConorM/GBSNov21/GBS_combined/totalVCFhets/*vcf
do
  x=${file%vcf*}
  cp $file totalVCFhets/SA_${x##*hets/}vcf
done

for file in /Volumes/archive/cancergeneticslab/ConorM/Oct2020_SSM3DataRAW/MyAnalysis/totalVCFhets/*vcf
do
  x=${file%vcf*}
  cp $file totalVCFhets/PIG_${x##*hets/}vcf
done
cd totalVCFhets
bgzip irrad_0D5.vcf
bcftools index irrad_0D5.vcf.gz
bgzip SA_OD5.vcf
bcftools index SA_OD5.vcf.gz
mkdir ../totalVCFnoGerm/
for i in `ls *vcf*`
do
  x=${i%.vcf*}
  bgzip $x.vcf
  bcftools index $x.vcf.gz
  echo y is $y
  echo  x is $x
  bcftools isec -w 1 -Ov -n~10 $x.vcf.gz $GERMVCF > ../totalVCFnoGerm/$x.vcf
  bcftools isec -w 1 -Ov -n~1000 $x.vcf.gz irrad_0D5.vcf.gz SA_OD5.vcf.gz /Volumes/scratch/cancergeneticslab/SSM3/totalVCFhets/SRR2142076.vcf.gz > ../uniqueVCFhets/$x.vcf
done
cd ..
gunzip uniqueVCFhets/*gz
gunzip totalVCFhets/*gz
rm uniqueVCFhets/*gz*
rm totalVCFhets/*gz*

