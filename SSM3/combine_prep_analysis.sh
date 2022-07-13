#usr/bin/env bash
rm -r germremoved_SA_SigProfileranalysis*
mkdir germremoved_SA_SigProfileranalysis
cp totalVCFnoGerm/*vcf germremoved_SA_SigProfileranalysis
#cp totalVCFnoGerm/irrad_0*vcf germremoved_SA_SigProfileranalysis
cp /Volumes/scratch/cancergeneticslab/SSM3/SNPs/tmp/germremovedSSM3.vcf .
bcftools view -i 'FORMAT/DP>20 & FORMAT/AF<0.7 & FORMAT/AF>0.25' germremovedSSM3.vcf > germremoved_SA_SigProfileranalysis/SSM3_clonal.vcf
bcftools view -i 'FORMAT/DP>20 & FORMAT/AF<=0.25' germremovedSSM3.vcf > germremoved_SA_SigProfileranalysis/SSM3_subclonal.vcf
mkdir germremoved_SA_SigProfileranalysis_Hets_sep
# cp totalVCFnoGerm/irrad_0*vcf germremoved_SA_SigProfileranalysis_Hets_sep/
cp  germremoved_SA_SigProfileranalysis/SSM3* germremoved_SA_SigProfileranalysis_Hets_sep/
for i in `ls totalVCFnoGerm/*vcf`
do
bcftools view -i 'FORMAT/DP>20 & FORMAT/AF<0.7 & FORMAT/AF>0.2' $i > germremoved_SA_SigProfileranalysis_Hets_sep/clonal_${i##*/}
bcftools view -i 'FORMAT/DP>20 & FORMAT/AF<=0.2' $i > germremoved_SA_SigProfileranalysis_Hets_sep/subclonal_${i##*/}
done
# python3 SigProfiler.py
