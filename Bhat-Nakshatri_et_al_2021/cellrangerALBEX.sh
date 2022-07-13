#!/usr/bin/env bash
for sample in `cat SRR_Acc_list`
do
  ~/cellranger-6.1.1/cellranger count --nosecondary --id=ALBEX1_$sample \
  --fastqs=outs/$sample \
  --transcriptome=/Network/Servers/biocldap.otago.ac.nz/Volumes/userdata/student_users/conormcguinness/refdata-gex-GRCh38-2020-ALBEX1
done
