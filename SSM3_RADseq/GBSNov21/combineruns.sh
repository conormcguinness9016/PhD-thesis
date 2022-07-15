#usr/bin/env bash
mkdir GBS_combined
mkdir GBS_combined/adaptertrimmed
for i in $(cat nov21_GBS/MyAnalysis/samples.txt)
do
 cat nov21_GBS/MyAnalysis/adaptertrimmed/$i.fastq.gz dec21_GBS/MyAnalysis/adaptertrimmed/$i.fastq.gz > GBS_combined/adaptertrimmed/$i.fastq.gz
 done
