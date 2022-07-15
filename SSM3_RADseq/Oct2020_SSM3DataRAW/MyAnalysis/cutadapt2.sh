#usr/bin/env bash
mkdir cutadapt2
mkdir cutadapt2/lane1
mkdir cutadapt2/lane2
mkdir adaptertrimmed
cutadapt --no-indels -e 0 -g file:bc.fasta -o "cutadapt2/{name}.fastq.gz" SQ1429_CECCRANXX_s_4.fastq
# for i in `cat samples.txt`
# do
#   cat cutadapt2/lane1/$i*fastq.gz > cutadapt2/lane1_$i.fastq.gz
#   cat cutadapt2/lane2/$i*fastq.gz > cutadapt2/lane2_$i.fastq.gz
#   cat cutadapt2/lane*_$i.fastq.gz > adaptertrimmed/$i.fastq.gz
#   done
