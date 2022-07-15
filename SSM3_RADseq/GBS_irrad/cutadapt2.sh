#usr/bin/env bash
mkdir cutadapt2
mkdir cutadapt2/lane1
mkdir adaptertrimmed
cutadapt --no-indels -e 0 -g file:bc.fasta -o "cutadapt2/lane1/{name}.fastq.gz" SQ0818_CD04EANXX_s_5.fastq
for i in `cat samples.txt`
do
  cat cutadapt2/lane1/$i*fastq.gz > adaptertrimmed/$i.fastq.gz
  done
