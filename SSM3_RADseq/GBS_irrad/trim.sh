#usr/bin/env bash
#Script to align reads to index
PATH2JAR=~/Trimmomatic-0.36/trimmomatic-0.36.jar
PATH2FQ=adaptertrimmed
mkdir trimmedfiles
for sample in `ls $PATH2FQ/*.fastq.gz`
do
dir=$PATH2FQ
x=$sample
y=${x%%.fastq.gz*}
z=${y##*/}
echo found $z
##single end
java -jar $PATH2JAR SE -phred33 $x trimmedfiles/$z.trimmed.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

done
#mv $dir/*.trimmed.fq /trimmedfiles
echo “files trimmed, moved to folder trimmedfiles within $PATH2FQ"
##paired end

#java -jar ~/Trimmomatic-0.36/trimmomatic-0.36.jar PE $BASE1 $BASE2 -baseout SSM3.fq.gz SRR2142076_1.fastq.gz SRR2142076_2.fastq.gz ILLUMINACLIP:TruSeq3-PE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
