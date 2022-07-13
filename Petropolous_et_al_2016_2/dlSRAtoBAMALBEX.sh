#!/usr/bin/env bash

#Script to align reads to index
rm bwalog.txt
PATH2FQ=fastq
rm -r BAMfilesALBEX
mkdir fastq
cd $PATH2FQ
mkdir ../BAMfilesALBEX
mkdir ../BAMfilesALBEX/STAR
mkdir ../transcripts
mkdir ../A3Bintrons
STAR --genomeLoad LoadAndExit --genomeDir /Volumes/scratch/cancergeneticslab/ConorM_scratch_Nov21/ref_genomes/hg38ALBEX/
for sample in `cat ../SRR_Acc_list`
do
  x=$sample
  y=${x%.*}
  z=${y##*/}
  p=${z%.*}
  if [[ -f ../BAMfilesALBEX/STAR/${p}Aligned.sortedByCoord.out.bam ]]
  then
    echo skipping $x as already processed
  else
    echo downloading and aligning $z
    #fasterq-dump $sample
    mkdir ../BAMfilesALBEX/STAR/$p
    date +"%T" >> ../bwalog.txt
    echo processing $z with BWA >> ../bwalog.txt
    pigz -f ${x}.fastq
    STAR --genomeLoad LoadAndKeep \
          --limitBAMsortRAM 10000000000 \
           --outFilterScoreMinOverLread 0.1 \
           --outFilterMatchNminOverLread 0.1 \
           --runThreadN 4 \
           --genomeDir /Volumes/scratch/cancergeneticslab/ConorM_scratch_Nov21/ref_genomes/hg38ALBEX/ \
           --readFilesIn $p.fastq \
           --readFilesCommand zcat \
           --quantMode GeneCounts \
           --outSAMstrandField intronMotif \
           --outSAMtype BAM SortedByCoordinate  \
           --outFileNamePrefix ../BAMfilesALBEX/STAR/$p
    date +"%T" >> ../bwalog.txt
    echo BAMfile generated for $p >> ../bwalog.txt
    pwd
    BAM=../BAMfilesALBEX/STAR/${p}Aligned.sortedByCoord.out.bam
    samtools index -@ 3 $BAM
    date +"%T" >> ../bwalog.txt
    echo bam file generated for $z >> ../bwalog.txt
  fi
  if [[ -d ../transcripts/$z/ ]]
  then
    echo skipping $x as albex1 already quantified
  else
    echo $p
    # ~/cufflinks/cufflinks -G ../../ALBEX.gtf $BAM -o ../transcripts/$z/
    # awk -v var=$z ' NR>1 { print $0, var } ' OFS='\t' ../transcripts/$z/isoforms.fpkm_tracking >> ../isoforms.fpkm_tracking.1
  fi
  if [ -d ./A3Bintrons/$p/ ]
  then
    echo skipping $x as already processed
  else
    echo $p
    # ~/cufflinks/cufflinks -G ../../a3bintrons.gtf $BAM -o ../A3Bintrons/$p/
    # awk -v var=$p ' NR>1 { print $0, var } ' OFS='\t' ../A3Bintrons/$p/isoforms.fpkm_tracking >> ../introns.fpkm_tracking.1
  fi
done
echo “SAM files created, moved to folder SAMfiles within $PATH2FQ.”
