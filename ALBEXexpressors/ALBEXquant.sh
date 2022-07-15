#!/usr/bin/env bash
mkdir transcripts
mkdir STAR
mkdir fastq
mkdir A3BAMS
mkdir A3FQS
mkdir A3FAS
echo "date +"%T"" 1 thread start >> time.txt
STAR --genomeLoad LoadAndExit --genomeDir /Volumes/scratch/cancergeneticslab/ConorM/ref_genomes/hg38/
for sample in `cat ALBEXexpressors.txt`
do
  x=$sample
  y=${x%.*}
  z=${y##*/STAR/}
  p=${z%/*}
  d=${sample%BAMfiles/*}fastq/
  f=${sample%BAMfiles/*}fastq/$p
  paired=${f}_1.fastq.gz
  single=${f}.fastq.gz
  pairedunzip=${f}_1.fastq
  singleunzip=${f}.fastq
  # echo $single
  # echo z is $z
  # echo p is $p
  # echo d is $d
  echo running ALBEXquant.sh
  # if [[ -e $single ]] || [[ -e $paired ]] || [[ -e $pairedunzip ]] || [[ -e $singleunzip ]]
  # then
  #   echo "fastq present ($d/${p})"
  #   mkdir STAR/$p
  #   ls $d/${p}*
  #   #$d/${p}*
  #   [ -e $single ] &&  STAR --genomeLoad LoadAndKeep \
  #       --limitBAMsortRAM 10000000000 \
  #        --outFilterScoreMinOverLread 0.1 \
  #        --outFilterMatchNminOverLread 0.1 \
  #        --runThreadN 4 \
  #        --genomeDir /Volumes/scratch/cancergeneticslab/ConorM/ref_genomes/hg38/ \
  #        --readFilesIn $single \
  #        --readFilesCommand zcat \
  #        --quantMode GeneCounts \
  #        --outSAMstrandField intronMotif \
  #        --outSAMtype BAM SortedByCoordinate  \
  #        --outFileNamePrefix STAR/$p/
  # else
  #   cd fastq
  #   [ ! -e ${p}.fastq ] && fasterq-dump $p
  #   cd ..
  #   mkdir STAR/$p
  #   STAR --genomeLoad LoadAndKeep \
  #       --limitBAMsortRAM 10000000000 \
  #        --outFilterScoreMinOverLread 0.1 \
  #        --outFilterMatchNminOverLread 0.1 \
  #        --runThreadN 4 \
  #        --genomeDir /Volumes/scratch/cancergeneticslab/ConorM/ref_genomes/hg38/ \
  #        --readFilesIn fastq/$p.fastq \
  #        --quantMode GeneCounts \
  #        --outSAMstrandField intronMotif \
  #        --outSAMtype BAM SortedByCoordinate  \
  #        --outFileNamePrefix STAR/$p/
  # fi
  # echo indexing STAR/$p/Aligned.sortedByCoord.out.bam
  # samtools index STAR/$p/Aligned.sortedByCoord.out.bam
  # echo subsetting bam to apobec region
  # samtools view -h -b STAR/$p/Aligned.sortedByCoord.out.bam "chr22:38957504-39000467" > A3BAMS/$p.bam
  # echo indexing subset bam
  samtools reheader STAR/$p/Aligned.sortedByCoord.out.bam A3BAMS/$p.bam > tst_$p.bam
  mv tst_$p.bam A3BAMS/$p.bam
  samtools index A3BAMS/$p.bam
  # # ~/cufflinks/cufflinks A3BAMS/$p.bam \
  # #   -o transcripts/$p/ \
  # #   --library-type fr-unstranded
  # echo converting bam to fastq
  # bedtools bamtofastq -i A3BAMS/$p.bam -fq A3FQS/$p.fq
done
cd A3BAMS
find . -type f -empty -print -delete
ls -1 ERR104*bam > bamlist.txt
ls -1 SRR4*bam >> bamlist.txt
samtools merge -b bamlist.txt merge.bam
samtools sort merge.bam -o merge.sorted.bam
samtools index merge.sorted.bam
Trinity --genome_guided_bam merge.sorted.bam --max_memory 50G --genome_guided_max_intron 50000 --CPU 10 --no_salmon
cp trinity_out_dir ..
