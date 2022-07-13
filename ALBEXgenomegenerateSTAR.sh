#!/usr/bin/env bash
STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir /Volumes/archive/cancergeneticslab/ConorM_scratch_Nov21/ref_genomes/hg38ALBEX/ \
--genomeFastaFiles /Volumes/archive/cancergeneticslab/ConorM_scratch_Nov21/ref_genomes/hg38ALBEX/hg38.fa  \
--sjdbGTFfile /Volumes/archive/cancergeneticslab/ConorM_scratch_Nov21/ref_genomes/hg38ALBEX/hg38.ncbiRefSeq.gtf  \
--sjdbOverhang 100
