
library(BSgenome.Mmusculus.UCSC.mm10)

library("biomaRt")
library(SomaticSignatures)
library(tidyverse)
library(reshape2)
library(factoextra)
library(stringr)
library(MutationalPatterns)
mygenome <- FaFile("/Volumes/archive/cancergeneticslab/ConorM_scratch_Nov21/ref_genomes/mouse2/mouse.fa")

SSM3genoms<-read.delim("AllSNPsnoheader.txt", sep =" ", header = FALSE)

names(SSM3genoms)<-c("chr", "start", "ref", "alt", "sample", "FMT")
SSM3genoms$freqs<-str_extract(SSM3genoms$FMT, pattern = "(?<=:)\\w{1,},\\w{1,}")
SSM3genoms<-SSM3genoms[-6]
SSM3genoms$RefFreq<-as.numeric(str_extract(SSM3genoms$freqs, pattern = "\\w{1,}(?=,)"))
SSM3genoms$AltFreq<-as.numeric(str_extract(SSM3genoms$freqs, "(?<=,)\\w{1,}"))
SSM3genoms$cov<-SSM3genoms$RefFreq + SSM3genoms$AltFreq
SSM3genoms$VAF<-SSM3genoms$AltFreq/SSM3genoms$cov


hist(SSM3genoms$VAF[SSM3genoms$VAF!=1])
hist(SSM3genoms$VAF[SSM3genoms$VAF<0.9])
SSM3genoms$end<-SSM3genoms$start
SSM3genoms$group<-NA
SSM3genoms$group[str_detect(SSM3genoms$sample, "^A")]<-"APOBEC'd"
SSM3genoms$group[str_detect(SSM3genoms$sample, "^E")]<-"E"
SSM3genoms.gr<-makeGRangesFromDataFrame(SSM3genoms, keep.extra.columns = TRUE)
table(SSM3genoms$sample)

SSM3_vr <- VRanges(
  seqnames = seqnames(SSM3genoms.gr),
  ranges = ranges(SSM3genoms.gr),
  ref = SSM3genoms.gr$ref,
  alt = SSM3genoms.gr$alt,
  sampleNames = SSM3genoms.gr$sample,
  group = SSM3genoms.gr$group)
SSM3_motifs<-mutationContext(SSM3_vr, mygenome)
SSM3_mm = motifMatrix(SSM3_motifs, group = "sampleNames", normalize = TRUE)
png("Rplots/96profilebysample.png")
plot_96_profile(SSM3_mm)
dev.off()

SSM3_mm_group = motifMatrix(SSM3_motifs, group = "group", normalize = TRUE)
png("Rplots/96profilebygroup.png")
plot_96_profile(SSM3_mm_group)
dev.off()
n_sigs=4
#' set the number of signatures to 4, and plot the signatures for both NMF and PCA
sigs_nmf = identifySignatures(SSM3_mm, n_sigs, nmfDecomposition)
plotSignatures(sigs_nmf) + ggtitle("NMF barchart")

plot_96_profile(SSM3_mm)

EditGroups_mm = motifMatrix(SSM3_motifs, group = "group", normalize = TRUE)
EditGroups_mm
plot_96_profile(EditGroups_mm)

##

SSM3uniquegenoms<-read.delim("uniqueSNPsnoheader.txt", sep =" ", header = FALSE)

names(SSM3uniquegenoms)<-c("chr", "start", "ref", "alt", "sample", "FMT")



SSM3uniquegenoms$end<-SSM3uniquegenoms$start
SSM3uniquegenoms$group<-NA
SSM3uniquegenoms$group[str_detect(SSM3uniquegenoms$sample, "^A")]<-"APOBEC'd"
SSM3uniquegenoms$group[str_detect(SSM3uniquegenoms$sample, "^E")]<-"E"
SSM3uniquegenoms<-SSM3uniquegenoms[which(SSM3uniquegenoms$sample != c("E1D12") & SSM3uniquegenoms$sample != c("A1A12")),]
SSM3uniquegenoms$sample<-as.character(SSM3uniquegenoms$sample)
SSM3uniquegenoms$freqs<-str_extract(SSM3uniquegenoms$FMT, pattern = "(?<=:)\\w{1,},\\w{1,}")
SSM3uniquegenoms<-SSM3uniquegenoms[-6]
SSM3uniquegenoms$RefFreq<-as.numeric(str_extract(SSM3uniquegenoms$freqs, pattern = "\\w{1,}(?=,)"))
SSM3uniquegenoms$AltFreq<-as.numeric(str_extract(SSM3uniquegenoms$freqs, "(?<=,)\\w{1,}"))
SSM3uniquegenoms$cov<-SSM3uniquegenoms$RefFreq + SSM3uniquegenoms$AltFreq
mean(SSM3uniquegenoms$cov)
SSM3uniquegenoms$VAF<-SSM3uniquegenoms$AltFreq/SSM3uniquegenoms$cov
SSM3uniquegenoms.gr<-makeGRangesFromDataFrame(SSM3uniquegenoms, keep.extra.columns = TRUE)




SSM3unique_vr <- VRanges(
  seqnames = seqnames(SSM3uniquegenoms.gr),
  ranges = ranges(SSM3uniquegenoms.gr),
  ref = SSM3uniquegenoms.gr$ref,
  alt = SSM3uniquegenoms.gr$alt,
  sampleNames = SSM3uniquegenoms.gr$sample,
  group = SSM3uniquegenoms.gr$group,
  VAF=SSM3uniquegenoms.gr$VAF)
SSM3unique_motifs<-mutationContext(SSM3unique_vr, mygenome)

SSM3unique_mm = motifMatrix(SSM3unique_motifs, group = "sampleNames", normalize = TRUE)
plot_96_profile(SSM3unique_mm)
n_sigs=8
#' set the number of signatures to 4, and plot the signatures for both NMF and PCA
sigs_nmf = identifySignatures(SSM3unique_mm, n_sigs, nmfDecomposition)
plotSignatures(sigs_nmf) + ggtitle("NMF barchart")
plotSamples(sigs_nmf)
cosmic<-read.csv("~/Desktop/COSMIC/sigProfiler_SBS_signatures_2018_03_28.csv", stringsAsFactors = FALSE)

cosmicsigs<-data.matrix(cosmic[,3:37])
plot_96_profile(cosmicsigs[c(2,3,17)])
plot_cosine_heatmap(cos_sim_matrix(data.matrix(SSM3unique_mm), cosmicsigs))
EditGroups_mm = motifMatrix(SSM3unique_motifs, group = "group", normalize = TRUE)
EditGroups_mm
plot_96_profile(EditGroups_mm)
depthinf<-read.delim("depthtab.txt")
depthinf<-depthinf[depthinf$sample %in% names(table(SSM3uniquegenoms$sample)),]
nmuts<-table(SSM3uniquegenoms$sample)

depthinf$nSNVs<-nmuts[match(data.frame(nmuts)$Var1, depthinf$sample)]
depthinf$mutrate<-depthinf$nSNVs/depthinf$greater20*1000000
depthinf$mutrate

depthinf
depthinf$group[str_detect(depthinf$sample, "^A")]<-"APOBEC'd"
depthinf$group[str_detect(depthinf$sample, "^E")]<-"control"
depthinf$mutrate
group[str_detect(SSM3uniquegenoms$sample, "^A")]<-"APOBEC'd"
png("Rplots/mutratesallSNVs.png")
ggplot(depthinf, aes(group, mutrate))  + geom_boxplot(aes(alpha=0.5)) + geom_jitter(aes(colour = sample)) + ylab("Mutation rate (muts/Mb, all mutations)")
dev.off()


###split into clonals
hist(SSM3uniquegenoms$cov[SSM3uniquegenoms$sample == "A1A2"], binwidth =0.5)
png("Rplots/VAFuniquevars.png", width = 1500)
ggplot(SSM3uniquegenoms[which(SSM3uniquegenoms$VAF!=1),], aes(VAF))+geom_density(aes(colour = sample), bw =0.05, alpha =0.05)+geom_density(bw=0.05) +
  +geom_rect()
  theme_classic()
dev.off()
?geom_rect()
png("Rplots/VAFuniquevarsnoHoms.png", width =150)
ggplot(SSM3uniquegenoms[SSM3uniquegenoms$VAF!=1,], aes(VAF)) + geom_density(bw=0.025) +theme_classic()
dev.off()

png("Rplots/VAFuniquevarsnoHoms.png")
ggplot(SSM3uniquegenoms, aes(VAF, colour = group)) + geom_density(bw=0.025) +theme_classic()
dev.off()

png("Rplots/VAFbygroup.png")
ggplot(SSM3uniquegenoms[SSM3uniquegenoms$VAF!=1,], aes(VAF, colour =group)) + geom_density(bw=0.025) +theme_classic()
dev.off()

png("Rplots/VAFbysample.png")
ggplot(SSM3uniquegenoms[SSM3uniquegenoms$VAF!=1,], aes(VAF, colour =sample)) + geom_density(bw=0.025) +theme_classic()
dev.off()

hets<-SSM3uniquegenoms[c(which(SSM3uniquegenoms$VAF<0.65 & SSM3uniquegenoms$VAF>0.35)),]
depthinf$nSNVshets<-table(hets$sample)
depthinf$mutratehets<-as.numeric(depthinf$nSNVshets)*10^6/depthinf$greater20
png("Rplots/mutratehets.png")
ggplot(depthinf, aes(group, mutratehets)) + geom_jitter(aes(colour = sample)) +
  geom_boxplot(alpha = 0.5) + ylab("Heterozygous mutation rate (muts/Mb)")
dev.off()
hist(hets$VAF)
hetsgenoms.gr<-makeGRangesFromDataFrame(hets, keep.extra.columns = TRUE)
hets_vr <- VRanges(
  seqnames = seqnames(hetsgenoms.gr),
  ranges = ranges(hetsgenoms.gr),
  ref = hetsgenoms.gr$ref,
  alt = hetsgenoms.gr$alt,
  sampleNames = hetsgenoms.gr$sample,
  group = hetsgenoms.gr$group)
hets_motifs<-mutationContext(hets_vr, mygenome)
hets_mm_sample = motifMatrix(hets_motifs, group = "sampleNames", normalize = TRUE)
png("Rplots/96profilebysamplehetsonly.png")
plot_96_profile(hets_mm_sample)
dev.off()



hets_mm = motifMatrix(hets_motifs, group = "group", normalize = TRUE)

png("Rplots/96profilebygrouphetsonly.png")
plot_96_profile(hets_mm)
dev.off()


hets_mm_sample_nn = motifMatrix(hets_motifs, group = "sampleNames", normalize = FALSE)
rownames(hets_mm_sample_nn)
A3Btargets<-colSums(hets_mm_sample_nn[c(13:16,29:32,45:48),])
depthinf$A3Btargetshets<-A3Btargets
depthinf$A3Btargethetsrate<-depthinf$A3Btargetshets/depthinf$greater20
depthinf
png("Rplots/A3btargethetrate.png")
ggplot(depthinf, aes(group, A3Btargethetsrate)) + geom_jitter(aes(colour = sample)) + geom_boxplot(alpha = 0.1)+ylab("Rate of APOBEC muts in TCN context (muts/Mb \n normalised to number of bp sequenced)")
plot_96_profile(hets_mm_sample_nn)
dev.off()

subclonal<-SSM3uniquegenoms[c(which(SSM3uniquegenoms$VAF<0.35)),]
depthinf$nSNVssubclonal<-table(subclonal$sample)
depthinf$mutratesubclonal<-as.numeric(depthinf$nSNVssubclonal)*10^6/depthinf$greater20
png("Rplots/mutratesubclonal.png")
ggplot(depthinf, aes(group, mutratesubclonal)) + geom_boxplot(alpha = 0.5)+ geom_jitter(aes(colour = sample)) +
    ylab("subclonal mutation rate (muts/Mb)")
dev.off()
hist(subclonal$VAF)
subclonalgenoms.gr<-makeGRangesFromDataFrame(subclonal, keep.extra.columns = TRUE)
subclonal_vr <- VRanges(
  seqnames = seqnames(subclonalgenoms.gr),
  ranges = ranges(subclonalgenoms.gr),
  ref = subclonalgenoms.gr$ref,
  alt = subclonalgenoms.gr$alt,
  sampleNames = subclonalgenoms.gr$sample,
  group = subclonalgenoms.gr$group)
subclonal_motifs<-mutationContext(subclonal_vr, mygenome)

subclonal_mm_sample_nn = motifMatrix(subclonal_motifs, group = "sampleNames", normalize = FALSE)
A3BtargetsSubclone<-colSums(subclonal_mm_sample_nn[c(29:32,45:48),])
depthinf$A3BtargetsSubclone<-A3BtargetsSubclone
depthinf$A3BtargetSubclonerate<-depthinf$A3BtargetsSubclone*1000000/depthinf$greater20
png("Rplots/A3btargetsubclonerate.png")

ggplot(depthinf, aes(group, A3BtargetSubclonerate)) + geom_jitter(aes(colour = sample)) + geom_boxplot(alpha = 0.1)+ylab("Rate of APOBEC muts in TCN context (muts/Mb)")
dev.off()
subclonal_mm = motifMatrix(subclonal_motifs, group ="group", normalize = FALSE)
png("Rplots/96profilesubclone.png")
plot_96_profile(subclonal_mm)
dev.off()


homs<-SSM3uniquegenoms[c(which(SSM3uniquegenoms$VAF>0.65)),]
depthinf$nSNVshoms<-table(homs$sample)
depthinf$mutratehoms<-as.numeric(depthinf$nSNVshoms)*10^6/depthinf$greater20
png("Rplots/mutratehoms.png")
ggplot(depthinf, aes(group, mutratehoms)) + geom_jitter(aes(colour = sample)) +
  geom_boxplot(alpha = 0.5) + ylab("homs mutation rate (muts/Mb)")
dev.off()
hist(homs$VAF)
homsgenoms.gr<-makeGRangesFromDataFrame(homs, keep.extra.columns = TRUE)
homs_vr <- VRanges(
  seqnames = seqnames(homsgenoms.gr),
  ranges = ranges(homsgenoms.gr),
  ref = homsgenoms.gr$ref,
  alt = homsgenoms.gr$alt,
  sampleNames = homsgenoms.gr$sample,
  group = homsgenoms.gr$group)
homs_motifs<-mutationContext(homs_vr, mygenome)

homs_mm_sample_nn = motifMatrix(homs_motifs, group = "sampleNames", normalize = FALSE)
A3Btargetshoms<-colSums(homs_mm_sample_nn[c(29:32,45:48),])
depthinf$A3Btargetshoms<-A3Btargetshoms
depthinf$A3Btargethomsrate<-depthinf$A3Btargetshoms/depthinf$greater20
ggplot(depthinf, aes(group, A3Btargethomsrate)) + geom_jitter(aes(colour = sample)) + geom_boxplot(alpha = 0.1)+ylab("Rate of APOBEC muts in TCN context (muts/Mb)")
homs_mm = motifMatrix(homs_motifs, group ="group", normalize = FALSE)
plot_96_profile(homs_mm)
plot_96_profile(hets_mm)
plot_96_profile(subclonal_mm)


allsub<-SSM3uniquegenoms[c(which(SSM3uniquegenoms$VAF<0.3)),]
depthinf$nSNVsallsub<-table(allsub$sample)
depthinf$mutrateallsub<-as.numeric(depthinf$nSNVsallsub)*10^6/depthinf$greater20
png("Rplots/mutrateallsub.png")
ggplot(depthinf, aes(group, mutrateallsub)) + geom_boxplot(alpha = 0.5)+geom_jitter(aes(colour = sample)) +
  ylab("allsub mutation rate (muts/Mb)")
dev.off()
hist(allsub$VAF)
ggplot(allsub,aes(VAF)) + geom_density(aes(colour = sample))
allsubgenoms.gr<-makeGRangesFromDataFrame(allsub, keep.extra.columns = TRUE)
allsub_vr <- VRanges(
  seqnames = seqnames(allsubgenoms.gr),
  ranges = ranges(allsubgenoms.gr),
  ref = allsubgenoms.gr$ref,
  alt = allsubgenoms.gr$alt,
  sampleNames = allsubgenoms.gr$sample,
  group = allsubgenoms.gr$group)
allsub_motifs<-mutationContext(allsub_vr, mygenome)
allsub_mm_sample_nn = motifMatrix(allsub_motifs, group = "sampleNames", normalize = FALSE)
A3Btargetsallsub<-colSums(allsub_mm_sample_nn[c(29:32,45:48),])
depthinf$A3Btargetsallsub<-A3Btargetsallsub
depthinf$A3Btargetallsubrate<-depthinf$A3Btargetsallsub/depthinf$greater20
ggplot(depthinf, aes(group, A3Btargetallsubrate)) + geom_jitter(aes(colour = sample)) + geom_boxplot(alpha = 0.1)+ylab("Rate of APOBEC muts in TCN context (muts/Mb)")
allsub_mm = motifMatrix(allsub_motifs, group ="group", normalize = FALSE)
plot_96_profile(allsub_mm)
plot_96_profile(hets_mm)
plot_96_profile(subclonal_mm)
depthinf
###split by group
SSM3groupgenoms<-read.delim("SNPsbygroupnoheader.txt", sep =" ", header = FALSE)

names(SSM3groupgenoms)<-c("chr", "start", "ref", "alt", "sample", "FMT")



SSM3groupgenoms$end<-SSM3groupgenoms$start
SSM3groupgenoms$group<-NA
SSM3groupgenoms$group[str_detect(SSM3groupgenoms$sample, "^A")]<-"APOBEC'd"
SSM3groupgenoms$group[str_detect(SSM3groupgenoms$sample, "^E")]<-"E"
SSM3groupgenoms<-SSM3groupgenoms[which(SSM3groupgenoms$sample != c("E1D12") & SSM3groupgenoms$sample != c("A1A12")),]
SSM3groupgenoms$sample<-as.character(SSM3groupgenoms$sample)
SSM3groupgenoms$freqs<-str_extract(SSM3groupgenoms$FMT, pattern = "(?<=:)\\w{1,},\\w{1,}")
SSM3groupgenoms<-SSM3groupgenoms[-6]
SSM3groupgenoms$RefFreq<-as.numeric(str_extract(SSM3groupgenoms$freqs, pattern = "\\w{1,}(?=,)"))
SSM3groupgenoms$AltFreq<-as.numeric(str_extract(SSM3groupgenoms$freqs, "(?<=,)\\w{1,}"))
SSM3groupgenoms$cov<-SSM3groupgenoms$RefFreq + SSM3groupgenoms$AltFreq
(SSM3groupgenoms$cov)
SSM3groupgenoms$VAF<-SSM3groupgenoms$AltFreq/SSM3groupgenoms$cov
hist(SSM3groupgenoms$VAF)
SSM3groupgenoms.gr<-makeGRangesFromDataFrame(SSM3groupgenoms, keep.extra.columns = TRUE)




SSM3group_vr <- VRanges(
  seqnames = seqnames(SSM3groupgenoms.gr),
  ranges = ranges(SSM3groupgenoms.gr),
  ref = SSM3groupgenoms.gr$ref,
  alt = SSM3groupgenoms.gr$alt,
  sampleNames = SSM3groupgenoms.gr$sample,
  group = SSM3groupgenoms.gr$group)
SSM3group_motifs<-mutationContext(SSM3group_vr, mygenome)
SSM3group_motifs
SSM3group_mm = motifMatrix(SSM3group_motifs, group = "sampleNames", normalize = FALSE)
plot_96_profile(SSM3group_mm)
A3Btargets<-colSums(SSM3group_mm[c(29:32,45:48),])
data.frame(A3Btargets, names(A3Btargets))
SSM3groupgenoms$VAF
n_sigs=8
#' set the number of signatures to 4, and plot the signatures for both NMF and PCA
sigs_nmf = identifySignatures(SSM3group_mm, n_sigs, nmfDecomposition)
plotSignatures(sigs_nmf) + ggtitle("NMF barchart")
plotSamples(sigs_nmf)
cosmic<-read.csv("~/Desktop/COSMIC/sigProfiler_SBS_signatures_2018_03_28.csv", stringsAsFactors = FALSE)

cosmicsigs<-data.matrix(cosmic[,3:37])
plot_96_profile(cosmicsigs[c(2,3,17)])
plot_cosine_heatmap(cos_sim_matrix(data.matrix(SSM3group_mm), cosmicsigs))
EditGroups_mm = motifMatrix(SSM3group_motifs, group = "group", normalize = TRUE)
EditGroups_mm
plot_96_profile(EditGroups_mm)
depthinf<-read.delim("depthtab.txt")
depthinf<-depthinf[depthinf$sample %in% names(table(SSM3groupgenoms$sample)),]
nmuts<-table(SSM3groupgenoms$sample)

depthinf$nSNVs<-nmuts[match(data.frame(nmuts)$Var1, depthinf$sample)]
depthinf$mutrate<-depthinf$nSNVs/depthinf$greater20*1000000
depthinf$mutrate

depthinf
depthinf$group[str_detect(depthinf$sample, "^A")]<-"APOBEC'd"
depthinf$group[str_detect(depthinf$sample, "^E")]<-"control"
depthinf$mutrate
group[str_detect(SSM3groupgenoms$sample, "^A")]<-"APOBEC'd"
png("Rplots/mutratesallSNVs.png")
ggplot(depthinf, aes(group, mutrate))  + geom_jitter() + geom_boxplot(aes(alpha=0.5)) + ylab("Mutation rate (muts/Mb, all mutations)")
dev.off()

SuM<-data.frame(SSM3unique_motifs)

ggplot(SuM[which(SuM$VAF<0.9),], aes(VAF)) + geom_density(aes(colour = sampleNames), bw = 0.05, alpha = 0.5) + geom_density(bw = 0.05) + theme_classic()
A3Btargets<-SuM[which(SuM$alteration == "CT" | SuM$alteration == "CG"),]
hist(SuM$VAF[SuM$sampleNames == "A2D3"])
hist(SuM$VAF[SuM$sampleNames == "E1D7"])
A3Btargets<-A3Btargets[str_detect(string = A3Btargets$context, pattern = "^T"),]
A3Btargets<-A3Btargets[which(A3Btargets$VAF !=1),]
A3Btargets
table(A3Btargets$sampleNames)/depthinf$greater20
ggplot(A3Btargets, aes(VAF)) + geom_density(aes(colour = sampleNames), stat = "count")
A3Bspl<-split(A3Btargets, A3Btargets$sampleNames)
hist(A3Bspl$A2D3$VAF, bw= 0.05)
normfactor<-depthinf$greater20
hist


SSM3indels<-read.delim("uniqueindelsnoheader.txt", sep =" ", header = FALSE)
names(SSM3indels)<-c("chr", "start", "ref", "alt", "sample", "FMT")
SSM3indels$end<-SSM3indels$start
SSM3indels$group<-NA
SSM3indels$group[str_detect(SSM3indels$sample, "^A")]<-"APOBEC'd"
SSM3indels$group[str_detect(SSM3indels$sample, "^E")]<-"E"
SSM3indels<-SSM3indels[which(SSM3indels$sample != c("E1D12") & SSM3indels$sample != c("A1A12")),]
SSM3indels$sample<-as.character(SSM3indels$sample)
SSM3indels$freqs<-str_extract(SSM3indels$FMT, pattern = "(?<=:)\\w{1,},\\w{1,}")
SSM3indels<-SSM3indels[-6]
SSM3indels$RefFreq<-as.numeric(str_extract(SSM3indels$freqs, pattern = "\\w{1,}(?=,)"))
SSM3indels$AltFreq<-as.numeric(str_extract(SSM3indels$freqs, "(?<=,)\\w{1,}"))
SSM3indels$cov<-SSM3indels$RefFreq + SSM3indels$AltFreq
mean(SSM3indels$cov)
SSM3indels$VAF<-SSM3indels$AltFreq/SSM3indels$cov
ggplot(SSM3indels[SSM3indels$VAF!=1,], aes(VAF)) + geom_density(aes(colour = sample))+ geom_density()
SSM3indels.gr<-makeGRangesFromDataFrame(SSM3indels, keep.extra.columns = TRUE)
nindels<-table(SSM3indels$sample)
nindels
depthinf$nindels<-nindels[match(data.frame(nindels)$Var1, depthinf$sample)]
depthinf$indelmutrate<-depthinf$nindels*1000000/depthinf$greater20
depthinf$indelmutrate
depthinf
nindelshet<-table(SSM3indels[SSM3indels$VAF<0.65 & SSM3indels$VAF> 0.35,]$sample)
depthinf$nindelshet<-nindelshet[match(data.frame(nindelshet)$Var1, depthinf$sample)]
depthinf$nindelshetmutrate<-depthinf$nindelshet*(10^6)/depthinf$greater20
ggplot(depthinf, aes(group, indelmutrate, colour = sample))+ geom_boxplot() + geom_point()
ggplot(depthinf, aes(group, nindelshetmutrate))+ geom_boxplot() + geom_jitter(aes(colour =sample))
