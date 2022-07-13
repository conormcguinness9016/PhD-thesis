library(Biostrings)
library("BSgenome.Hsapiens.UCSC.hg38")
library(plyr)
library(ggplot2)
library(reshape2)
library(scales)
library(tidyverse)
Hsapiens[[1]]

simseq <- sim.DNAseq(size=100000, GCfreq=0.51)
simseq<-DNAString(simseq)
RElist$Enzyme<-as.character(RElist$Enzyme)
REfrags<-list()
revcomp<-list()
REfraglibs<-list()
for(i in 1:nrow(RElist)){
  m<-matchPattern(RElist$Seq[i], simseq, fixed = FALSE)
  ends<-end(gaps(m))
  temp_df<-data.frame(end=ends + RElist$FiveCutSite[i] + 1,chr="simseq")
  temp_df$start<-c(temp_df$end[1:nrow(temp_df)] - 1)
  if(RElist$Seq[i] != RElist$revcompSeq[i]){
    m<-matchPattern(RElist$revcompSeq[i], simseq, fixed = FALSE)
    ends<-end(gaps(m))
    revtemp_df<-data.frame(end=ends - RElist$ThreeCutSite[i] - 1,chr="simseq")
    revtemp_df$start<-c(revtemp_df$end[1:nrow(revtemp_df)] - 1)
    revtemp_df<-revtemp_df[c("chr","start","end")]
    revtemp_df$end<-replace(revtemp_df$end, which(revtemp_df$end>length(simseq)), length(simseq))
    revtemp_df$start<- replace(revtemp_df$start, which(revtemp_df$start<0), 1)
    forandrev<-rbind(temp_df, revtemp_df)
    forandrev<-arrange(forandrev, end)
    forandrev$start<-c(1, forandrev$end[1:nrow(forandrev)-1] + 1)
    forandrev<-forandrev[c("chr","start","end")]
    forandrev$end<-replace(forandrev$end, which(forandrev$end>length(simseq)), length(simseq))
    forandrev$start<- replace(forandrev$start, which(forandrev$start<0), 1)
    revcomp[[RElist$Enzyme[i]]]<-forandrev
  }
  temp_df$start<-c(1, temp_df$end[1:nrow(temp_df)-1] + 1)
  temp_df<-temp_df[c("chr","start","end")]
  temp_df$end<-replace(temp_df$end, which(temp_df$end>length(simseq)), length(simseq))
  temp_df$start<- replace(temp_df$start, which(temp_df$start<0), 1)
  if(RElist$Seq[i] != RElist$revcompSeq[i]){
    REfrags[[RElist$Enzyme[i]]]<-forandrev
  } else {
    REfrags[[RElist$Enzyme[i]]]<-temp_df
  }
  REfrags[[i]][,"width"]<-REfrags[[i]][,3]-REfrags[[i]][,2]
  REfraglibs[[RElist$Enzyme[i]]]<-REfrags[[i]][which(REfrags[[i]][,4]>39&REfrags[[i]][,4]<600),]
  REfraglibs[[i]][,"Enzyme"]<-rep(RElist$Enzyme[i], nrow(REfraglibs[[i]]))
}
REfraglibs
library(plyr)

ggplot(aes(x=width, col = Enzyme), data = AllLibs) + geom_density(show.legend = FALSE)
library(tidyverse)
fragcounts<-sapply(REfrags, nrow)
fraglibcounts<-sapply(REfraglibs, nrow)
library(GenomicRanges)
library(SomaticSignatures)

res <- lapply(split(AllLibs, AllLibs$Enzyme), function(i){
    GRanges(seqnames = i$chr,
            ranges = IRanges(start = i$start,
                             end = i$end,
                             names = i$Enzyme))})

abasi<-REfraglibs$AbaSI
i = abasi
start=i$start
end=i$end
head(i)
abasi<-GRanges( seqnames="Abasi", ranges = IRanges(start = start,
        width = 3))
gr1 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:4, width = 3))
seqinfo(simseq)

library(BSgenome.Hsapiens.UCSC.hg38)
?writeXStringSet
ss<-DNAStringSet(simseq)
writeXStringSet(ss, "simseq.fa")
ss<-FaFile("simseq.fa")
kmerFrequency(ss)
simseq
Hsapiens
abasi
?FaFile
?kmerFrequency
start
?granges
starts<-abasi$start
ends<-abasi$end
names<-abasi$Enzyme
??seqinfo()
GRanges(seqnames = Rle(abasi$chr),
        ranges = IRanges(start = starts,
                         end = ends,
                         names = names))

res
simseq
kmerFrequency(simseq,ranges = res)
