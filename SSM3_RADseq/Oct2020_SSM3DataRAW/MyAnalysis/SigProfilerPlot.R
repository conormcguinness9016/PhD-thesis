library(reticulate)
library(SigProfilerExtractorR)
library(SigProfilerMatrixGeneratorR)
library(SigProfilerPlottingR)
b<-SigProfilerMatrixGeneratorR("BRCA", "mm10", "/Volumes/archive/cancergeneticslab/ConorM/Oct2020_SSM3DataRAW/MyAnalysis/totalVCFhets/", plot=T, exome=F, bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F, cushion=100)
plotSBS(b + "96", "sigprofiler/","BRCA", "96",FALSE)
b<-SigProfilerMatrixGeneratorR("BRCA", "mm10", "/Volumes/archive/cancergeneticslab/ConorM/Oct2020_SSM3DataRAW/MyAnalysis/uniqueVCFhets/", plot=T, exome=F, bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F, cushion=100)
plotSBS(b + "96", "sigprofiler/","BRCA", "96",FALSE)
d<-read.delim("unSNPS.txt", sep = " ", header = FALSE)
library(stringr)
library(devtools)
library(tidyverse)
library(forcats)
library(devtools)
library(SomaticSignatures)
d$V8<-as.numeric(str_extract(d$V7, "0\\.\\d{1,2}"))

colnames(d)<-c("CHROM", "POS", "REF", "ALT", "sample", "DP", "FORMAT", "AF")
dir.create("Rplotsgood")
png("Rplotsgood/varhist.png")
hist(d$AF)
dev.off()
table()

d$group_code<-str_extract(d$sample, "^\\w{1,1}")
d$group[d$group_code == "E"]<-"Empty Vector+GFP"
d$group[d$group_code == "A"]<-"APOBEC3B + GFP"


d<-subset(d, AF>0.30)
nmuts<-table(d$sample)

depth<-read.delim("depthtab.txt")
depth<-depth[!str_detect(depth$sample, "GBSN"),]

depth
depth<-depth[match(names(nmuts), depth$sample),]
depth$nmuts<-nmuts
depth<-na.omit(depth)
depth$mutrate<-(depth$nmuts*10^6)/depth$greater20
depth$group_code<-str_extract(depth$sample, "^\\w{1,1}")
depth$group[depth$group_code == "E"]<-"Empty Vector+GFP"
depth$group[depth$group_code == "A"]<-"APOBEC3B + GFP"
depth
depth$group<-fct_relevel(depth$group,c("Empty Vector+GFP","APOBEC3B + GFP"))
dir.create("Rplotsgood")
png("Rplotsgood/mutationrate.png", width =1500, height=1500)
ggplot(depth, aes(group, mutrate)) + geom_boxplot()+
  geom_point(aes(col = sample,size= 40)) +
  xlab("Construct Transfected")+ylab("Mutation Rate \n (mutations/Mb sequenced)")+
  guides(color = guide_legend(override.aes = list(size = 10) ) )+
  theme(axis.title = element_text(size = 50), axis.text= element_text(size = 35),
        axis.text.x=element_text(angle = 45, hjust=1),
        legend.text=element_text(size =45), legend.title = element_text(size=50))
dev.off()
summary(aov(mutrate ~ group, data = depth))
d$start<-d$POS
d$end<-d$POS
mygenome <- FaFile("/Volumes/archive/cancergeneticslab/ConorM_scratch_Nov21/ref_genomes/mouse2/mouse.fa")
depth<-depth[!str_detect(depth$sample, "_D"),]
dA<-subset(d, sample == "A2H1")
dG<-GRanges(d)
d_vr <- VRanges(
  seqnames = seqnames(dG),
  ranges = ranges(dG),
  ref = dG$REF,
  alt = dG$ALT,
  sampleNames = dG$sample,
  AF = dG$AF,
  group =dG$group)
mut_context<-mutationContext(d_vr, mygenome)
Mc<-data.frame(mut_context)
Mc$Context<-"Other"
Mc$Context[str_detect(Mc$context, "^T") & Mc$alteration=="CT"]<-"T(C>T)N"
chisq.test(Mc$Context, Mc$sampleNames)
dfprops<-data.frame(prop.table(table(Mc$Context, Mc$sampleNames),2))
dfprops<-subset(dfprops, Var1 == "T(C>T)N")
expected<-median(dfprops$Freq)
dfcounts<-data.frame(table(Mc$Context, Mc$sampleNames))
samples<-levels(dfcounts$Var2)
t<-lapply(seq(samples),function(x){
  s<-samples[x]
  df<-subset(dfcounts, Var2 ==s)
  muts<-sum(df$Freq)
  A3muts<-df$Freq[2]
  binom<-binom.test(A3muts, muts, expected, alternative = "greater")
  res<-data.frame(sampleNames=s,
                  p.val=binom$p.value,
                  lwr_95=binom$conf.int[1],
                  upr_95=binom$conf.int[2], 
                  estimate=as.numeric(binom$estimate),
                  mutations=as.numeric(binom$parameter))
})
binomtab<-bind_rows(t)
binomtab$p.adj<-p.adjust(binomtab$p.val)
binomtab$y<-1.01
ggsave("Rplotsgood/mutproportions.png")
ggplot(Mc, aes(sampleNames))+geom_bar(aes(fill=Context),position = "fill")+
  geom_text(binomtab, mapping =aes(x=sampleNames, y=Inf,vjust=-1, label=signif(p.adj,3)), size =3)+
  annotate("text",x=-Inf,y=Inf,hjust=1,vjust =-1,label ="Adjusted p-value", size =3 )+
  geom_text(binomtab, mapping =aes(x=sampleNames, y=Inf,vjust=-3, label=mutations), size =3)+
  annotate("text",x=-Inf,y=Inf,hjust=1,vjust =-3,label ="Number of mutations", size =3 )+
  xlab("Sample")+ylab("Proportion of mutations \n in context")+

  coord_cartesian(ylim = c(0, 1), clip = 'off')+
  theme(plot.margin = unit(c(3,1,1,1.8), "lines"),axis.title = element_text(size = 15), axis.text= element_text(size = 10),
        axis.text.x=element_text(angle = 45, hjust=1),
        legend.text=element_text(size =10), legend.title = element_text(size=10.0))

ggplot(Mc, aes(sampleNames, AF, col=Context))+geom_boxplot()
summary(aov(AF~Context*sampleNames, Mc))
library(ggstatsplot)
png("Rplotsgood/fancybarplot.png", width =1500, height =1500)
ggbarstats(Mc, 
           Context,
           sampleNames)
?ggbarstats
dev.off()