library(stringr)
library(devtools)
library(tidyverse)
library(forcats)
library(devtools)
library(SomaticSignatures)
library(MutationalPatterns)
library(reshape)
library(cowplot)
library(ggpubr)
mygenome <- FaFile("/Volumes/scratch/cancergeneticslab/ConorM_scratch_Nov21/ref_genomes/mouse2/mouse.fa")
SSM3vars<-read.delim("/Volumes/scratch/cancergeneticslab/SSM3/AllSNPsnoheadernoGERM.txt", sep = " ", header = FALSE)



colnames(SSM3vars)<-c("CHROM", "POS", "REF", "ALT", "sample", "FORMAT")
SSM3vars$AF<-as.numeric(str_extract(SSM3vars$FORMAT, "0\\.\\d{0,}"))
hist(SSM3vars$AF)
SSM3vars<-SSM3vars[SSM3vars$AF<0.7,]
SSM3vars$sample[SSM3vars$AF<0.7 & SSM3vars$AF>0.25]<-"SSM3_clonal"
SSM3vars$sample[SSM3vars$AF<=0.25]<-"SSM3_subclonal"
table(SSM3vars$sample)
SSM3vars$start<-SSM3vars$POS
SSM3vars$end<-SSM3vars$POS
SSM3vars$group<-"SSM3"
SSM3vars$experiment_code<-"WGS"
SSM3vars$subexp<-"WGS"
SSM3varsG<-GRanges(SSM3vars)
SSM3vars_vr <- VRanges(
  seqnames = seqnames(SSM3varsG),
  ranges = ranges(SSM3varsG),
  ref = SSM3varsG$REF,
  alt = SSM3varsG$ALT,
  sampleNames = SSM3varsG$sample,
  AF = SSM3varsG$AF,
  group =SSM3varsG$group)



d<-read.delim("totSNPS.txt", sep = " ", header = FALSE)


d$V8<-as.numeric(str_extract(d$V7, "0\\.\\d{1,2}"))

colnames(d)<-c("CHROM", "POS", "REF", "ALT", "sample", "DP", "FORMAT", "AF")
d$clonality<-paste(str_extract(d$sample, "\\w{0,}clonal"), "mutations")

d$experiment_code<-str_extract(d$sample, "(?<=_)[:alpha:]{0,}(?=_)")
d$experiment_code[grep("SSM3_subclonal",d$sample)]<-"SSM3_WGS"
d$experiment_code[grep("SSM3_clonal",d$sample)]<-"SSM3_WGS"

table(d$experiment_code)
table(d$sample)
table(d$sample)
d$Sample<-str_extract(d$sample, "(?<=_)[:alpha:]{0,}_\\w{0,}")





d$subexp<-str_extract(d$Sample, "(?<=_)\\w{1,1}")

d$group[d$experiment_code == "irrad" & d$subexp =="0"]<-"SSM3 subclones"
d$group[d$experiment_code == "irrad" & d$subexp =="1"]<-"10cGy irradiated"
d$group[d$experiment_code == "irrad" & d$subexp =="5"]<-"5cGy irradiated"
d$group[d$experiment_code == "PIG" & d$subexp =="E"]<-"GFP transfected"
d$group[d$experiment_code == "PIG" & d$subexp =="A"]<-"A3B+GFP transfected"
d$group[d$experiment_code == "SA" & d$subexp =="A"]<-"A3A+UGI transfected"
d$group[d$experiment_code == "SA" & d$subexp =="B"]<-"A3B+UGI transfected"
d$group[d$experiment_code == "SA" & d$subexp =="D"]<-"GFP+mCherry transfected"
d$group[d$experiment_code == "SA" & d$subexp =="O"]<-"0D5 resequenced"
d$group<-as.factor(d$group)
d$group<-fct_relevel(d$group,levels(d$group)[c(4,2,5,1,3)])
d$sampleplot<-str_extract(d$sample, "(?<=_\\w{0,6}_)\\w{0,}")

d$start<-d$POS
d$end<-d$POS
d<-subset(d, experiment_code=="SA")

d<-subset(d, !Sample %in% c("SA_DT1A4_DT1B5", "SA_DT1C2_DT1D12", "SA_B2F11"))
samps<-unique(d$sampleplot)[c(12,8:11,1:7)]
d$sampleplot<-factor(d$sampleplot, levels=samps)
table(d$Sample)
histplots<-ggplot(d, aes(AF))+geom_density()+facet_wrap(sampleplot~.)+
  theme(plot.margin = unit(c(5,5,1,1.8), "lines"), axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5,hjust=1, angle=90),
        axis.title = element_text(size=25), strip.text = element_text(size = 15),
        legend.text = element_text(size=20),legend.title = element_text(size=22))+scale_fill_manual(values = c("blue","red"))+
  xlab("Allele frequency")
ggsave(plot=histplots, file="SA_strat/AFdist.png",width=10,height=10)
histplots<-ggplot(d, aes(AF))+geom_density()+facet_wrap(sampleplot~.)+
  theme(plot.margin = unit(c(5,5,1,1.8), "lines"), axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5,hjust=1, angle=90),
        axis.title = element_text(size=25), strip.text = element_text(size = 15),
        legend.text = element_text(size=20),legend.title = element_text(size=22))+scale_fill_manual(values = c("blue","red"))+
  xlab("Allele frequency")+scale_y_continuous(trans="log10")+ylab("Density (log scale)")
ggsave(plot=histplots, file="SA_strat/AFdist.png",width=10,height=10)
nmuts<-table(d$sample)

nmuts

depth<-read.delim("depthtab.txt")
depth<-depth[,1:4]
depth<-depth[!str_detect(depth$sample, "DT.{0,}\\d_DT"),]
depth
depth<-depth[match(str_extract(names(nmuts),"(?<=_)\\w{0,}"), depth$sample),]
depth$nmuts<-as.numeric(nmuts)
depth$sample<-names(nmuts)
depth$Sample<-str_extract(names(nmuts),"(?<=_)\\w{0,}")
SSM3depth<-read.delim("SSM3depthtab.txt", header =FALSE, sep =" ")
SSM3depth<-data.frame(sample=unique(SSM3vars$sample)[c(2,1)],Sample=unique(SSM3vars$sample)[c(2,1)],depth=rep(30,2),stdev=rep(20,2), greater20=rep(SSM3depth$V2,2), nmuts=as.vector(table(SSM3vars$sample)))
depth<-bind_rows(SSM3depth,depth)
depth$experiment_code<-str_extract(depth$Sample, "^\\w{1,2}")
depth$subexp<-str_extract(depth$Sample, "(?<=_)\\w{1,1}")
depth$group[depth$experiment_code == "ir" & depth$subexp =="0"]<-"SSM3 subclones"
depth$group[depth$experiment_code == "ir" & depth$subexp =="1"]<-"10cGy irradiated"
depth$group[depth$experiment_code == "ir" & depth$subexp =="5"]<-"5cGy irradiated"
depth$group[depth$experiment_code == "PI" & depth$subexp =="E"]<-"GFP transfected"
depth$group[depth$experiment_code == "PI" & depth$subexp =="A"]<-"A3B+GFP transfected"
depth$group[depth$experiment_code == "SA" & depth$subexp =="A"]<-"A3A+UGI transfected"
depth$group[depth$experiment_code == "SA" & depth$subexp =="B"]<-"A3B+UGI transfected"
depth$group[depth$experiment_code == "SA" & depth$subexp =="D"]<-"GFP+mCherry transfected"
depth$group[depth$experiment_code == "SA" & depth$subexp =="O"]<-"0D5 resequenced"
depth$group[depth$experiment_code == "SS"]<-"WGS SSM3"
depth$mutrate<-(depth$nmuts*10^6)/depth$greater20
depth$clonality<-paste(str_extract(depth$sample, "\\w{0,}clonal"), "mutations")
depth$clonality[c(1,2)]<-c("clonal mutations","subclonal mutations")
depth<-subset(depth, log(greater20)>14)
depth$sampleplot<-str_extract(depth$Sample, "(?<=_)\\w{0,}")

depth$sampleplot[1]<-"SSM3"
depth$sampleplot[2]<-"SSM3"
depth$group[1:2]<-"WGS SSM3"
depth$sampleplot<-as.factor(depth$sampleplot)
depth$sampleplot<-fct_relevel(depth$sampleplot, levels(depth$sampleplot)[c(13,12,8:11,1:7)])
depth$group<-factor(depth$group)
depth$group<-fct_relevel(depth$group,levels(depth$group)[c(5,4,2,1,3)])

depth

depthplot<-ggplot(depth, aes(sampleplot,))


mutrateplot<-ggplot(depth, aes(sampleplot,mutrate,fill=group))+
  geom_bar(stat="identity")+
  theme(axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5,hjust=1, angle=90),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),legend.title = element_text(size=22), strip.text = element_text(size=15))+
  guides(fill=guide_legend(title="Group"))+ylab("Mutation rate (mutations/Mb,log2 scale)")+xlab("Sample")+
  coord_cartesian(clip="off")+facet_grid(rows=vars(clonality))+scale_y_continuous(trans='log2', breaks = c(0,10,100,200))
ggsave(plot=mutrateplot, file="SA_strat/mutrateplot.png",width=10,height=10)
mutrateplot
sub_mut<-subset(depth, clonality =="subclonal mutations" & subexp %in% c("A","B","D"))
a<-aov(mutrate ~ subexp,sub_mut)
summary(a)
clonal_mut<-subset(depth, clonality =="clonal mutations" & subexp %in% c("A","B","D"))
a<-aov(mutrate ~ subexp,clonal_mut)
summary(a)
TukeyHSD(a)
unique(d$experiment_code)
allvars<-bind_rows(d,SSM3vars)
allvars$clonality[grep("SSM3_clonal", allvars$sample)]<-"clonal mutations"
allvars$clonality[grep("SSM3_subclonal", allvars$sample)]<-"subclonal mutations"
allvars$sampleplot<-str_extract(allvars$sample, "(?<=_)\\w{0,}")
allvars$sampleplot[allvars$sampleplot=="clonal mutations"]<-"SSM3_clonal"
allvars$sampleplot[allvars$sampleplot=="subclonal mutations"]<-"SSM3_subclonal"
allvars$sampleplot[allvars$sampleplot=="OD5"]<-"0D5 resequenced"
dG<-GRanges(allvars)
d_vr <- VRanges(
  seqnames = seqnames(dG),
  ranges = ranges(dG),
  ref = dG$REF,
  alt = dG$ALT,
  sampleNames = dG$sample,
  AF = dG$AF,
  clonality=dG$clonality,
  group =dG$group,
  sampleplot=dG$sampleplot)
mut_context<-mutationContext(d_vr, mygenome)
d_mm<-motifMatrix(mut_context,normalize = TRUE, group="sampleNames")
mut_context
mat<-cos_sim_matrix(d_mm,d_mm)
library(dplyr)
library(cowplot)
mat<-melt(mat)
mat$X1<-factor(mat$X1)
mat$X2<-factor(mat$X2)
levels(mat$X2)
levels(mat$X1)
mat$X1<-fct_relevel(mat$X1, levels(mat$X1)[c(13,14,12,8:11,1:7,26,22:25,15:21)])
mat$X2<-fct_relevel(mat$X2, levels(mat$X2)[c(13,14,12,8:11,1:7,26,22:25,15:21)])
cosmat0<-ggplot(mat, aes(X1,X2, fill=value))+geom_tile()+geom_text(aes(label=signif(value,2)),size=5)+
  theme(axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5,hjust=1, angle=90),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),legend.title = element_text(size=22))+
  scale_fill_gradient(low="white",high="red",guide_legend(title="Cosine \nsimilarity"))+
  xlab("Sample")+ylab("Sample")
ggsave(plot=cosmat0, file="SA_strat/cos_heatmap0.png",width=20,height=20)
##evaluating SigProf results
solstats<-NULL
for(i in 1:10){
  path=paste0("sigprofresultsSA_hets_sep/SBS96/All_Solutions/SBS96_",i,"_Signatures/Solution_Stats/SBS96_S",i,"_Samples_stats.txt")
  ss<-read.delim(path)
  ss$n_sigs<-i
  solstats[[i]]<-ss
}
solstats<-bind_rows(solstats)

solstats<-solstats[solstats$Sample.Names %in% depth$sample,]
head(solstats)
nsigscosim<-ggplot(solstats, aes(n_sigs, Cosine.Similarity, col=Sample.Names))+geom_line()+
  scale_x_continuous(breaks=1:10)+
  theme(axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),legend.title = element_text(size=22))+
  scale_fill_gradient(low="white",high="red",guide_legend(title="Cosine \nsimilarity"))+
  xlab("Number of Signatures")+ylab("Cosine similarity of \nreconstructed ptofile \nto original profile")
ggsave(plot=nsigscosim, file="SA_strat/nsigscosim.png",width=20,height=12)
nsigscosim+theme(legend.position = "none")
solstats3<-subset(solstats, n_sigs==3)
cossim3<-ggplot(solstats3, aes(Total.Mutations, Cosine.Similarity, col=Sample.Names))+
  geom_point(size=5)+scale_x_continuous(trans='log2')+
  theme(axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),legend.title = element_text(size=22))+
  scale_fill_gradient(low="white",high="red",guide_legend(title="Cosine \nsimilarity"))+
  geom_text(data=subset(solstats3, Sample.Names=="clonal_SA_A2H1"), aes(label=Sample.Names),size=10,vjust=-1)+
  xlab("Total number of mutations in sample")+ylab("Cosine similarity of reconstructed ptofile \nto original profile (at 3 extracted signatures)")
ggsave(plot=cossim3, file="SA_strat/cossim3.png",width=20,height=12)
sigprofres<-read.delim("sigprofresultsSA_hets_sep/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt")
cossim3+theme(legend.position = "none")
#sigprofres<-sigprofres[grep("subclonal_",sigprofres$Samples,invert = TRUE),]

sigprofres<-sigprofres[sigprofres$Samples %in% depth$sample,]

sigprofresnorm<-data.frame(Samples=sigprofres$Samples,sigprofres[,2:4]/rowSums(sigprofres[,2:4]))
spn<-melt(sigprofresnorm)
spn$clonality<-paste(str_extract(spn$Samples, "\\w{0,}clonal"), "mutations")
spn$clonality[spn$clonality=="SSM3_clonal mutations"]<-"clonal mutations"
spn$clonality[spn$clonality=="SSM3_subclonal mutations"]<-"subclonal mutations"
spn$sampleplot<-str_extract(spn$Samples, "(?<=_\\w{0,3}_)\\w{0,}$")
spn$sampleplot[grep("SSM3",spn$Samples)]<-"SSM3"
spn$sampleplot<-factor(spn$sampleplot,levels(depth$sampleplot))
spn

sigprofplot<-ggplot(spn,aes(sampleplot, fill=variable,value))+geom_bar(stat="identity")+
  theme(axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5,hjust=1, angle=90),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),legend.title = element_text(size=22),
        strip.text = element_text(size=20))+
  ylab("Proportion of mutations \nattributed to signature")+xlab("Sample")+scale_fill_brewer(palette = "Set2",guide_legend(title="SigProfiler \nextracted \nsignature"))+facet_grid(rows=vars(clonality))
ggsave(plot=sigprofplot, file="SA_strat/sigprofplot.png",width=10,height=10)
sigprofplot


Mc<-data.frame(mut_context)
Mc$Context<-"Other"
Mc$Context[Mc$context =="T.A" & Mc$alteration=="CT"]<-"T(C>T)W"
Mc$Context[Mc$context =="T.T" & Mc$alteration=="CT"]<-"T(C>T)W"

Mc_clonal<-subset(Mc, clonality=="clonal mutations")
Mc_clonal$sampleNames<-as.character(Mc_clonal$sampleNames)
Mc_clonal$sampleplot<-str_extract(Mc_clonal$sampleNames,"(?<=_\\w{0,3}_)\\w{0,}")
Mc_clonal$sampleplot[is.na(Mc_clonal$sampleplot)]<-"SSM3"

Mc_clonal$sampleplot<-factor(Mc_clonal$sampleplot, levels(depth$sampleplot))
levels(Mc_clonal$sampleplot)

dfprops<-data.frame(prop.table(table(Mc_clonal$Context, Mc_clonal$sampleplot),2))
dfprops<-subset(dfprops, Var1 == "T(C>T)W")
expected<-dfprops$Freq[dfprops$Var2=="OD5"]
dfcounts<-data.frame(table(Mc_clonal$Context, Mc_clonal$sampleplot))
dfcounts
samples<-levels(dfcounts$Var2)[2:length(levels(dfcounts$Var2))]
samples
x=1
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
binomtab$p.adj<-p.adjust(binomtab$p.val, method ="holm")
binomtab$sig<-"ns"
binomtab$sig[binomtab$p.adj<0.05]<-"*"
binomtab$sig[binomtab$p.adj<0.01]<-"**"
binomtab$sig[binomtab$p.adj<0.001]<-"***"
binomtab$sig[1]<-"Ref"
TCWplot_clonal<-ggplot(Mc_clonal, aes(sampleplot))+geom_bar(aes(fill=Context),position = "fill")+
  geom_text(binomtab, mapping =aes(x=sampleNames, y=Inf,vjust=-1.5, label=sig), size =5)+
  xlab("Sample")+ylab("Proportion of mutations \n in context")+
  coord_cartesian(ylim = c(0, 1), clip = 'off')+
  theme(plot.margin = unit(c(3,1,1,1), "lines"), axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5,hjust=1, angle=90),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),legend.title = element_text(size=22))+scale_fill_manual(values = c("blue","red"))
ggsave(plot=TCWplot_clonal, file="SA_strat/TCWplot_clonal.png",width=10,height=5)
TCWplot_clonal

Mc_subclonal<-subset(Mc, clonality=="subclonal mutations")

Mc_subclonal$sampleNames<-as.character(Mc_subclonal$sampleNames)
Mc_subclonal$sampleplot<-str_extract(Mc_subclonal$sampleNames,"(?<=_\\w{0,3}_)\\w{0,}")
Mc_subclonal$sampleplot[is.na(Mc_subclonal$sampleplot)]<-"SSM3"

Mc_subclonal$sampleplot<-factor(Mc_subclonal$sampleplot, levels(depth$sampleplot))
levels(Mc_subclonal$sampleplot)

dfprops<-data.frame(prop.table(table(Mc_subclonal$Context, Mc_subclonal$sampleplot),2))
dfprops<-subset(dfprops, Var1 == "T(C>T)W")
expected<-dfprops$Freq[dfprops$Var2=="OD5"]
dfcounts<-data.frame(table(Mc_subclonal$Context, Mc_subclonal$sampleplot))
dfcounts
samples<-levels(dfcounts$Var2)[2:length(levels(dfcounts$Var2))]
samples
x=1
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
binomtab$p.adj<-p.adjust(binomtab$p.val, method ="holm")
binomtab$sig<-"ns"
binomtab$sig[binomtab$p.adj<0.05]<-"*"
binomtab$sig[binomtab$p.adj<0.01]<-"**"
binomtab$sig[binomtab$p.adj<0.001]<-"***"
binomtab$sig[1]<-"Ref"
TCWplot_subclonal<-ggplot(Mc_subclonal, aes(sampleplot))+geom_bar(aes(fill=Context),position = "fill")+
  geom_text(binomtab, mapping =aes(x=sampleNames, y=Inf,vjust=-1.5, label=sig), size =5)+
  xlab("Sample")+ylab("Proportion of mutations \n in context")+
  coord_cartesian(ylim = c(0, 1), clip = 'off')+
  theme(plot.margin = unit(c(3,1,1,1), "lines"), axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5,hjust=1, angle=90),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),legend.title = element_text(size=22))+scale_fill_manual(values = c("blue","red"))
ggsave(plot=TCWplot_subclonal, file="SA_strat/TCWplot_subclonal.png",width=10,height=5)
TCWplot_subclonal
TCWplot<-ggarrange(TCWplot_clonal+theme(axis.text.x = element_blank(), axis.title.x = element_blank()), TCWplot_subclonal,nrow=2,ncol=1, heights=c(1,1.5),common.legend = TRUE)
TCWplot
ggsave(plot=TCWplot, file="SA_strat/TCWplot.png",width=10,height=10)

