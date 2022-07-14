library(reticulate)
library(SigProfilerExtractorR)
library(SigProfilerMatrixGeneratorR)
library(SigProfilerPlottingR)
use_python("/usr/bin/python3.6")

b<-SigProfilerMatrixGeneratorR("BRCA", "mm10", "/Volumes/archive/cancergeneticslab/ConorM/SigsCombinedAnalysis/germremoved_SA_SigProfileranalysis_Hets_sep", plot=T, exome=F, bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F, cushion=100)
plotSBS(b + "96", "sigprofiler/","BRCA", "96",FALSE)




library(stringr)
library(devtools)
library(tidyverse)
library(forcats)
library(devtools)
library(SomaticSignatures)
library(MutationalPatterns)
library(reshape)
library(cowplot)
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
nmuts<-table(d$sample[d$AF>0.2])
nmuts
d$experiment_code<-str_extract(d$sample, "^\\w{1,1}")
d$subexp<-str_extract(d$sample, "(?<=_)\\w{1,1}")

d$group[d$experiment_code == "i" & d$subexp =="0"]<-"SSM3 subclones"
d$group[d$experiment_code == "i" & d$subexp =="1"]<-"10cGy irradiated"
d$group[d$experiment_code == "i" & d$subexp =="5"]<-"5cGy irradiated"
d$group[d$experiment_code == "P" & d$subexp =="E"]<-"GFP transfected"
d$group[d$experiment_code == "P" & d$subexp =="A"]<-"A3B+GFP transfected"
d$group[d$experiment_code == "S" & d$subexp =="A"]<-"A3A+UGI transfected"
d$group[d$experiment_code == "S" & d$subexp =="B"]<-"A3B+UGI transfected"
d$group[d$experiment_code == "S" & d$subexp =="D"]<-"GFP+mCherry transfected"
d$group[d$experiment_code == "S" & d$subexp =="O"]<-"0D5 resequenced"
d$group<-as.factor(d$group)
d$group<-fct_relevel(d$group,levels(d$group)[c(4,2,5,1,3)])
d$sampleplot<-str_extract(d$sample, "(?<=_)\\w{0,}")
d$sampleplot[d$sampleplot=="OD5"]<-"0D5_resequenced"
d$start<-d$POS
d$end<-d$POS
depth<-read.delim("depthtab.txt")
depth<-depth[,1:4]
depth<-depth[!str_detect(depth$sample, "DT.{0,}\\d_DT"),]
depth<-depth[match(names(nmuts), depth$sample),]
depth$nmuts<-as.numeric(nmuts)
depth<-na.omit(depth)
SSM3depth<-read.delim("SSM3depthtab.txt", header =FALSE, sep =" ")
SSM3depth<-data.frame(sample=unique(SSM3vars$sample)[c(2,1)],depth=rep(30,2),stdev=rep(20,2), greater20=rep(SSM3depth$V2,2), nmuts=as.vector(table(SSM3vars$sample)))
depth<-bind_rows(SSM3depth,depth)
depth$mutrate<-(depth$nmuts*10^6)/depth$greater20
depth
depth$experiment_code<-str_extract(depth$sample, "^\\w{1,2}")
depth$subexp<-str_extract(depth$sample, "(?<=_)\\w{1,1}")
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
depth$group<-as.factor(depth$group)
depth<-subset(depth, log(greater20)>14)
depth$sampleplot<-str_extract(depth$sample, "(?<=_)\\w{0,}")
depth$sampleplot[1]<-"SSM3_clonal"
depth$sampleplot[2]<-"SSM3_subclonal"
depth$sampleplot[depth$group=="0D5 resequenced"]<-"0D5 resequenced"
depth$sampleplot<-as.factor(depth$sampleplot)
depth$sampleplot<-fct_relevel(depth$sampleplot, levels(depth$sampleplot)[c(33,34,1:3,5:6,4,24:27,13:23,28:32,7:12)])

depth$group<-fct_relevel(depth$group,levels(depth$group)[c(10,9,1,8,4,6,3,2,7,5)])
depthSS<-depth[depth$experiment_code=="SS",]
depth0<-depth[depth$experiment_code=="ir" & depth$group == "SSM3 subclones",]
depthAPO<-depth[depth$experiment_code=="SA",]
depthreseq<-depth[depth$group == "0D5 resequenced",]
depth_prim<-bind_rows(depthSS,depth0,depthreseq)
levels(depth_prim$group)
levels(depth_prim$sampleplot)
mutrateplot0<-ggplot(depth_prim, aes(sampleplot,mutrate,fill=group))+
  geom_bar(stat="identity")+
  theme(axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5,hjust=1, angle=90),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),legend.title = element_text(size=22))+
  guides(fill=guide_legend(title="Group"))+ylab("Mutation rate (mutations/Mb)")+
  coord_cartesian(clip="off")
ggsave(plot=mutrateplot0, file="Rplotsgood/mutrateplot0.png",width=10,height=10)
unirrad<-subset(d, d$group=="SSM3 subclones" | d$group == "0D5 resequenced")
gghistogram(unirrad[unirrad$AF>0.18,],x="AF",fill="sampleplot",merge=FALSE,bins=10,facet.by = "sampleplot")

# d<-d[d$sample %in% depth_prim$sample,]
# unique(d$sample)

unique(d$experiment_code)
allvars<-bind_rows(subset(d,group=="0D5 resequenced" & AF>0.2),subset(d,experiment_code=="i" & group == "SSM3 subclones" ),SSM3vars)
allvars$sampleplot<-str_extract(allvars$sample, "(?<=_)\\w{0,}")
allvars$sampleplot[allvars$sampleplot=="clonal"]<-"SSM3_clonal"
allvars$sampleplot[allvars$sampleplot=="subclonal"]<-"SSM3_subclonal"
allvars$sampleplot[allvars$sampleplot=="OD5"]<-"0D5 resequenced"
dG<-GRanges(allvars)
d_vr <- VRanges(
  seqnames = seqnames(dG),
  ranges = ranges(dG),
  ref = dG$REF,
  alt = dG$ALT,
  sampleNames = dG$sample,
  AF = dG$AF,
  group =dG$group,
  sampleplot=dG$sampleplot)
mut_context<-mutationContext(d_vr, mygenome)
d_mm<-motifMatrix(mut_context,normalize = TRUE, group="sampleplot")
mat<-cos_sim_matrix(d_mm,d_mm)
library(dplyr)
library(cowplot)
mat<-melt(mat)
mat$X1<-factor(mat$X1)
mat$X2<-factor(mat$X2)
levels(mat$X1)
mat$X1<-fct_relevel(mat$X1, levels(mat$X1)[c(7,8,1:3,5,6,4)])
mat$X2<-fct_relevel(mat$X2, levels(mat$X2)[c(7,8,1:3,5,6,4)])
cosmat0<-ggplot(mat, aes(X1,X2, fill=value))+geom_tile()+geom_text(aes(label=signif(value,2)),size=5)+
  theme(axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5,hjust=1, angle=90),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),legend.title = element_text(size=22))+
  scale_fill_gradient(low="white",high="red",guide_legend(title="Cosine \nsimilarity"))+
  xlab("Sample")+ylab("Sample")
ggsave(plot=cosmat0, file="Rplotsgood/cos_heatmap0.png",width=10,height=10)

sigprofres<-read.delim("sigprofresultsSA_hets_sep/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt")
sigprofres$X96AB<-sigprofres$X96A+sigprofres$X96B
sigprofres<-sigprofres[grep("subclonal_",sigprofres$Samples,invert = TRUE),]
sigprofres<-sigprofres[,c(1,4,5,6)]

sigprofres$Samples<-str_extract(sigprofres$Samples, "(?<=_)\\w{0,}")
sigprofres$Samples[sigprofres$Samples=="clonal"]<-"SSM3_clonal"
sigprofres$Samples[sigprofres$Samples=="subclonal"]<-"SSM3_subclonal"
sigprofres$Samples[sigprofres$Samples=="SA_OD5"]<-"0D5 resequenced"
sigprofresnorm<-data.frame(Samples=sigprofres$Samples,sigprofres[,2:4]/rowSums(sigprofres[,2:4]))
spn<-melt(sigprofresnorm)
spn$sampleplot<-spn$Samples
spn$sampleplot[grep("SA",spn$Samples)]<-str_extract(spn$Samples[grep("SA",spn$Samples)], "(?<=_)\\w{0,}")

spn0<-spn[spn$sampleplot %in% mat$X1,]

spn0$sampleplot<-factor(spn0$sampleplot)
levels(spn0$sampleplot)
spn0$sampleplot<-fct_relevel(spn0$sampleplot, levels(spn0$sampleplot)[c(7,8,1:3,5:6,4)])
sigprofplot0<-ggplot(spn0,aes(sampleplot, fill=variable,value))+geom_bar(stat="identity")+
  theme(axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5,hjust=1, angle=90),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),legend.title = element_text(size=22))+
  ylab("Proportion of mutations \nattributed to signature")+xlab("Sample")+scale_fill_brewer(palette = "Set2",guide_legend(title="SigProfiler \nextracted \nsignature"))
ggsave(plot=sigprofplot0, file="Rplotsgood/sigprofplot0.png",width=10,height=10)


#scale_fill_manual(values=c("black","dark green"),guide_legend(title="Cosine \nsimilarity"))
grid<-plot_grid(mutrateplot0+theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                cosmat0+theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                             axis.title.y = element_blank()),
                sigprofplot0,
                ncol=1,nrow=4, align="v",axis="rl")
grid
save_plot( grid,file="Rplotsgood/Fig7.pdf",base_width=20,base_height=24)

depth_APO<-bind_rows(depthSS,depthAPO)
levels(depth_APO$group)
levels(depth_APO$sampleplot)
mutrateplotAPO<-ggplot(depth_APO, aes(sampleplot,mutrate,fill=group))+
  geom_bar(stat="identity")+
  theme(axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5,hjust=1, angle=90),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),legend.title = element_text(size=22))+
  guides(fill=guide_legend(title="Group"))+ylab("Mutation rate (mutations/Mb)")+
  coord_cartesian(clip="off")
ggsave(plot=mutrateplotAPO, file="Rplotsgood/mutrateplotAPO.png",width=10,height=10)


allvars<-bind_rows(subset(d,group=="0D5 resequenced" & AF>0.2),d[d$sampleplot %in% depth_APO$sampleplot & d$AF>0.2,],SSM3vars)
allvars$sampleplot<-str_extract(allvars$sample, "(?<=_)\\w{0,}")
allvars$sampleplot[allvars$sampleplot=="clonal"]<-"SSM3_clonal"
allvars$sampleplot[allvars$sampleplot=="subclonal"]<-"SSM3_subclonal"
allvars$sampleplot[allvars$sampleplot=="OD5"]<-"0D5 resequenced"
dG<-GRanges(allvars)
d_vr <- VRanges(
  seqnames = seqnames(dG),
  ranges = ranges(dG),
  ref = dG$REF,
  alt = dG$ALT,
  sampleNames = dG$sample,
  AF = dG$AF,
  group =dG$group,
  sampleplot=dG$sampleplot)
mut_context<-mutationContext(d_vr, mygenome)
d_mm<-motifMatrix(mut_context,normalize = TRUE, group="sampleplot")
mat<-cos_sim_matrix(d_mm,d_mm)

mat<-melt(mat)
mat$X1<-factor(mat$X1)
mat$X2<-factor(mat$X2)
levels(mat$X1)
levels(mat$X2)
mat$X1<-fct_relevel(mat$X1, levels(mat$X1)[c(13,14,1,9:12,2:8)])
mat$X2<-fct_relevel(mat$X2, levels(mat$X2)[c(13,14,1,9:12,2:8)])
cosmatAPO<-ggplot(mat, aes(X1,X2, fill=value))+geom_tile()+geom_text(aes(label=signif(value,2)),size=5)+
  theme(axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5,hjust=1, angle=90),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),legend.title = element_text(size=22))+
  scale_fill_gradient(low="white",high="red",guide_legend(title="Cosine \nsimilarity"))+
  xlab("Sample")+ylab("Sample")
ggsave(plot=cosmatAPO, file="Rplotsgood/cos_heatmapAPO.png",width=10,height=10)


spnAPO<-spn[spn$sampleplot %in% mat$X1,]
spnAPO
spnAPO$sampleplot<-factor(spnAPO$sampleplot, levels =levels(mat$X1))
levels(spnAPO$sampleplot)
sigprofplotAPO<-ggplot(spnAPO,aes(sampleplot, fill=variable,value))+geom_bar(stat="identity")+
  theme(axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5,hjust=1, angle=90),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),legend.title = element_text(size=22))+
  ylab("Proportion of mutations \nattributed to signature")+xlab("Sample")+scale_fill_brewer(palette = "Set2",guide_legend(title="SigProfiler \nextracted \nsignature"))
ggsave(plot=sigprofplotAPO, file="Rplotsgood/sigprofplotAPO.png",width=10,height=10)
Mc<-data.frame(mut_context)
Mc$Context<-"Other"
Mc$Context[Mc$context =="T.A" & Mc$alteration=="CT"]<-"T(C>T)W"
Mc$Context[Mc$context =="T.T" & Mc$alteration=="CT"]<-"T(C>T)W"
Mc$sampleNames<-as.character(Mc$sampleNames)



Mc$sampleplot<-factor(Mc$sampleplot)
levels(Mc$sampleplot)
Mc$sampleplot<-fct_relevel(Mc$sampleplot, levels(Mc$sampleplot)[c(13:14,1,9:12,2:8)])


dfprops<-data.frame(prop.table(table(Mc$Context, Mc$sampleplot),2))
dfprops<-subset(dfprops, Var1 == "T(C>T)W")
expected<-dfprops$Freq[dfprops$Var2=="0D5 resequenced"]
dfcounts<-data.frame(table(Mc$Context, Mc$sampleplot))

samples<-levels(dfcounts$Var2)[3:length(levels(dfcounts$Var2))]
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
TCWplot<-ggplot(Mc, aes(sampleplot))+geom_bar(aes(fill=Context),position = "fill")+
  geom_text(binomtab, mapping =aes(x=sampleNames, y=Inf,vjust=-1.5, label=sig), size =10)+
  # annotate("text",x=-Inf,y=Inf,hjust=1,vjust =-0.15,label ="Significance \n (Holm adjusted)", size =5 )+
  # geom_text(binomtab, mapping =aes(x=sampleNames, y=Inf,vjust=-6.25, label=mutations), size =5)+
  # annotate("text",x=-Inf,y=Inf,hjust=1,vjust =-2,label ="Number of unique \n mutations", size =3 )+
  xlab("Sample")+ylab("Proportion of mutations \n in context")+
  coord_cartesian(ylim = c(0, 1), clip = 'off')+
  theme(plot.margin = unit(c(5,5,1,1.8), "lines"), axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5,hjust=1, angle=90),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),legend.title = element_text(size=22))+scale_fill_manual(values = c("blue","red"))

TCWplot

grid<-plot_grid(mutrateplotAPO+theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                cosmatAPO+theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                              axis.title.y = element_blank()),
                sigprofplotAPO+theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                TCWplot,
                ncol=1,nrow=4, align="v",axis="rl")
grid
save_plot( grid,file="Rplotsgood/Fig8.pdf",base_width=20,base_height=24)




sigprofprob<-read.delim("sigprofresultsSA_hets_sep/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/De_Novo_Mutation_Probabilities_refit.txt")
probs<-subset(sigprofprob, Sample.Names=="irrad_0D5" & X96D >0.9)
sigprofprob
probs
#depth<-depth[depth$experiment_code != "i",]


depthSA<-bind_rows(depthSA,depth0)
depthSA$sampleplot<-as.factor(depthSA$sampleplot)
levels(depthSA$sampleplot)
depthSA$sampleplot<-fct_relevel(depthSA$sampleplot, levels(depthSA$sampleplot)[c(13,14,12,8:11,1:3,4:7)])
dir.create("Rplotsgood")
mutrateplot<-ggplot(depthSA, aes(sampleplot,mutrate,group=group,fill=group))+
  geom_bar(stat="identity")+
  theme(axis.text=element_text(size =20), axis.text.x = element_text(vjust=0.5,hjust=1, angle=90),
        axis.title = element_text(size=25),
        legend.text = element_text(size=20),legend.title = element_text(size=22))+
  guides(fill=guide_legend(title="Group"))+ylab("Mutation rate (mutations/Mb)")+
  coord_cartesian(clip="off")

ggsave(plot=mutrateplot, file="Rplotsgood/mutrateplot.png",width=10,height=10)
probs
OD5mut_context<-data.frame(mut_context[sampleNames(mut_context)=="SA_OD5",])

Fiveprimecontext<-str_extract(OD5mut_context$context,"\\w(?=\\.)")
threeprimecontext<-str_extract(OD5mut_context$context,"(?<=\\.)\\w")
ref<-str_extract(OD5mut_context$alteration,"\\w(?=\\w)")
alt<-str_extract(OD5mut_context$alteration, "(?<=\\w)\\w")
Mut_type<-paste0(Fiveprimecontext,"[",ref,">",alt,"]",threeprimecontext)
OD5mut_context$totpop<-"Rest"
OD5mut_context[Mut_type %in% probs$MutationTypes,]

ggplot(OD5mut_context, aes(totpop, AF))+geom_boxplot()+geom_jitter()
ggsave(plot=TCWplot, file="Rplotsgood/TCWplot.png",width=10,height=10)
grid<-plot_grid(mutrateplot+theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
          cosmat+theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                       axis.title.y = element_blank()),
          sigprofplot+theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
          TCWplot,
          ncol=1,nrow=4, align="v",axis="rl")
grid
save_plot( grid,file="Rplotsgood/Fig7.pdf",base_width=26,base_height=24)
