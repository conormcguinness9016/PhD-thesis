

library(CytoExploreR)
gs <- cyto_setup("./comps/",
                 gatingTemplate = "Activation-gatingTemplate.csv")
gs <- cyto_transform(gs,
                     type = "logicle")

##transforms the fluorescent gates to logicle format

cyto_gatingTemplate_apply(gs)
##draw the FSC-H and SSC-H plot
cyto_gate_draw(gs,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"))

cyto_gate_draw(gs,
               parent = "Cells",
               alias = "Single Cells",
               channels = c("FSC-A","FSC-H"))
spill <- cyto_spillover_compute(gs,
                                parent = "Single Cells",
                                spillover = "Spillover-Matrix.csv")
cyto_plot_compensation(gs,
                       parent = "Single Cells",
                       spillover = "Spillover-Matrix.csv",
                       compensate = TRUE)
cyto_plot_compensation(gs,
                       parent = "Single Cells")
rm(gs)
setwd("fcsfiles/")
gs <- cyto_setup(".",
                 gatingTemplate = "Activation-gatingTemplate.csv")
gs <- cyto_compensate(gs,
                      spillover = "../Spillover-Matrix.csv")
gs <- cyto_transform(gs,
                     type = "logicle")
cyto_details(gs)$mutant[cyto_details(gs)$mutant == "NA"]<-""

cyto_details(gs)$Enzyme<-paste0(cyto_details(gs)$enzyme, cyto_details(gs)$mutant)
cyto_details(gs)$EnzymeUGI<-paste0(cyto_details(gs)$Enzyme,"_", cyto_details(gs)$UGI)
cyto_gatingTemplate_apply(gs)
cyto_gate_draw(gs,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"))

cyto_gate_draw(gs,
               parent = "Cells",
               alias = "Single Cells",
               channels = c("FSC-A","FSC-H"))
cyto_names(gs)
LDFMO<-cyto_extract(gs, parent = "Single Cells")[c(17,29)]
cyto_gate_draw(gs,
               parent = "Single Cells",
               alias = c("Dead cells", "Live cells"),
               channels = c("R_780/60-A"),
               overlay = LDFMO,
               type = "interval")
Unst<-cyto_extract(gs, parent = "Single Cells")[29]
GFP<-cyto_extract(gs, parent = "Single Cells")[21]
mCherry<-cyto_extract(gs, parent = "Single Cells")[22]
DT<-cyto_extract(gs, parent = "Single Cells")[19]
cyto_names(gs)
cyto_gate_edit(gs, parent = "Live cells",
               alias = c("mCherry+","mCherry+GFP+"),
               channels = c("YG_610/20-A", "B_530/30-A"),
               overlay = c(DT),
               density_fill_alpha=0.5,
               type = c("rectangle", "polygon"))
cyto_gate_remove(gs, parent = "mCherry+",
                 alias = "GFP+",
                 channels = c("B_530/30-A"),
                 overlay = c(Unst, GFP),
                 density_fill_alpha=0.5,
                 type = "threshold")


cyto_names(gs)
##cyto_extract "pulls out" one sample
cyto_details_edit(gs)
##LD gate draw
cyto_details(gs)$Plotorder
plotorder<-order(na.omit(as.numeric(cyto_details(gs)$plotorder)))
plotorder
labels<-cyto_details(gs)$Enzyme[rev(plotorder)]
cyto_plot(gs[rev(plotorder)],
          parent = "Single Cells",
          alias = "Live cells",
          channels = c("R_780/60-A"),
          #select = list(Sample=c("Yes")),
          density_modal = FALSE, 
          legend = TRUE,
          legend_text =  labels,
          density_fill_alpha = 1,
          xlab = "Zombie NIRdye Live Dead Stain R_780/60-A",
          density_stack = 1.3, popup = TRUE)
names<-paste0(rev(labels), cyto_details(gs)$UGI[plotorder])
cyto_plot_save("../figs/SourApplesplot.png", width =15, height = 10)
cyto_plot(gs[plotorder],
          parent = "Live cells",
          alias = "",
          channels = c("YG_610/20-A", "B_530/30-A"),
          # select = list(Comp=c("No")),
          density_modal = FALSE,
          # legend_text =  legendtextLD,
          density_fill_alpha = 1,
          xlab = "mCherry",
          ylab = "GFP",
          title = names,
          density_stack = 1.3, popup = TRUE)
dev.off()
cyto_plot_save("../Rplots/mCherry+plot.png", width =15)
cyto_plot(gs[rev(plotorder)],
          parent = "mCherry+",
          alias = "",
          channels = c("B_530/30-A"),
          #select = list(Comp=c("No")),
          density_modal = FALSE,
          legend_text =  rev(names),
          density_fill_alpha = 1,
          xlab = "GFP",
          legend = TRUE,
          density_stack = 1.3, popup = TRUE, label_position = "manual")
cyto_gate_edit(gs, parent = "mCherry+",
               alias = c("GFP+"),
               channels = c("B_530/30-A"),
               overlay = c(Unst, DT),
               density_fill_alpha=0.5,
               type = c("interval"))
cyto_plot(gs,
          parent = "mCherry+",
          alias = "",
          channels = c("B_530/30-A"),
          #select = list(Sample=c("Yes")),
          density_modal = FALSE, 
          legend = TRUE,
          #legend_text =  labels,
          density_fill_alpha = 1,
          xlab = "GFP+",
          density_stack = 1.3, popup = TRUE)
cyto_plot(gs,alias="",parent="root",channels = c("FSC-A","SSC-A"), group_by = "all")
cyto_plot(gs,alias="",parent="Cells",channels = c("FSC-A","FSC-H"), group_by = "all")
cyto_plot(gs,alias="",parent="Cells",channels = c("FSC-A","FSC-H"), group_by = "all")
cyto_plot(gs,alias="",parent="Single Cells",channels = c("R_780/60-A"), group_by = "all")
cyto_plot(gs,alias="",parent="Live cells",channels = c("YG_610/20-A", "B_530/30-A"), group_by = "EnzymeUGI")
library(tidyverse)

GFPstats<-cyto_stats_compute(gs, channels = "B_530/30-A", stat = c("median"),alias = "GFP+", parent = "mCherry+")


GFPfreqs<-cyto_stats_compute(gs, channels = "B_530/30-A", stat = c("freq"), alias = "GFP+", parent = "mCherry+")


LivecellfreqsmChez<-cyto_stats_compute(gs, channels = c("R_780/60-A"), 
                                       stat = c("freq"), 
                                       alias = "Live cells", parent = "Single Cells")
GFPfreqsmChez<-cyto_stats_compute(gs, channels = c("YG_610/20-A", "B_530/30-A"), 
                                  stat = c("freq"), 
                                  alias = c("mCherry+GFP+","mCherry+"), parent = "Live cells")
mCherryfreqs<-cyto_stats_compute(gs, channels = c("YG_610/20-A", "B_530/30-A"), 
                                 stat = c("freq"), 
                                 alias = "mCherry+", parent = "Live cells")
saveRDS(LivecellfreqsmChez, "../LivecellfreqsmChez.RDS")
saveRDS(GFPstats, "../GFPstats.RDS")
saveRDS(GFPfreqsmChez, "../GFPfreqsmChez.RDS")
saveRDS(GFPfreqs, "../GFPfreqs.RDS")

png("../figs/Livecellsnormalised.png")
ggplot(LivecellfreqsmChez[!LivecellfreqsmChez$Enzyme == "NA" & 
                            !LivecellfreqsmChez$Enzyme == "None"& 
                            !LivecellfreqsmChez$Enzyme == "DT",], aes(Enzyme, colour = UGI, Frequency)) + 
  geom_point(aes(group = UGI),  position = position_jitterdodge(jitter.width = .1, dodge.width = 1))+ geom_boxplot(alpha = 0.5)+ylab("Normalised % live cells")
dev.off()



png("../figs/cornormalisedLivecells.png")
plot(LivecellfreqsmChez$Frequency[!LivecellfreqsmChez$Enzyme == c("NA") &
                                    !LivecellfreqsmChez$Enzyme == "DMSO" &
                                    !LivecellfreqsmChez$Enzyme == "PMA"], GFPfreqsmChez$Frequency[!GFPfreqsmChez$Enzyme == c("NA") &
                                                                                                    !GFPfreqsmChez$Enzyme == "DMSO" &
                                                                                                    !GFPfreqsmChez$Enzyme == "PMA"], xlab= "Normalised percentage live cells", ylab="Percentage mCherry+GFP+ cells (live cell population)")

dev.off()
png("../figs/corLivecells.png")
plot(LivecellfreqsmChez$Frequency[!LivecellfreqsmChez$Enzyme == c("NA") &
                                    !LivecellfreqsmChez$Enzyme == "DMSO" &
                                    !LivecellfreqsmChez$Enzyme == "PMA"], GFPfreqsmChez$Frequency[!GFPfreqsmChez$Enzyme == c("NA") &
                                                                                                    !GFPfreqsmChez$Enzyme == "DMSO" &
                                                                                                    !GFPfreqsmChez$Enzyme == "PMA"], xlab= "Percentage live cells", ylab="Percentage mCherry+GFP+ cells (live cell population)")

dev.off()
png("../figs/GFPpercentmCherrycellpop.png")
ggplot(GFPfreqs[!GFPfreqs$Enzyme == "NA" & !GFPfreqs$Enzyme == "None",], aes(Enzyme, colour = UGI,Frequency)) + 
  geom_point(aes(group = UGI),  position = position_jitterdodge(jitter.width = .1, dodge.width = 1))+ geom_boxplot(alpha = 0.5)+ylab(" %GFP+ cells (out of mCherry cells)")
dev.off()

png("../figs/GFPpercentlivecellpop.png")
ggplot(GFPfreqsmChez[!GFPfreqsmChez$Enzyme == c("NA"),], aes(Enzyme, colour = UGI, Frequency)) + 
  geom_point(aes(group = UGI),  position = position_jitterdodge(jitter.width = .1, dodge.width = 1))+ geom_boxplot(alpha = 0.5)+ylab(" %mCherry+ GFP+ cells (out of live cells)")
dev.off()
png("../figs/mCherryfreqs.png")
ggplot(mCherryfreqs[!GFPfreqsmChez$Enzyme == c("NA"),], aes(Enzyme, colour = UGI, Frequency)) + 
  geom_point(aes(group = UGI),  position = position_jitterdodge(jitter.width = .1, dodge.width = 1))+ geom_boxplot(alpha = 0.5)+ylab(" %mCherry+ cells (out of live cells)")
dev.off()
cyto_details_edit(gs)
