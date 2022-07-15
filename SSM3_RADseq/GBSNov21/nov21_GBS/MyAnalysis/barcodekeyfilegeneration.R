####fast_x_barcode_splitter only accepts barcodes of the same length as input, so this script splits up the barcodes file into those of different lengths


library(stringr)
barcodes<-read.delim("MyAnalysis/key/translationtable.txt", header = TRUE)
bc<-data.frame(barcodes$uidtag,as.character(barcodes$barcode))
seq<-c(0,rep(c(1,2,3), 15))
bc$fn<-paste0(bc$barcodes.uidtag, "_", seq)
cat(unlist(lapply(seq(nrow(bc)), function(x){
  paste0(">",bc[x,3], "\n","^",bc[x,2], "\n")
})), file = "MyAnalysis/bc.fasta")
write(levels(bc$barcodes.uidtag), file = "MyAnalysis/samples.txt")
