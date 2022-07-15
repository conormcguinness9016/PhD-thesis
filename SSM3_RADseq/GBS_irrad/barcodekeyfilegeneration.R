####fast_x_barcode_splitter only accepts barcodes of the same length as input, so this script splits up the barcodes file into those of different lengths


library(stringr)
barcodes<-read.delim("key_set1.txt", header = TRUE)
bc<-data.frame(barcodes$sample,barcodes$barcode)
bc$barcodes.sample<-str_replace(bc$barcodes.sample, " ", "_")
bc$sample<-str_extract(bc$barcodes.sample, "\\w{0,5}(?=_)")
seq<-c(rep(c(1,2), 15))
bc$fn<-paste0(bc$sample, "_", seq)
cat(unlist(lapply(seq(nrow(bc)), function(x){
  paste0(">",bc[x,4], "\n","^",bc[x,2], "\n")
})), file = "bc.fasta")
write(unique(bc$sample), file = "samples.txt")
