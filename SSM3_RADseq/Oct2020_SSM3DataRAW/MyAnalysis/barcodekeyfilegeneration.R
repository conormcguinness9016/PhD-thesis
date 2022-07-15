library(stringr)
barcodes<-read.delim("key.txt", header = TRUE)
barcodes$V1<-as.character(barcodes$V1)



bc<-as.matrix(barcodes)
bc<-data.frame(v1=paste(">",barcodes[,4], sep = ""), v2=paste("^", barcodes[,3], sep = ""))
bc<-paste(bc[,1], bc[,2], sep = "\n")
cat(as.matrix(bc), sep = "\n", file = "bc.fasta")


