
library(rtracklayer)
te<-read.csv("/data1/linhua/QIANLAB/PROJECT/lnc_rna_list_ready.txt",header = T)
gr<-GRanges(seqnames = te$Chr,ranges = IRanges(start = te$Start,end = te$End,names = te$Loci_ID),strand = te$strand)
export.bed(object = gr,con = "test.bed")
export.gff3(object = gr,con = "test.gff3")

