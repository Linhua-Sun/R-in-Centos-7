
library(rtracklayer)
te<-read.csv("/data1/linhua/QIANLAB/PROJECT/lnc_rna_list_ready.txt",header = T)
gr<-GRanges(seqnames = te$Chr,ranges = IRanges(start = te$Start,end = te$End,names = te$Loci_ID),strand = te$strand)
export.bed(object = gr,con = "test.bed")
export.gff3(object = gr,con = "test.gff3")

##
dir()
Ap<-import.gff3(con = "/data1/linhua/Rstudio/Araport11_GFF3_genes_transposons.201606.gff")

setwd("/data1/linhua/QIANLAB/PROJECT/Long-Noncoding-RNA-project/QC_RNA_SEQ/batch_com/XLS")
dir()

tophat_unstranded<-read.table("Tophat_unstrand_Col-0.tin.xls",header=T)
tophat_637_stranded<-read.table("tophat_637_stranded.tin.xls",header = T)
star_637_stranded<-read.table("star_637_stranded.tin.xls",header = T)


all<-rbind(tophat_unstranded,tophat_637_stranded,star_637_stranded)
head(all)
tail(all)
library()

boxplot(tophat_unstranded$TIN,tophat_637_stranded$TIN,star_637_stranded$TIN)

head(t1)
library(ggpubr)
ggdensity(data = t1,x = t1$TIN)
?ggdensity


all<-read.table("merge.txt",header = F)
head(all)
colnames(all)<-c("group","gene_ID","chr","start","end","tin")
head(all)
ggviolin(data = all,x = "",y = all$tin)

# Load data
data("ToothGrowth")
df <- ToothGrowth
ggviolin(df, "dose", "len",  color = "dose",
         palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot")

head(df)


