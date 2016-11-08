####

setwd("/data1/linhua/QIANLAB/PROJECT/lncRNA_annotation_merge/")
dir()
##

library(rtracklayer)

NONCODE_lnc<-import.bed("NONCODE2016_tair10.lncAndGene.bed")
head(NONCODE_lnc)
levels(NONCODE_lnc@seqnames)

zhu<-import.gff3("zhu_im-ncRNA.gff3")
levels(zhu@seqnames)

levels(NONCODE_lnc@seqnames)<-c("Chr1","Chr2","Chr3","Chr4","Chr5" ,"ChrM" ,"ChrC")
head(NONCODE_lnc)

length(unique(sort(substr(NONCODE_lnc@elementMetadata@listData$name,start = 1,stop = 13)))) == length(unique(sort(NONCODE_lnc@elementMetadata@listData$name)))


new<-GRanges(seqnames = NONCODE_lnc@seqnames,IRanges(start = NONCODE_lnc@ranges@start,width=NONCODE_lnc@ranges@width),strand = NONCODE_lnc@strand,ID=NONCODE_lnc@elementMetadata@listData$name)

export.gff3(new,con="Chr_NONCODE2016_tair10.lncAndGene.gff")

head(new)

plot(density.default(new@ranges@width))


######################################

dir()

Araport_extracted_antisense_lncRNA <-import.gff3("Araport_extracted_antisense_lncRNA.gff3")
Araport_extracted_antisense_RNA    <-import("Araport_extracted_antisense_RNA.gff3")
Araport_extracted_lnc_RNA          <-import("Araport_extracted_lnc_RNA.gff3")
Araport_extracted_ncRNA            <-import("Araport_extracted_ncRNA.gff3")
Chr_NONCODE2016_tair10             <-import("Chr_NONCODE2016_tair10.lncAndGene.gff")
zhu_im_ncRNA                       <-import("zhu_im-ncRNA.gff3")


Araport_extracted_antisense_lncRNA_width<-Araport_extracted_antisense_lncRNA@ranges@width 
Araport_extracted_antisense_RNA_width   <-Araport_extracted_antisense_RNA   @ranges@width
Araport_extracted_lnc_RNA_width         <-Araport_extracted_lnc_RNA         @ranges@width
Araport_extracted_ncRNA_width           <-Araport_extracted_ncRNA           @ranges@width
Chr_NONCODE2016_tair10_width            <-Chr_NONCODE2016_tair10            @ranges@width
zhu_im_ncRNA_width                      <-zhu_im_ncRNA                      @ranges@width


length(Araport_extracted_antisense_lncRNA_width)<-6330
length(Araport_extracted_antisense_RNA_width   )<-6330
length(Araport_extracted_lnc_RNA_width         )<-6330
length(Araport_extracted_ncRNA_width           )<-6330
length(Chr_NONCODE2016_tair10_width            )<-6330
length(zhu_im_ncRNA_width                      )<-6330

annotation<-data.frame(ID=seq(1,6330),Araport_extracted_antisense_lncRNA_width,Araport_extracted_antisense_RNA_width,Araport_extracted_lnc_RNA_width,Araport_extracted_ncRNA_width,Chr_NONCODE2016_tair10_width,zhu_im_ncRNA_width)
colnames(annotation)<-gsub(colnames(annotation),pattern = "_width",replacement = "")
colnames(annotation)<-gsub(colnames(annotation),pattern = "_extracted",replacement = "")
colnames(annotation)<-gsub(colnames(annotation),pattern = "Chr_",replacement = "")

library(reshape2)
M_annotation<-melt(data = annotation,variable.name = "Group",id.vars = "ID")
head(M_annotation)
table(M_annotation$Group)
ggviolin(data = M_annotation,x = "Group",y = "value")


#####################################################################

export.gff3(Araport_extracted_antisense_lncRNA ,con = "modi_Araport_extracted_antisense_lncRNA.gff3")
export.gff3(Araport_extracted_antisense_RNA    ,con = "modi_Araport_extracted_antisense_RNA.gff3")
export.gff3(Araport_extracted_lnc_RNA          ,con = "modi_Araport_extracted_lnc_RNA.gff3")
export.gff3(Araport_extracted_ncRNA            ,con = "modi_Araport_extracted_ncRNA.gff3")
export.gff3(Chr_NONCODE2016_tair10             ,con = "modi_Chr_NONCODE2016_tair10.gff3")
export.gff3(zhu_im_ncRNA                       ,con = "modi_zhu_im_ncRNA.gff3")

#### bed ####

export.bed(Araport_extracted_antisense_lncRNA,con = "modi_Araport_extracted_antisense_lncRNA.bed")
export.bed(Araport_extracted_antisense_RNA   ,con = "modi_Araport_extracted_antisense_RNA.bed")
export.bed(Araport_extracted_lnc_RNA         ,con = "modi_Araport_extracted_lnc_RNA.bed")
export.bed(Araport_extracted_ncRNA           ,con = "modi_Araport_extracted_ncRNA.bed")
export.bed(Chr_NONCODE2016_tair10            ,con = "modi_Chr_NONCODE2016_tair10.bed")
export.bed(zhu_im_ncRNA                      ,con = "modi_zhu_im_ncRNA.bed")
