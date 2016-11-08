## analysis pca Tue Nov  8 09:55:59 2016

## 

setwd("/data1/linhua/Hang_Rice_Data/6K_DOWNLOAD_DATA/2016-11-07-download")
dir()

k1<-read.table("new_1K.txt",header = F,sep = "\t")
head(k1)

fam<-read.table("merged_K_all_chr_geno_0.1_maf_0.001_prunedData.fam",header = F,sep = " ")
head(fam)
k1$V1 %in% fam$V2


head(k1)
table(k1$V3)

levels(k1$V3)[8]<-"weedy_rice"
levels(k1$V3)[7]<-"japonica"
levels(k1$V3)[6]<-"indica"
head(k1)
table(k1$V3)
colnames(k1)<-c("VCF_ID","Species","subspecies")
write.table(k1,file = "Last_1K_samples_info.txt",quote = F,col.names = T,row.names = F,sep = "\t")


#######################################

K1<-read.table("Last_1K_samples_info.txt",header = T,sep = "\t")

head(K1)
table(K1$subspecies)
table(K1$Species)
levels(K1$Species)<-c("O.glaberrima","O.barthii","O.rufipogon","O.sativa")

New_K1<-K1[order(K1$Species),]

write.table(New_K1,file = "Last_1K_samples_info_1.txt",quote = F,col.names = T,row.names = F,sep = "\t")


K1<-read.table("Last_1K_samples_info_1.txt",header = T)

head(K1)



dataset<-read.csv("/data1/linhua/Hang_Rice_Data/All_5K_data_2016-10-24/LD_OKED/5152_meta_info_2016_10_21.csv",header = T)

K5<-data.frame(dataset$Sample_ID_IN_VCF,dataset$Species,dataset$Subpopulation_1)


colnames(K5)<-colnames(K1)

head(K5)
head(K1)

K6<-rbind(K5,K1)
colnames(K6)<-c("VCF_ID","Species","Subspecies")

K6[order(K6$Species,K6$Subspecies,K6$VCF_ID),]



#########################################################################################





library(ggpubr)
dataset<-read.table("6K_samples_meta_info.txt",header = T)

PCA<-read.table("merged_K_all_chr_geno_0.1_maf_0.001.pca.evec",header = F,comment.char = "#")
PCA$V12<-NULL
colnames(PCA)<-c("Taxa","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
PCA$Taxa<-gsub(PCA$Taxa,pattern = "0:",replacement = "")
head(PCA)

num<-length(PCA$Taxa)
temp<-rep(x = "TEMP",num)
for (i in 1:num) {
    temp[i]<-as.character(dataset[match(PCA$Taxa[i],dataset$VCF_ID),]$Species)
}
attach(PCA)
pca<-data.frame(Taxa,temp,PC1,PC2,PC3,PC4,PC5)
detach(PCA)
colnames(pca)<-c("Taxa","Group","PC1","PC2","PC3","PC4","PC5")

head(pca)

p<-ggscatter(pca,x ="PC1",y = "PC2",color = "Group",shape = "Group",size = 1,alpha=0.5)+theme_pubr(base_size = 16)

p+scale_color_manual(values = c(rainbow(n=6)))

ggscatter(pca,x ="PC2",y = "PC3",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times",base_size = 16)

ggscatter(pca,x ="PC3",y = "PC4",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times",base_size = 16)

a<-read.table("merged_K_all_chr_geno_0.1_maf_0.001.eval")

