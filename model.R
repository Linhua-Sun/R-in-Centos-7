## mdoel 
dataset<-read.csv("/data1/linhua/Hang_Rice_Data/All_5K_data_2016-10-24/LD_OKED/5152_meta_info_2016_10_21.csv",header = T)

Osa<-subset(dataset,dataset$Species=="O.sativa")

head(Osa)

Osa<-subset(dataset,dataset$Subpopulation_1)

table(dataset$Subpopulation_1)

osa<-subset(dataset,dataset$Species=="O.sativa")
ind<-subset(dataset,dataset$Subpopulation_1=="Indica")
jap<-subset(dataset,dataset$Subpopulation_1=="japonica")
rufi<-subset(dataset,dataset$Subpopulation_1=="O.rufipogon")

write(as.character(osa$Sample_ID_IN_VCF),file = "osa_sample_names.txt")
write(as.character(ind$Sample_ID_IN_VCF),file = "ind_sample_names.txt")
write(as.character(jap$Sample_ID_IN_VCF),file = "jap_sample_names.txt")
write(paste("0",as.character(rufi$Sample_ID_IN_VCF)),file = "rufi_sample_names.txt")


## test

library(adegenet)
rufi_raw<-read.PLINK("rufi.raw")


#!/usr/bin/env Rscript
#######
Args<- commandArgs(trailingOnly = TRUE)
input<-Args[1]

library(adegenet)

rufi_raw<-read.PLINK(input)
png(filename = paste(input,".png",sep = ""))
plot(rufi_raw)
dev.off()



#######
all<-read.table("merged_5K_all_chr_missing.imiss",header = T)
plot(density.default(all$N_MISS),xlab="Number of missing genotype",ylab="Density",main="",lty=1,lwd=3)
abline(v = 33559583,lwd=2,col="red")
text(x=33559583,y=4e-08,adj =1,labels = "Only 10^6 SNP",cex = 1,col="blue",lwd=2)

subset(all,all$N_GENO > 34459583)

34559583-100000



#####  


NUM<-read.table("SNP_NUM_SAMPLES.txt")
head(NUM)



plot(density.default(NUM$V2),xlab="Number of called SNPs",ylab="Density",main="",lty=1,lwd=3)
abline(v = 10^5,lwd=2,col="red")
text(x=10^5,y=4e-08,adj =0,labels = "Only 10^5 SNP",cex = 1,col="blue",lwd=2)

filtered_samples<-subset(NUM,NUM$V2 < 10^5)

head(filtered_samples)
write(paste("0",as.character(filtered_samples$V1)),file = "filtered_samples_names_10_5.txt")

#####



af_NUM<-setdiff(NUM$V1,filtered_samples$V1)

head(af_NUM)
head(dataset)
colnames(dataset)


osa<-subset(dataset,dataset$Species=="O.sativa" & dataset$Sample_ID_IN_VCF %in% af_NUM)

ind<-subset(dataset,dataset$Subpopulation_1=="Indica" & dataset$Sample_ID_IN_VCF %in% af_NUM)

jap<-subset(dataset,dataset$Subpopulation_1=="japonica" & dataset$Sample_ID_IN_VCF %in% af_NUM)

rufi<-subset(dataset,dataset$Subpopulation_1=="O.rufipogon" & dataset$Sample_ID_IN_VCF %in% af_NUM)


write(paste("0",as.character(osa$Sample_ID_IN_VCF)),file =  "filtered_osa_sample_names.txt")
write(paste("0",as.character(ind$Sample_ID_IN_VCF)),file =  "filtered_ind_sample_names.txt")
write(paste("0",as.character(jap$Sample_ID_IN_VCF)),file =  "filtered_jap_sample_names.txt")
write(paste("0",as.character(rufi$Sample_ID_IN_VCF)),file = "filtered_rufi_sample_names.txt")

O.nivara <- subset(dataset,dataset$Species=="O.nivara" & dataset$Sample_ID_IN_VCF %in% af_NUM)
write(paste("0",as.character(O.nivara$Sample_ID_IN_VCF)),file =  "filtered_onivara_sample_names.txt")


###################################################################################################
## analysis rufi_PCA
library(ggpubr)

## input subset_rufi_samples.pca.evec

dataset<-read.csv("/data1/linhua/Hang_Rice_Data/All_5K_data_2016-10-24/LD_OKED/5152_meta_info_2016_10_21.csv",header = T)

PCA<-read.table("subset_rufi_samples.pca.evec",header = F,comment.char = "#")
PCA$V12<-NULL
colnames(PCA)<-c("Taxa","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
PCA$Taxa<-gsub(PCA$Taxa,pattern = "0:",replacement = "")
num<-length(PCA$Taxa)
temp<-rep(x = "TEMP",num)
for (i in 1:num) {
    temp[i]<-as.character(dataset[match(PCA$Taxa[i],dataset$Sample_ID_IN_VCF),]$Subpopulation)
}
attach(PCA)
pca<-data.frame(Taxa,temp,PC1,PC2,PC3,PC4,PC5)
detach(PCA)
colnames(pca)<-c("Taxa","Group","PC1","PC2","PC3","PC4","PC5")

ggscatter(pca,x ="PC1",y = "PC2",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)

ggscatter(pca,x ="PC2",y = "PC3",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)

ggscatter(pca,x ="PC3",y = "PC4",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)

## input
## without_maf_geno_subset_rufi_samples.pca.evec

dataset<-read.csv("/data1/linhua/Hang_Rice_Data/All_5K_data_2016-10-24/LD_OKED/5152_meta_info_2016_10_21.csv",header = T)
PCA<-read.table("without_maf_geno_subset_rufi_samples.pca.evec",header = F,comment.char = "#")
PCA$V12<-NULL
colnames(PCA)<-c("Taxa","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
PCA$Taxa<-gsub(PCA$Taxa,pattern = "0:",replacement = "")
num<-length(PCA$Taxa)
temp<-rep(x = "TEMP",num)
for (i in 1:num) {
    temp[i]<-as.character(dataset[match(PCA$Taxa[i],dataset$Sample_ID_IN_VCF),]$Subpopulation)
}
attach(PCA)
pca<-data.frame(Taxa,temp,PC1,PC2,PC3,PC4,PC5)
detach(PCA)
colnames(pca)<-c("Taxa","Group","PC1","PC2","PC3","PC4","PC5")
ggscatter(pca,x ="PC1",y = "PC2",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)

ggscatter(pca,x ="PC2",y = "PC3",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)

ggscatter(pca,x ="PC3",y = "PC4",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)


## input 

## subset_osa__samples.pca.evec
dataset<-read.csv("/data1/linhua/Hang_Rice_Data/All_5K_data_2016-10-24/LD_OKED/5152_meta_info_2016_10_21.csv",header = T)
PCA<-read.table("subset_osa__samples.pca.evec",header = F,comment.char = "#")
PCA$V12<-NULL
colnames(PCA)<-c("Taxa","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
PCA$Taxa<-gsub(PCA$Taxa,pattern = "0:",replacement = "")
num<-length(PCA$Taxa)
temp<-rep(x = "TEMP",num)
for (i in 1:num) {
    temp[i]<-as.character(dataset[match(PCA$Taxa[i],dataset$Sample_ID_IN_VCF),]$Subpopulation_1)
}
attach(PCA)
pca<-data.frame(Taxa,temp,PC1,PC2,PC3,PC4,PC5)
detach(PCA)
colnames(pca)<-c("Taxa","Group","PC1","PC2","PC3","PC4","PC5")
ggscatter(pca,x ="PC1",y = "PC2",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)

ggscatter(pca,x ="PC2",y = "PC3",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)

ggscatter(pca,x ="PC3",y = "PC4",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)

## subset_ind__samples.pca.evec
dataset<-read.csv("/data1/linhua/Hang_Rice_Data/All_5K_data_2016-10-24/LD_OKED/5152_meta_info_2016_10_21.csv",header = T)
PCA<-read.table("subset_ind__samples.pca.evec",header = F,comment.char = "#")
PCA$V12<-NULL
colnames(PCA)<-c("Taxa","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
PCA$Taxa<-gsub(PCA$Taxa,pattern = "0:",replacement = "")
num<-length(PCA$Taxa)
temp<-rep(x = "TEMP",num)
for (i in 1:num) {
    temp[i]<-as.character(dataset[match(PCA$Taxa[i],dataset$Sample_ID_IN_VCF),]$Subpopulation_1)
}
attach(PCA)
pca<-data.frame(Taxa,temp,PC1,PC2,PC3,PC4,PC5)
detach(PCA)
colnames(pca)<-c("Taxa","Group","PC1","PC2","PC3","PC4","PC5")
ggscatter(pca,x ="PC1",y = "PC2",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)

ggscatter(pca,x ="PC2",y = "PC3",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)

ggscatter(pca,x ="PC3",y = "PC4",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)
## subset_jap__samples.pca.evec

dataset<-read.csv("/data1/linhua/Hang_Rice_Data/All_5K_data_2016-10-24/LD_OKED/5152_meta_info_2016_10_21.csv",header = T)
PCA<-read.table("subset_jap__samples.pca.evec",header = F,comment.char = "#")
PCA$V12<-NULL
colnames(PCA)<-c("Taxa","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
PCA$Taxa<-gsub(PCA$Taxa,pattern = "0:",replacement = "")
num<-length(PCA$Taxa)
temp<-rep(x = "TEMP",num)
for (i in 1:num) {
    temp[i]<-as.character(dataset[match(PCA$Taxa[i],dataset$Sample_ID_IN_VCF),]$Subpopulation_1)
}
attach(PCA)
pca<-data.frame(Taxa,temp,PC1,PC2,PC3,PC4,PC5)
detach(PCA)
colnames(pca)<-c("Taxa","Group","PC1","PC2","PC3","PC4","PC5")
ggscatter(pca,x ="PC1",y = "PC2",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)

ggscatter(pca,x ="PC2",y = "PC3",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)

ggscatter(pca,x ="PC3",y = "PC4",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)


## merged_5K_all_chr_geno_0.1_maf_0.001_prunedData.pca.evec

dataset<-read.csv("/data1/linhua/Hang_Rice_Data/All_5K_data_2016-10-24/LD_OKED/5152_meta_info_2016_10_21.csv",header = T)
PCA<-read.table("merged_5K_all_chr_geno_0.1_maf_0.001_prunedData.pca.evec",header = F,comment.char = "#")
PCA$V12<-NULL
colnames(PCA)<-c("Taxa","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
PCA$Taxa<-gsub(PCA$Taxa,pattern = "0:",replacement = "")
num<-length(PCA$Taxa)
temp<-rep(x = "TEMP",num)
for (i in 1:num) {
    temp[i]<-as.character(dataset[match(PCA$Taxa[i],dataset$Sample_ID_IN_VCF),]$Species)
}
attach(PCA)
pca<-data.frame(Taxa,temp,PC1,PC2,PC3,PC4,PC5)
detach(PCA)
colnames(pca)<-c("Taxa","Group","PC1","PC2","PC3","PC4","PC5")
ggscatter(pca,x ="PC1",y = "PC2",color = "Group",palette= "lancet",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)

ggscatter(pca,x ="PC2",y = "PC3",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)

ggscatter(pca,x ="PC3",y = "PC4",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times New Roman",base_size = 16)


### subset_rice_samples.pca.evec

dataset<-read.csv("/data1/linhua/Hang_Rice_Data/All_5K_data_2016-10-24/LD_OKED/5152_meta_info_2016_10_21.csv",header = T)

length(dataset$Papers._ID)
test<-rep(x = "TEMP",5152)

for (i in 1:5152) {
    if (dataset$Subpopulation_1[i] == "Indica")   {test[i]<-"indica"} 
    if (dataset$Subpopulation_1[i] == "japonica") {test[i]<-"japonica"} 
    if (dataset$Subpopulation_1[i] == "O.rufipogon") {test[i]<-"O.rufipogon"}
    if (dataset$Species[i] == "O.nivara") {test[i]<-"O.nivara"}
    if (dataset$Species[i] == "O.barthii") {test[i]<-"O.barthii"}
    if (dataset$Species[i] == "O.Glaberrima") {test[i]<-"O.Glaberrima"}
    if (dataset$Species[i] == "O.longistaminata") {test[i]<-"O.longistaminata"}
    if (dataset$Species[i] == "O.meridionalis") {test[i]<-"O.meridionalis"}
}
for (i in 1:5152) {
    if (test[i]=="TEMP") {test[i]<-"admixed_samples"}
}

New_dataset<-data.frame(dataset,test)
colnames(New_dataset)[9]="group"

PCA<-read.table("subset_rice_samples.pca.evec",header = F,comment.char = "#")
PCA$V12<-NULL
colnames(PCA)<-c("Taxa","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
PCA$Taxa<-gsub(PCA$Taxa,pattern = "0:",replacement = "")
num<-length(PCA$Taxa)
temp<-rep(x = "TEMP",num)
for (i in 1:num) {
    temp[i]<-as.character(New_dataset[match(PCA$Taxa[i],New_dataset$Sample_ID_IN_VCF),]$group)
}
attach(PCA)
pca<-data.frame(Taxa,temp,PC1,PC2,PC3,PC4,PC5)
detach(PCA)
colnames(pca)<-c("Taxa","Group","PC1","PC2","PC3","PC4","PC5")

ggscatter(pca,x ="PC1",y = "PC2",color = "Group",palette= "lancet",size = 1,alpha=0.7)+theme_pubr(base_family = "Times",base_size = 16)
ggsave(filename = "subset_rice_samples.pca.evec_pc1_pc2.pdf",device = "pdf")
ggscatter(pca,x ="PC2",y = "PC3",color = "Group",palette= "lancet",size = 1,alpha=0.7)+theme_pubr(base_family = "Times",base_size = 16)
ggsave(filename = "subset_rice_samples.pca.evec_pc2_pc3.pdf",device = "pdf")
ggscatter(pca,x ="PC3",y = "PC4",color = "Group",palette= "lancet",size = 1,alpha=0.7)+theme_pubr(base_family = "Times",base_size = 16)
ggsave(filename = "subset_rice_samples.pca.evec_pc3_pc4.pdf",device = "pdf")
