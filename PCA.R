#######
### analysis 5152 samples PCA tree and admixtures

## PCA
setwd("/data1/linhua/Hang_Rice_Data/All_5K_data_2016-10-24/LD_OKED")

dataset<-read.csv("/data1/linhua/Hang_Rice_Data/All_5K_data_2016-10-24/LD_OKED/5152_meta_info_2016_10_21.csv",header = T)

PCA<-read.table("merged_5K_all_chr_geno_0.1_maf_0.001_prunedData_pca.eigenvec",header = F)

PCA$V1<-NULL

colnames(PCA)<-c("Taxa","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")

num<-length(PCA$Taxa)

temp<-rep(x = "TEMP",num)

for (i in 1:num) {
    temp[i]<-as.character(dataset[match(PCA$Taxa[i],dataset$Sample_ID_IN_VCF),]$Subpopulation_1)
}

head(PCA)
attach(PCA)
pca<-data.frame(Taxa,temp,PC1,PC2,PC3,PC4,PC5)
detach(PCA)
colnames(pca)<-c("Taxa","Group","PC1","PC2","PC3","PC4","PC5")
head(pca,n = 10)
#############################################################################

plot(x = pca$PC1,y=pca$PC2)
plot(x=pca$PC2,y=pca$PC3)
ggscatter(pca,x ="PC1",y = "PC2",color = "Group",palette= "ucscgb",size = 4,alpha=0.6,label = "Taxa",repel = T)



## analysis missing data situation

imiss<-read.table("merged_5K_all_chr_geno_0.1_maf_0.001_prunedData.imiss")