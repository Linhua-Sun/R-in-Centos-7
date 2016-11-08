library(ggpubr)
dataset<-read.table("6K_samples_meta_info.txt",header = T)

head(dataset)

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

p<-ggscatter(subset(pca,Group=="O.rufipogon" ),x ="PC1",y = "PC2",color = "Group",shape = "Group",size = 1,alpha=0.5)+theme_pubr(base_size = 16)

p+scale_color_manual(values = c(rainbow(n=6)))

ggscatter(pca,x ="PC2",y = "PC3",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times",base_size = 16)

ggscatter(pca,x ="PC3",y = "PC4",color = "Group",palette= "aaas",size = 1,alpha=0.7)+theme_pubr(base_family = "Times",base_size = 16)

a<-read.table("merged_K_all_chr_geno_0.1_maf_0.001.eval")


library(ggsci)
#################################################################################

p <- ggscatter(pca,x ="PC1",y = "PC2",palette = "npg",size = 1,alpha=0.3)

h1 <- p + geom_point(data = subset(pca,pca$Group == levels(pca$Group)[1] ),color=pal_npg()(6)[1],size=1,alpha=1)+ggtitle(label = levels(pca$Group)[1])
h2 <- p + geom_point(data = subset(pca,pca$Group == levels(pca$Group)[2] ),color=pal_npg()(6)[2],size=1,alpha=1)+ggtitle(label = levels(pca$Group)[2])
h3 <- p + geom_point(data = subset(pca,pca$Group == levels(pca$Group)[3] ),color=pal_npg()(6)[3],size=1,alpha=1)+ggtitle(label = levels(pca$Group)[3])
h4 <- p + geom_point(data = subset(pca,pca$Group == levels(pca$Group)[4] ),color=pal_npg()(6)[4],size=1,alpha=1)+ggtitle(label = levels(pca$Group)[4])
h5 <- p + geom_point(data = subset(pca,pca$Group == levels(pca$Group)[5] ),color=pal_npg()(6)[5],size=1,alpha=1)+ggtitle(label = levels(pca$Group)[5])
h6 <- p + geom_point(data = subset(pca,pca$Group == levels(pca$Group)[6] ),color=pal_npg()(6)[6],size=1,alpha=1)+ggtitle(label = levels(pca$Group)[6])

ggsave(filename = paste(levels(pca$Group)[1],"_PCA_points.png",sep=""),device = "png",plot = h1,width = 6,height = 4)
ggsave(filename = paste(levels(pca$Group)[2],"_PCA_points.png",sep=""),device = "png",plot = h2,width = 6,height = 4)
ggsave(filename = paste(levels(pca$Group)[3],"_PCA_points.png",sep=""),device = "png",plot = h3,width = 6,height = 4)
ggsave(filename = paste(levels(pca$Group)[4],"_PCA_points.png",sep=""),device = "png",plot = h4,width = 6,height = 4)
ggsave(filename = paste(levels(pca$Group)[5],"_PCA_points.png",sep=""),device = "png",plot = h5,width = 6,height = 4)
ggsave(filename = paste(levels(pca$Group)[6],"_PCA_points.png",sep=""),device = "png",plot = h6,width = 6,height = 4)


h<-multiplot(h1,h2,h3,h4,h5,h6,ncol=2)
show_col(pal_npg()(6))
pal_npg()(6)[1]

png(filename = "test_1.png")
grid.arrange(h1,h2,h3,h4,h5,h6,ncol=1)
dev.off()




fam<-read.table("merged_K_all_chr_geno_0.1_maf_0.001.fam",header = F)

length(fam$V2)

nrow(PCA)

unique(sort(dataset$VCF_ID %in% fam$V2))

##############################
dataset<-read.table("6K_samples_meta_info.txt",header = T)
pl<-read.table("plink.eigenvec",header = F)
head(pl)
attach(pl)
PCA<-data.frame(pl[,2:6])
detach(pl)
colnames(PCA)<-c("Taxa","PC1","PC2","PC3","PC4")
head(PCA)

num<-length(PCA$Taxa)
temp<-rep(x = "TEMP",num)
for (i in 1:num) {
    temp[i]<-as.character(dataset[match(PCA$Taxa[i],dataset$VCF_ID),]$Species)
}

attach(PCA)
pca<-data.frame(Taxa,temp,PC1,PC2,PC3,PC4)
detach(PCA)
colnames(pca)<-c("Taxa","Group","PC1","PC2","PC3","PC4")

tail(pca)

g1<-ggscatter(pca,x ="PC1",y = "PC2",color = "Group",palette = "aaas",size = 1,alpha=1)
g2<-ggscatter(pca,x ="PC2",y = "PC3",color = "Group",palette = "aaas",size = 1,alpha=1)
g3<-ggscatter(pca,x ="PC3",y = "PC4",color = "Group",palette = "aaas",size = 1,alpha=1)

ggsave("6K_g1.png",plot = g1,device = "png")

ggsave("6K_g2.png",plot = g2,device = "png")

ggsave("6K_g3.png",plot = g3,device = "png")

table(pca$Group)
