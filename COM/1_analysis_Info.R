## important project in analysis the data from bcftools

library(data.table)

colnames_info<-c("CHROM","POS","REF,ALT","AC","AF","AN","BaseQRankSum","ClippingRankSum","DP","DS","END","ExcessHet","FS","HaplotypeScore","InbreedingCoeff","MLEAC","MLEAF","MQ","MQRankSum","QD","RAW_MQ","ReadPosRankSum","SOR")

info_01<-fread("/data1/linhua/Rstudio/INFO_COM/chr01_rufipogon_genotyping.vcf.gz_BIALLELIC_SNPs.vcf.gz_INFO.txt",sep = "\t",header = F,na.strings = ".")
colnames(info_01)<-colnames_info


info_01<-fread("/data1/linhua/Rstudio/INFO_COM/chr01_rufipogon_genotyping.vcf.gz_BIALLELIC_SNPs.vcf.gz_INFO.txt",sep = "\t",header = F,na.strings = ".")
info_02<-fread("/data1/linhua/Rstudio/INFO_COM/chr02_rufipogon_genotyping.vcf.gz_BIALLELIC_SNPs.vcf.gz_INFO.txt",sep = "\t",header = F,na.strings = ".")
info_03<-fread("/data1/linhua/Rstudio/INFO_COM/chr03_rufipogon_genotyping.vcf.gz_BIALLELIC_SNPs.vcf.gz_INFO.txt",sep = "\t",header = F,na.strings = ".")
info_04<-fread("/data1/linhua/Rstudio/INFO_COM/chr04_rufipogon_genotyping.vcf.gz_BIALLELIC_SNPs.vcf.gz_INFO.txt",sep = "\t",header = F,na.strings = ".")
info_05<-fread("/data1/linhua/Rstudio/INFO_COM/chr05_rufipogon_genotyping.vcf.gz_BIALLELIC_SNPs.vcf.gz_INFO.txt",sep = "\t",header = F,na.strings = ".")
info_06<-fread("/data1/linhua/Rstudio/INFO_COM/chr06_rufipogon_genotyping.vcf.gz_BIALLELIC_SNPs.vcf.gz_INFO.txt",sep = "\t",header = F,na.strings = ".")
info_07<-fread("/data1/linhua/Rstudio/INFO_COM/chr07_rufipogon_genotyping.vcf.gz_BIALLELIC_SNPs.vcf.gz_INFO.txt",sep = "\t",header = F,na.strings = ".")
info_08<-fread("/data1/linhua/Rstudio/INFO_COM/chr08_rufipogon_genotyping.vcf.gz_BIALLELIC_SNPs.vcf.gz_INFO.txt",sep = "\t",header = F,na.strings = ".")
info_09<-fread("/data1/linhua/Rstudio/INFO_COM/chr09_rufipogon_genotyping.vcf.gz_BIALLELIC_SNPs.vcf.gz_INFO.txt",sep = "\t",header = F,na.strings = ".")
info_10<-fread("/data1/linhua/Rstudio/INFO_COM/chr10_rufipogon_genotyping.vcf.gz_BIALLELIC_SNPs.vcf.gz_INFO.txt",sep = "\t",header = F,na.strings = ".")
info_11<-fread("/data1/linhua/Rstudio/INFO_COM/chr11_rufipogon_genotyping.vcf.gz_BIALLELIC_SNPs.vcf.gz_INFO.txt",sep = "\t",header = F,na.strings = ".")
info_12<-fread("/data1/linhua/Rstudio/INFO_COM/chr12_rufipogon_genotyping.vcf.gz_BIALLELIC_SNPs.vcf.gz_INFO.txt",sep = "\t",header = F,na.strings = ".")

all_info<-rbind(info_01,
                info_02,
                info_03,
                info_04,
                info_05,
                info_06,
                info_07,
                info_08,
                info_09,
                info_10,
                info_11,
                info_12)

colnames(all_info)<-colnames_info

