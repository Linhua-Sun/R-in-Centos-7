#!/usr/bin/env Rscript

## important project in analysis the data from bcftools
## Linhua Sun
## Mon Oct 10 15:06:50 CST 2016
################################################################################################################
#http://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
#http://gatkforums.broadinstitute.org/gatk/discussion/6925/understanding-and-adapting-the-generic-hard-filtering-recommendations
## --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || 00M0QRankSum < -12.5 || ReadPosRankSum < -8.0" \ ##
################################################################################################################


Args<- commandArgs()

input<-Args[6]

library(data.table)

library(ggpubr)

colnames_info<-c("CHROM","POS","REF,ALT","AC","AF","AN","BaseQRankSum","ClippingRankSum","DP","DS","END","ExcessHet","FS","HaplotypeScore","InbreedingCoeff","MLEAC","MLEAF","MQ","MQRankSum","QD","RAW_MQ","ReadPosRankSum","SOR")

info_01<-fread(input,sep = "\t",header = F,na.strings = ".")

colnames(info_01)<-colnames_info

pathname<-paste("./",input,"_OUTPUT",sep = "")

print(pathname)

dir.create(path = pathname)

print("Create a new dir!")

#1
filepath<-paste(pathname,"/",input, "_density_plot.pdf" ,sep="")

pdf(file = filepath)

ggdensity(info_01, x = "QD",
          add = "mean",
          palette ="npg" ,rug = T)
dev.off()

print("PDF is OK !")
#2
pngfilepath<-paste(pathname,"/",input, "_density_plot.png" ,sep="")
png(filename = pngfilepath)
ggdensity(info_01, x = "QD",
          add = "mean",
          palette ="npg" ,rug = T)
dev.off()

print("png is OK !")
