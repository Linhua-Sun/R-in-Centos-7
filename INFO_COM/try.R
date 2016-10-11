#!/usr/bin/env Rscript

Args<- commandArgs()
input<-Args[6]

library(data.table)
library(ggpubr)

colnames_info<-c("CHROM","POS","REF,ALT","AC","AF","AN","BaseQRankSum","ClippingRankSum","DP","DS","END","ExcessHet","FS","HaplotypeScore","InbreedingCoeff","MLEAC","MLEAF","MQ","MQRankSum","QD","RAW_MQ","ReadPosRankSum","SOR")

info_01<-fread(input,sep = "\t",header = F,na.strings = ".")

colnames(info_01)<-colnames_info

pathname<-paste("./",basename(input),"_DENSITY_PLOT_OUTPUT",sep = "")

print(pathname)

dir.create(path = pathname)

print("Create a new dir!")

for (i in colnames(info_01)[4:length(colnames(info_01))]) {
    print(i)
    pngfilepath<-paste(pathname,"/",basename(input), "_density_plot_",i,".png" ,sep="")
    print(pngfilepath)
    png(filename = pngfilepath)
    ggdensity(info_01, x = i,
              add = "mean",
              palette ="npg" ,rug = T)
    dev.off()
    print("png is OK !")
}
getwd()
