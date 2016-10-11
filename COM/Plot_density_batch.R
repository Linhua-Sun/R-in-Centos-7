#!/usr/bin/env Rscript

Args<- commandArgs(trailingOnly = TRUE)

input<-Args[1]
input="/data1/linhua/Rstudio/INFO_COM/10_4.txt"
library(data.table)

library(ggpubr)

colnames_info<-c("CHROM","POS","REF,ALT","AC","AF","AN","BaseQRankSum","ClippingRankSum","DP","DS","END","ExcessHet","FS","HaplotypeScore","InbreedingCoeff","MLEAC","MLEAF","MQ","MQRankSum","QD","RAW_MQ","ReadPosRankSum","SOR")

info_01<-fread(input,sep = "\t",header = F,na.strings = ".")

colnames(info_01)<-colnames_info

pathname<-paste("./",basename(input),"_DENSITY_PLOT_OUTPUT",sep = "")

cat(pathname,"\n")

if (file.exists(pathname)==T) {
    cat(green("dir is exit! makedir another new one"))
    new_pathname<-paste(pathname,"_new",sep="")
    dir.create(new_pathname)
} else {
    dir.create(path = pathname) 
    cat(green("Create a new dir!"))
}






for (i in colnames(info_01)[4:length(colnames(info_01))]) {
    cat(i,fill = T)
    pdffilepath<-paste(pathname,"/",basename(input), "_density_plot_",i,".pdf" ,sep="")
    cat(pdffilepath,fill=T)
    logout<-paste(i," pdf process is OK !",sep = "")
    tryCatch({
        p<-ggdensity(info_01, x = i,
                     add = "mean",
                     palette ="npg" ,rug = T)
        ggsave(filename = pdffilepath,plot = p,device = "pdf")
    },
    error=function(e){cat("There are all NAs. ",conditionMessage(e),fill = T)},
    finally={cat(logout,fill=T)} 
    )
}

##Refer to R coding from http://www.cnblogs.com/weibaar/p/4382397.html;trycatch.

if (file.exists(TTT)==T) cat("1")
