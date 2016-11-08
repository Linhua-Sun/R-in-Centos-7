#!/usr/bin/env Rscript

Args<- commandArgs(trailingOnly = TRUE)
input<-Args[1]
terms<-Args[2]

library(data.table)
library(ggpubr)
library(ggsci)
library(crayon,quietly = T)

t1<-fread(input,header = T)

colnames(t1)
print((head(t1)))

p<-ggecdf(t1,x = terms,color ="red",size = 1.2) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

ggsave(filename = paste(basename(input),"_ECDF_curve_",terms,".png",sep = ""),plot = p,device = "png")

cat(yellow("ggdensity is finished\n"))