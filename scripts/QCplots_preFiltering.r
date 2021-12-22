#!/usr/bin/Rscript

################################
## QC Plots - qscore filtering
################################

args=(commandArgs(TRUE))
print(args)
BaseFolder<-args[1]

library(ggplot2)
library(dplyr)

setwd(paste0(BaseFolder, "/wat3r"))

ns<-try(read.table("./QC/ns.txt"))
if(class(ns)=="try-error"){
  ns<-NULL
}

qscore<-read.table("./QC/qscore.txt")

if(!is.null(ns)){
  pdf("./QC/QCplots_preFiltering.pdf", width = 15, height = 10, useDingbats = F)

  # Reads with >0 Ns => N number vs average Q score
  plotter<-data.frame(Ns=ns$V1, Qscore=qscore$V2[ns$V2])
  n<-nrow(plotter)
  p1<-ggplot(plotter, aes(Ns, Qscore)) + geom_point(size = 0.1) + theme_bw() +
    geom_hline(yintercept = 25, color = "red") + ggtitle(paste0(n, " Reads with Ns")) + theme(text=element_text(size=16))

  # Reads with no Ns => average Q score distribution
  plotter<-data.frame(Qscore=qscore$V2[-(ns$V2)])
  xmax<-round(max(plotter$Qscore))+1
  xmin<-round(min(plotter$Qscore))-1
  n<-nrow(plotter)
  p2<-ggplot(plotter, aes(Qscore)) + geom_histogram(binwidth = 1, color="black") + 
    scale_x_reverse(breaks=xmax:xmin) + theme_bw() + geom_vline(xintercept = 25, color = "red") + 
    ggtitle(paste0(n, " Reads with no Ns")) + theme(text=element_text(size=16))
  print(p1)
  print(p2)
  dev.off()
}

if(is.null(ns)){
  pdf("./QC/QCplots_preFiltering.pdf", width = 15, height = 10, useDingbats = F)
  # Reads with no Ns => average Q score distribution
  plotter<-data.frame(Qscore=qscore$V2)
  xmax<-round(max(plotter$Qscore))+1
  xmin<-round(min(plotter$Qscore))-1
  n<-nrow(plotter)
  p1<-ggplot(plotter, aes(Qscore)) + geom_histogram(binwidth = 1, color="black") + 
    scale_x_reverse(breaks=xmax:xmin) + theme_bw() + geom_vline(xintercept = 25, color = "red") + 
    ggtitle(paste0(n, " Reads with no Ns")) + theme(text=element_text(size=16))
  print(p1)
  dev.off()
}
