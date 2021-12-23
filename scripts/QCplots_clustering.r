#!/usr/bin/Rscript

############################
## QC Plots - clustering
############################


# args[1] = BaseFolder

args=(commandArgs(TRUE))
print(args)
BaseFolder<-args[1]

library(ggplot2)
library(dplyr)
library(RColorBrewer)

setwd(BaseFolder)

reads<-read.table("./fastq_processed/BCSeq_final_filtered_qfiltered_cluster.txt", stringsAsFactors=F)
x<-data.frame(BC.UMI=substr(x = reads$V1, start = 1, stop = 28),
              Cluster=substr(x = reads$V1, start = 30, stop = nchar(reads$V1)))
x$Cluster<-factor(x$Cluster, levels = 1:max(as.integer(x$Cluster)))

x<-x %>% count(BC.UMI, Cluster) %>% rename(Frequency=n) %>% arrange(BC.UMI, desc(Frequency)) 
x<-x %>% group_by(BC.UMI) %>% mutate(Frequency.BC.UMI=sum(Frequency))
x$Proportion<-x$Frequency/x$Frequency.BC.UMI
ratio<-sapply(unique(x$BC.UMI), function(ix){
  dat<-x[x$BC.UMI==ix,]
  if(nrow(dat)==1){return(0)}
  else{return(dat$Frequency[1]/dat$Frequency[2])}
})

i<-!duplicated(x$BC.UMI)
plotter<-data.frame(BC.UMI=x$BC.UMI[i],
                    Cluster=x$Cluster[i],
                    Frequency.BC.UMI=x$Frequency.BC.UMI[i],
                    Frequency.BC.UMI.Cluster=x$Frequency[i],
                    Proportion=x$Proportion[i],
                    Ratio=ratio)
plotter$logRatio<-log2(plotter$Ratio+0.05)
plotter$logFrequency.BC.UMI<-log2(plotter$Frequency.BC.UMI)
nread<-sum(plotter$Frequency.BC.UMI.Cluster[plotter$Proportion>0.5&plotter$Ratio>2])
pct.read<-round((nread/sum(plotter$Frequency.BC.UMI))*100, 2)

pdf("./wat3r/QC/QCplot_clusters.pdf", width = 15, height = 10)
ggplot(plotter, aes(Proportion, logRatio, color=logFrequency.BC.UMI)) + geom_point(size = 0.05) + 
  theme_bw() + scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral"))) + 
  geom_vline(xintercept = 0.5, colour="red") + geom_hline(yintercept = 1, colour="red") + ylim(c(0, round(max(plotter$logRatio))+1)) + 
  theme(text=element_text(size=20)) + ggtitle("TCR Sequences Clustering", subtitle=paste0(nread, " (", pct.read, "%)", " reads passing default filters (Proportion > 0.5, Ratio > 2)"))
dev.off()

write.table(plotter, file = "./wat3r/QC/BC_UMI_cluster_metrics.txt", quote = F)
