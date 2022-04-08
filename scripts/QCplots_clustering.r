#!/usr/bin/Rscript

############################
## QC Plots - clustering
############################


# args[1] = BaseFolder

args=(commandArgs(TRUE))
print(args)
BaseFolder<-args[1]
BClength<-as.integer(args[2])
UMIlength<-as.integer(args[3])

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(data.table)

setwd(BaseFolder)

reads<-fread("./fastq_processed/BCSeq_final_filtered_qfiltered_cluster.txt", stringsAsFactors=F, header=F)
x<-data.table(BC.UMI=substr(x = reads$V1, start = 1, stop = BClength+UMIlength),
              Cluster=substr(x = reads$V1, start = BClength+UMIlength+2, stop = nchar(reads$V1)))
x$Cluster<-factor(x$Cluster, levels = 1:max(as.integer(x$Cluster)))

x<-x %>% count(BC.UMI, Cluster) %>% rename(Frequency=n) %>% arrange(BC.UMI, desc(Frequency)) 
x<-x %>% group_by(BC.UMI) %>% mutate(Frequency.BC.UMI=sum(Frequency))
x$Proportion<-x$Frequency/x$Frequency.BC.UMI

freq.first<-x %>% group_by(BC.UMI) %>% summarise(freq.first = max(Frequency))
x<-left_join(x, freq.first, by = "BC.UMI")

ix<-names(table(x$BC.UMI))[table(x$BC.UMI)>1]
freq.second<-x[x$BC.UMI%in%ix,]
freq.second<-freq.second %>% group_by(BC.UMI) %>%  filter(rank(-Frequency, ties.method="first")==2) %>% select(BC.UMI, Frequency) %>% rename(freq.second=Frequency)
x<-left_join(x, freq.second, by = "BC.UMI", keep = FALSE)
x$Ratio<-x$freq.first/x$freq.second

i<-!duplicated(x$BC.UMI)
plotter<-data.table(BC.UMI=x$BC.UMI[i],
                    Cluster=x$Cluster[i],
                    Frequency.BC.UMI=x$Frequency.BC.UMI[i],
                    Frequency.BC.UMI.Cluster=x$Frequency[i],
                    Proportion=x$Proportion[i],
                    Ratio=x$Ratio[i])

plotter$Ratio[is.na(plotter$Ratio)]<-1
plotter$logRatio<-log2(plotter$Ratio)
plotter$logFrequency.BC.UMI<-log2(plotter$Frequency.BC.UMI)
nread<-sum(plotter$Frequency.BC.UMI.Cluster[plotter$Proportion>0.5&plotter$Ratio>2])
pct.read<-round((nread/sum(plotter$Frequency.BC.UMI))*100, 2)

pdf("./wat3r/QC/QCplot_clusters.pdf", width = 15, height = 10)
ggplot(plotter, aes(Proportion, logRatio, color=logFrequency.BC.UMI)) + geom_point(size = 0.05) + 
  theme_bw() + scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral"))) + 
  geom_vline(xintercept = 0.5, colour="red") + geom_hline(yintercept = 1, colour="red") + ylim(c(0, round(max(plotter$logRatio))+1)) + 
  theme(text=element_text(size=20)) + ggtitle("TCR Sequences Clustering", subtitle=paste0(nread, " (", pct.read, "%)", " reads passing default filters (Proportion > 0.5, Ratio > 2)"))
dev.off()

fwrite(plotter, file = "./wat3r/QC/BC_UMI_cluster_metrics.txt", quote = F)
