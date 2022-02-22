#!/usr/bin/Rscript

library(ggplot2)
library(reshape2)
library(cowplot)
library(pheatmap)
library(ComplexHeatmap)
library(pals)
library(gridExtra)
library(RColorBrewer)
library(e1071)
library(dplyr)
library(stringr)
#source("/broad/vangalenlab/mainciburu/TCRseq_scripts/colors.r")

# Input --------------------
args=(commandArgs(TRUE))

BaseFolder<-args[1]
MySample<-args[2]
scRNAannotation<-args[3]
BClength<-as.numeric(args[4])
UMIlength<-as.numeric(args[5])

print(paste0("BaseFolder: ", BaseFolder))
print(paste0("MySample: ", MySample))
print(paste0("scRNAannotation: ", scRNAannotation))
print(paste0("BClength: ", BClength))
print(paste0("UMIlength: ", UMIlength))

setwd(paste0(BaseFolder, "/downstream"))

barcodes<-read.csv(paste0(MySample, "_barcode_results.csv"), stringsAsFactors=FALSE)

bc.sc<-read.table(scRNAannotation, header = F, stringsAsFactors = F)
bc.sc$V1<-substr(x = bc.sc$V1, start = 1, stop = BClength)

# Colors -----------------------
col.tcr.recovery<-c("No Results"="#B8B8B8FF","No Recovery"="#357EBDFF","TRA and TRB"="#EEA236FF",
          "TRA only"="#9632B8FF","TRB only"="#5CB85CFF")

col.celltype<-c(kelly(), cols25())
col.celltype<-sample(col.celltype, size = length(unique(bc.sc$V2)), replace = F)
names(col.celltype)<-unique(bc.sc$V2)

# TCR recovery ----------------
# How many scRNAseq cells have TCRseq results?
bc.sc$TCR<-barcodes$TCR_Recovery[match(bc.sc$V1, barcodes$BC)]
bc.sc$TCR[is.na(bc.sc$TCR)]<-"No Results"
tt<-prop.table(table(bc.sc$V2, bc.sc$TCR), margin = 1)
df<-melt(tt)

pdf("./plots/scRNAseq_TCRrecovery_proportions.pdf", width = 7, height = 5)
ggplot(df, aes(Var1, value, fill = Var2)) + geom_bar(stat = "identity", position = position_stack()) + scale_fill_manual(values = col.tcr.recovery) +
      labs(x = "Cell Type", y = "Proportion", fill = "TCR Recovery") + theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

# UMI count distribution ----------
# UMI counts contributing to the selected CDR3 sequence in each barcode
dfa<- barcodes %>% filter(TCR_Recovery%in%c("TRA and TRB", "TRA only"))
dfa$TRA_CDR3_UMIcount<-factor(dfa$TRA_CDR3_UMIcount, levels = 1:max(dfa$TRA_CDR3_UMIcount))
p1<-ggplot(dfa, aes(TRA_CDR3_UMIcount, fill = InRNAseq)) + geom_bar(colour = "black", position = position_dodge()) + 
  scale_fill_manual(values=c("#69b3a2", "#404080"), labels = c("No", "Yes")) + 
  theme_bw() + labs(x = "TRA CDR3 UMI count", y = "Barcode Number", fill = "In scRNAseq") + ggtitle("TRA CDR3 UMI count")

dfb<- barcodes %>% filter(TCR_Recovery%in%c("TRA and TRB", "TRB only"))
dfb$TRB_CDR3_UMIcount<-factor(dfb$TRB_CDR3_UMIcount, levels = 1:max(dfb$TRB_CDR3_UMIcount))
p2<-ggplot(dfb, aes(TRB_CDR3_UMIcount, fill = InRNAseq)) + geom_bar(colour = "black", position = position_dodge()) + 
  scale_fill_manual(values=c("#69b3a2", "#404080"), labels = c("No", "Yes")) + 
  theme_bw() + labs(x = "TRB CDR3 UMI count", y = "Barcode Number", fill = "In scRNAseq") + ggtitle("TRB CDR3 UMI count")

pdf("./plots/CDR3_UMIcount_distribution.pdf", width = 12, height = 6)
plot_grid(p1, p2, nrow = 1)
dev.off()

# Valid reads ---------
df<-barcodes[,c("BC", "TCR_Recovery", "TRB_nReads", 
                "TRA_nReads", "TRA.2_nReads", "RNAannotation")]
df<-df%>%group_by(RNAannotation)%>%summarise(TRA=sum(TRA_nReads, na.rm = T)+sum(TRA.2_nReads, na.rm = T), 
                                             TRB=sum(TRB_nReads, na.rm = T))
df<-melt(df)
colnames(df)<-c("CellAnnotation", "Gene", "nReads")
#df$CellAnnotation <- factor(df$CellAnnotation, levels = c(celltypes.ch, "Not_annotated"))
pdf("./plots/valid_reads.pdf", width = 8, height = 5)
ggplot(df, aes(CellAnnotation, nReads, fill = Gene)) + geom_bar(stat = "identity", position = position_dodge()) +
  theme_bw() + scale_fill_manual(values = c(TRA="#9632B8FF", TRB="#5CB85CFF")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + labs(y = "nReads") + ggtitle(paste0(sum(df$nReads), " Reads in final TCR analysis"))
dev.off()

# TRA and TRB correspondence -------------
# Match TRA CDR3 and TRB CRD3 in RNAseq cells
barcodes.ov<-barcodes[barcodes$InRNAseq==TRUE,]
#sort(table(barcodes.ov$TRB_CDR3), decreasing = T)[1:10]
#sort(table(barcodes.ov$TRA_CDR3), decreasing = T)[1:10]
barcodes.ov.trab<-barcodes.ov[barcodes.ov$TCR_Recovery=="TRA and TRB",]

## Filter clones by cell number
nx<-0
barcodes.ov.trab.filter<-barcodes.ov.trab
barcodes.ov.trab.filter<-barcodes.ov.trab.filter %>% group_by(TRA_CDR3) %>% 
  mutate(freq = n()) %>% ungroup() %>% filter(freq > nx) %>% dplyr::select(-freq)
barcodes.ov.trab.filter<-barcodes.ov.trab.filter %>% group_by(TRB_CDR3) %>% 
  mutate(freq = n()) %>% ungroup() %>% filter(freq > nx) %>% dplyr::select(-freq)
dim(barcodes.ov.trab.filter)

tt<-table(barcodes.ov.trab.filter$TRA_CDR3, 
          barcodes.ov.trab.filter$TRB_CDR3)
col<-brewer.pal(n = 9, "YlOrRd")
col<-c("white", colorRampPalette(col)(max(tt)))
pheatmap::pheatmap(tt, color = col,
                   width = 8, height = 9,
                   fontsize = 10,
                   filename = "./plots/CDR3_clones_heatmap.pdf")

heatmap.tib <- barcodes.ov.trab.filter %>% group_by(TRB_CDR3) %>% mutate(trb_clone_size = n()) %>% #filter(clone_size > 10) %>%
  group_by(TRA_CDR3) %>% mutate(tra_clone_size = n()) %>%
  arrange(desc(trb_clone_size), desc(tra_clone_size), TRB_CDR3, TRA_CDR3) %>% 
  dplyr::select(TRB_CDR3, TRA_CDR3, RNAannotation) %>% na.omit()
heatmap.tib <- heatmap.tib %>%
  mutate(TRB_CDR3 = factor(TRB_CDR3, levels = unique(TRB_CDR3))) %>%
  mutate(TRA_CDR3 = factor(TRA_CDR3, levels = unique(TRA_CDR3))) 

# Plot
cols<-c(kelly(), cols25())
myreplace<-ifelse(length(unique(barcodes.ov.trab.filter$TRA_CDR3))>length(cols), T, F)
cols<-c(sample(x = cols, size = length(unique(barcodes.ov.trab.filter$TRB_CDR3)), replace = myreplace),
        sample(x = cols, size = length(unique(barcodes.ov.trab.filter$TRA_CDR3)), replace = myreplace),
        col.celltype)

names(cols)<-c(unique(barcodes.ov.trab.filter$TRB_CDR3), 
               unique(barcodes.ov.trab.filter$TRA_CDR3),
               names(col.celltype))
ht <- Heatmap(t(as.matrix(heatmap.tib)),
              col = cols,
              na_col = "#FFFFFF",
              cluster_rows = F,
              cluster_columns = F,
              show_row_names = T,
              show_column_names = T,
              show_heatmap_legend = F,
              column_title = str_c(nrow(heatmap.tib), " cells, clones >", nx, " cells"))

# Assemble three legends
# *** pending *** 
# Make legend work!!
#TRB.lgd <- Legend(at = as.character(unique(heatmap.tib$TRB_CDR3)),
#                  legend_gp = gpar(fill = cols[as.character(unique(heatmap.tib$TRB_CDR3))]), 
#                  title = "TRB variable region")
#TRA.lgd <- Legend(at = as.character(unique(heatmap.tib$TRA_CDR3)),
#                  legend_gp = gpar(fill = cols[as.character(unique(heatmap.tib$TRA_CDR3))]), title = "TRA variable region")
#CellType.lgd <- Legend(at = as.character(unique(heatmap.tib$RNAannotation)),
#                  legend_gp = gpar(fill = cols[as.character(unique(heatmap.tib$RNAannotation))]), title = "Cell type")
#legend_pack <- packLegend(TRB.lgd, TRA.lgd, CellType.lgd, direction = "horizontal" )

# Plot
grobs.ls <- list(grid.grabExpr(draw(ht)))#, grid.grabExpr(draw(legend_pack)))

pdf("./plots/TRB_TRA_correspondence.pdf", width = 8, height = 12)
grid.arrange(grobs = grobs.ls, ncol = 1, nrow = 2, heights = c(1,3))
dev.off()


# Clone size ranking -------------------
# TRA
dfa.clones <- barcodes %>% filter(InRNAseq==TRUE, TCR_Recovery%in%c("TRA and TRB", "TRA only")) %>% 
              group_by(TRA_CDR3) %>% distinct(TRA_CDR3, .keep_all = T)
# TRB 
dfb.clones <- barcodes %>% filter(InRNAseq==TRUE, TCR_Recovery%in%c("TRA and TRB", "TRB only")) %>% 
              group_by(TRB_CDR3) %>% distinct(TRB_CDR3, .keep_all = T)

pa1<-ggplot(dfa.clones, aes(TRA_CloneID, TRA_CDR3_CloneSize)) + geom_point(size = 0.5) + theme_bw() +
    labs(x = "TRA Clone Rank", y = "TRA Clone Size (cell number)") + ggtitle("TRA Clone Size")
pa2<-ggplot(dfa.clones, aes(TRA_CloneID, TRA_CDR3_CloneSize_Norm)) + geom_point(size = 0.5) + theme_bw() +
    labs(x = "TRA Clone Rank", y = "Normalized TRA Clone Size (% of cells)") + ggtitle("TRA Normalized Clone Size")
pb1<-ggplot(dfb.clones, aes(TRB_CloneID, TRB_CDR3_CloneSize)) + geom_point(size = 0.5) + theme_bw() +
    labs(x = "TRB Clone Rank", y = "TRB Clone Size (cell number)") + ggtitle("TRB Clone Size")
pb2<-ggplot(dfb.clones, aes(TRB_CloneID, TRB_CDR3_CloneSize_Norm)) + geom_point(size = 0.5) + theme_bw() + 
    labs(x = "TRB Clone Rank", y = "Normalized TRB Clone Size (% of cells)") + ggtitle("TRB Normalized Clone Size")

pdf("./plots/TRA_TRB_clone_size.pdf", width = 10, height = 8)
plot_grid(pa1, pb1, pa2, pb2, nrow = 2)
dev.off()

# Clone size vs cell type -------------
dfb <- barcodes %>% filter(InRNAseq==TRUE, TCR_Recovery%in%c("TRA and TRB", "TRB only"))
#dfb$RNAannotation<-factor(dfb$RNAannotation, levels = celltypes.ch)
pdf("./plots/TRB_clone_size_celltype.pdf", width = 8, height = 6)
ggplot(dfb, aes(RNAannotation, TRB_CDR3_CloneSize_Norm, color = RNAannotation)) + 
  geom_jitter(size = 1, width = 0.3, height = 0) + theme_bw() + 
  scale_color_manual(values = col.celltype) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

# Top 50 TRB clones ---------------------
dfb.top<-barcodes[barcodes$TRB_CloneID<=50&!is.na(barcodes$TRB_CloneID),]
#dfb.top$RNAannotation <- factor(dfb.top$RNAannotation, levels = celltypes.ch)

pp<-ggplot(dfb.top, aes(TRB_CloneID, fill = RNAannotation)) + geom_bar(position = "stack", colour = "black", size = 0.1) + 
    theme_bw()  + scale_fill_manual(values = col.celltype) +
    labs(x = "TRB Clone Rank", y = "TRB Clone Size (cell number)")
pdf("./plots/trb_top_clones.pdf", width = 10, height = 6)
print(pp)
dev.off()

# Normalized
plotter <- dfb.top %>% group_by(TRB_CDR3, RNAannotation) %>% mutate(CellType_Count = n()) %>% 
  mutate(Clone_Composition_Norm = (CellType_Count/TRB_CDR3_CloneSize)*TRB_CDR3_CloneSize_Norm) %>%
  distinct(TRB_CDR3, RNAannotation, .keep_all = T)
pp<-ggplot(plotter, aes(TRB_CloneID, y = Clone_Composition_Norm, fill = RNAannotation)) + 
    geom_bar(position = "stack", stat = "identity", colour = "black", size = 0.1) +
    theme_bw() + scale_fill_manual(values = col.celltype) +
    labs(x = "TRB Clone Rank", y = "Normalized TRB Clone Size (% of cells)") + ggtitle(MySample)
pdf("./plots/trb_top_clones_norm.pdf", width = 10, height = 6)
print(pp)
dev.off()


# CDR3 Hamming distances (Top 50 TRB) ----------
TRB_CDR3<-plotter$TRB_CDR3
TRB_CDR3<-unique(TRB_CDR3[!is.na(TRB_CDR3)])
TRB_CDR3<-strsplit(TRB_CDR3, split = "")
n<-max(sapply(TRB_CDR3, length))
TRB_CDR3<-lapply(TRB_CDR3, function(x){
  if(length(x)<n){x<-c(x, rep("", n-length(x)))}
  else{x}
})
x<-do.call("rbind", TRB_CDR3)
hd<-hamming.distance(x)
rownames(hd)<-colnames(hd)<-unique(plotter$TRB_CDR3)
ix<-which(!duplicated(plotter$TRB_CloneID))
annot<-data.frame(Ranking=plotter$TRB_CloneID[ix], 
                  NormSize=plotter$TRB_CDR3_CloneSize_Norm[ix],
                  row.names = plotter$TRB_CDR3[ix])
annot.col<-list(Ranking=colorRampPalette(rev(brewer.pal(11, "RdYlGn")))(50),
                NormSize = colorRampPalette(colors = c("lightgrey", "navy"))(200))
annot<-annot[order(annot$Ranking, decreasing = F),,drop=FALSE]
hd<-hd[rownames(annot),rownames(annot)]
pheatmap::pheatmap(hd, cluster_rows = F, cluster_cols = F,
         show_rownames = T, show_colnames = T, 
         annotation_row = annot, annotation_colors = annot.col, 
         filename = "./plots/TRB_distance_heatmap.pdf", width = 12, height = 12, main = MySample)


