########### Summarize results ###############

library(dplyr)
library(stringr)
library(mclust)
library(ggplot2)
library(cowplot)

args=(commandArgs(TRUE))
print(args)
BaseFolder<-args[1]
mysample<-args[2]
MinProportion<-args[3]
MinRatio<-args[4]
scRNAannotation<-args[5]

setwd(paste0(BaseFolder, "/results/"))

### Log file
logfile<-file(paste0(mysample, "_results_analysis.log"), open = "w")
writeLines("########## TCR results analysis ########## \n", logfile)
writeLines(paste0(date(), "\n"), logfile)
writeLines(paste0("Starting analysis for sample ", mysample, "\n"), logfile)

### Import TCR alingnment results
db<-read.table("file_igblast_db-pass.tsv", sep = "\t", stringsAsFactors = F, header = T, row.names = NULL)
writeLines(paste0("Imported ", nrow(db), " results \n"), logfile)

## Import consensus building error
x<-read.table("stats.log", fill = T, header = T, sep = "\t", row.names = NULL)
db$error_rate<-x$ERROR[match(db$sequence_id, x$BARCODE)]

## Entries with no CDR3
a<-sum(db$cdr3=="")
writeLines(paste0("Removing ", a, " rows with no CDR3 result \n"), logfile)
#Filter umi that gave no CDR3
db<-db%>%filter(cdr3!="")

### Separate BC and UMI
db$barcode<-substr(db$sequence_id, start = 1, stop = 16)
db$UMI<-substr(db$sequence_id, start = 17, stop = 28)
writeLines(paste0("Found ", length(unique(db$barcode)), " unique cellular barcodes \n"), logfile)

### Arrange
db$consensus_count<-as.integer(db$consensus_count)
db<-db%>%arrange(barcode, desc(consensus_count))

### Filter to selected clusters
writeLines("Filtering by selected TCR sequence clusters ---------------- \n", logfile)
cl<-read.table("BC_UMI_cluster_metrics.txt")
cl.selected<-plotter[plotter$Proportion>MinProportion & plotter$Ratio>=MinRatio,c("BC.UMI", "Cluster")]
cl.selected<-paste0(cl.selected$BC.UMI, "_", cl.selected$Cluster)

db<-db[db$sequence_id%in%cl,]
writeLines(paste0("Selected ", nrow(db), " results passing filter \n"), logfile)
writeLines(paste0("Found ", length(unique(db$barcode)), " unique cellular barcodes \n"), logfile)

### Quantify how many different UMI each BC and CDR3 has
db_UMI_count <- db %>% count(barcode, cdr3) %>% rename(UMI_count = n)
db<-left_join(db, db_UMI_count, by = c('barcode','cdr3'))

### Quantify how many different UMI each cell BC has
db_UMI_count_BC <- db %>% count(barcode) %>% rename(UMI_count_BC = n)
db<-left_join(db, db_UMI_count_BC, by = c('barcode'))

### Arrange
db<- db %>% dplyr::select(barcode,UMI,UMI_count,UMI_count_BC,consensus_count,everything())
db <- db[with(db,order(barcode,UMI_count,consensus_count,decreasing=TRUE)),]

### Error rate
# Summarize error rate per BC, weighted by consensus counts
db$error_rate_weighted<-db$error_rate*db$consensus_count
a<-db %>% group_by(barcode) %>% summarise(error_rate_sum_BC = sum(error_rate_weighted),
                                          counsensus_count_BC = sum(consensus_count)) %>% mutate(weighted_mean_error_rate_BC=error_rate_sum_BC/counsensus_count_BC)
db<-left_join(db, a, by = c('barcode')) %>% dplyr::select(-c(error_rate_sum_BC,counsensus_count_BC,error_rate_weighted))

# Summarize error rate per BC and CDR3, weighted by consensus counts
db$error_rate_weighted<-db$error_rate*db$consensus_count
a<-db %>% group_by(barcode, cdr3) %>% summarise(error_rate_sum_BC_cdr3 = sum(error_rate_weighted),
                                          counsensus_count_BC_cdr3 = sum(consensus_count), .groups = "keep") %>% 
          mutate(weighted_mean_error_rate_BC_cdr3=error_rate_sum_BC_cdr3/counsensus_count_BC_cdr3)
db<-left_join(db, a, by = c('barcode', 'cdr3')) %>% dplyr::select(-c(error_rate_sum_BC_cdr3,counsensus_count_BC_cdr3,error_rate_weighted))

# Save db table
write.csv(db, paste0(mysample, "_barcode_UMI_results.csv"), row.names = F)

### QC
# Distribution of consensus counts per BC + UMI
pa<-ggplot(db, aes(consensus_count)) + geom_histogram(colour = "black", bins = 50) + theme_bw() + ggtitle("Consensus Counts per BC + UMI") + labs(x = "Consensus Counts", y = "BC + UMI Frequency")
writeLines("Distribution of consensus counts per BC + UMI:", logfile)
capture.output(summary(db$consensus_count), file = logfile)
writeLines("\n", logfile)

# Distribution of consensus error rate per BC + UMI
pb<-ggplot(db, aes(error_rate)) + geom_histogram(colour = "black", bins = 50) + theme_bw() + ggtitle("Consensus Error Rate per BC + UMI") + labs(x = "Error Rate", y = "BC + UMI Frequency")
writeLines("Distribution of consensus error rate per BC + UMI:", logfile)
capture.output(summary(db$error_rate), file = logfile)
writeLines("\n", logfile)

pdf("./plots/db_histograms.pdf", width = 12, height = 6)
plot_grid(pa, pb, nrow = 1)
dev.off()

### Separate df into alpha and beta chain into two dataframes
TRA_df <- db %>% filter(locus == "TRA")
TRB_df <- db %>% filter(locus == "TRB")

writeLines(paste0("Found ", nrow(TRA_df), " TRA results and ", nrow(TRB_df), " TRB results \n"), logfile)

### Get nReads per CDR3 
TRA_df<-TRA_df%>%group_by(barcode, cdr3)%>%mutate(TRA_nReads=sum(consensus_count))%>%ungroup()
TRB_df<-TRB_df%>%group_by(barcode, cdr3)%>%mutate(TRB_nReads=sum(consensus_count))%>%ungroup()

### Build new table with one row per BC
writeLines("Building Barcodes table ------------- \n", logfile)
barcodes <- unique(db$barcode)
barcodes = data.frame(row.names = barcodes, BC = barcodes)

# Do we recover TRA, TRB or both?
TRA_BC <- !is.na(match(barcodes$BC, TRA_df$barcode))
TRB_BC <- !is.na(match(barcodes$BC, TRB_df$barcode))
TRA_BC & TRB_BC ->TRAB_BC

barcodes$TCR_Recovery <- "No Recovery"
barcodes$TCR_Recovery[c(TRB_BC)] <- "TRB only"
barcodes$TCR_Recovery[c(TRA_BC)] <- "TRA only"
barcodes$TCR_Recovery[c(TRAB_BC)] <- "TRA and TRB"

writeLines("TCR genes recovery frequencies:", logfile)
capture.output(table(barcodes$TCR_Recovery), file = logfile)
writeLines("\n", logfile)

# -------------------------
# TRB CDR3 sequence and V, D and J segments
barcodes$TRB_CDR3 <- TRB_df$junction_aa[match(barcodes$BC,TRB_df$barcode)]
barcodes$TRB_CDR3nuc <- TRB_df$cdr3[match(barcodes$BC,TRB_df$barcode)]
barcodes$TRB_CDR3_UMIcount <- TRB_df$UMI_count[match(barcodes$BC,TRB_df$barcode)]
barcodes$TRB_CDR3_error <- TRB_df$weighted_mean_error_rate_BC_cdr3[match(barcodes$BC,TRB_df$barcode)]
barcodes$TRB_nReads <- TRB_df$TRB_nReads[match(barcodes$BC,TRB_df$barcode)]

barcodes$TRBV <- TRB_df$v_call[match(barcodes$BC,TRB_df$barcode)]
barcodes$TRBD <- TRB_df$d_call[match(barcodes$BC,TRB_df$barcode)]
barcodes$TRBJ<- TRB_df$j_call[match(barcodes$BC,TRB_df$barcode)]

# TRA CDR3 sequence and V and J segmentes
barcodes$TRA_CDR3 <- TRA_df$junction_aa[match(barcodes$BC,TRA_df$barcode)]
barcodes$TRA_CDR3nuc <- TRA_df$cdr3[match(barcodes$BC,TRA_df$barcode)]
barcodes$TRA_CDR3_UMIcount <- TRA_df$UMI_count[match(barcodes$BC,TRA_df$barcode)]
barcodes$TRA_CDR3_error <- TRA_df$weighted_mean_error_rate_BC_cdr3[match(barcodes$BC,TRA_df$barcode)]
barcodes$TRA_nReads <- TRA_df$TRA_nReads[match(barcodes$BC,TRA_df$barcode)]

barcodes$TRAV <- TRA_df$v_call[match(barcodes$BC,TRA_df$barcode)]
barcodes$TRAJ <- TRA_df$j_call[match(barcodes$BC,TRA_df$barcode)]

## Second TRA allele
# Select 1 row per BC + CDR3 + V GENE
TRA_uniq_df <- TRA_df %>% distinct(barcode,cdr3,v_call, .keep_all = TRUE)

# Remove 1st allele
i<-match(barcodes$BC,TRA_uniq_df$barcode)
i<-i[!is.na(i)]
TRA_df_2<-TRA_uniq_df[-i,]

barcodes$TRA.2_CDR3 <- TRA_df_2$junction_aa[match(barcodes$BC,TRA_df_2$barcode)]
barcodes$TRA.2_CDR3nuc <- TRA_df_2$cdr3[match(barcodes$BC,TRA_df_2$barcode)]
barcodes$TRA.2_CDR3_UMIcount <- TRA_df_2$UMI_count[match(barcodes$BC,TRA_df_2$barcode)]
barcodes$TRA.2_CDR3_error <- TRA_df_2$weighted_mean_error_rate_BC_cdr3[match(barcodes$BC,TRA_df_2$barcode)]
barcodes$TRA.2_nReads <- TRA_df_2$TRA_nReads[match(barcodes$BC,TRA_df_2$barcode)]

barcodes$TRAV.2 <- TRA_df_2$v_call[match(barcodes$BC,TRA_df_2$barcode)]
barcodes$TRAJ.2 <- TRA_df_2$j_call[match(barcodes$BC,TRA_df_2$barcode)]

### Final reads used
nread<-sum(barcodes$TRB_nReads, na.rm = T)+sum(barcodes$TRA_nReads, na.rm = T)+sum(barcodes$TRA.2_nReads, na.rm = T)
writeLines(paste0("Final Read Number: ", nread, "\n"), logfile)

## Save table
write.csv(barcodes, paste0(mysample, "_barcode_results.csv"), row.names = F)

### QC 
# Distribution of CDR3 UMI counts
writeLines("Distribution of TRB CDR3 UMI counts per BC: \n", logfile)
capture.output(summary(barcodes$TRB_CDR3_UMIcount), file = logfile)
writeLines("\n", logfile)

writeLines("Distribution of TRA CDR3 UMI counts per BC: \n", logfile)
capture.output(summary(barcodes$TRA_CDR3_UMIcount), file = logfile)
writeLines("\n", logfile)

### Finish here if there are no RNAseq annotations
if(is.na(scRNAannotation)){
   writeLines("No RNAseq annotations found \n", logfile)
   writeLines("Closing \n", logfile)
   close(logfile)
   stop()
}

### Overlap with scRNAseq
writeLines("Overlapping with scRNAseq ------------ \n", logfile)
bc.sc<-read.table(scRNAannotation, header = T)
bc.sc$V1<-substr(x = bc.sc$V1, start = 1, stop = 16)
head(bc.sc)
writeLines(paste0("Imported ", nrow(bc.sc), " annotated cells \n"), logfile)
writeLines(paste0("Found ", sum(bc.sc$V1%in%barcodes$BC), " overlapping barcodes \n"), logfile)
barcodes$InRNAseq<-FALSE
barcodes$InRNAseq[barcodes$BC%in%bc.sc$V1]<-TRUE
barcodes$RNAannotation<-bc.sc$V3[match(barcodes$BC, bc.sc$V1)]
barcodes$RNAannotation[is.na(barcodes$RNAannotation)]<-"Not_annotated"
nread<-sum(barcodes$TRB_nReads[barcodes$InRNAseq==T], na.rm = T)+sum(barcodes$TRA_nReads[barcodes$InRNAseq==T], na.rm = T)+sum(barcodes$TRA.2_nReads[barcodes$InRNAseq==T], na.rm = T)
writeLines(paste0("Final Read Number overlapping with scRNAseq: ", nread, "\n"), logfile)


### TRA and TRB clones
writeLines("Identifying TRA and TRB clones ------------ \n", logfile)

df<-barcodes %>% filter(InRNAseq==TRUE)
## TRA
dfa<-df %>% filter(TCR_Recovery%in%c("TRA and TRB", "TRA only"))
dfa<-dfa %>% group_by(TRA_CDR3) %>% mutate(TRA_CDR3_CloneSize = n()) %>% ungroup()
# Normalize by total cell number with TRA 
dfa<-dfa %>% mutate(TRA_CDR3_CloneSize_Norm = (TRA_CDR3_CloneSize/nrow(dfa))*100)
# Rank / clone names
dfa <- dfa %>% arrange(desc(TRA_CDR3_CloneSize))
dfa$TRA_CDR3 <- factor(dfa$TRA_CDR3, levels = unique(dfa$TRA_CDR3))
dfa <- dfa %>% group_by(TRA_CDR3) %>% mutate(TRA_CloneID=cur_group_id())
# Add info to barcodes table
x<-dfa[,c("BC","TRA_CDR3", "TRA_CDR3_CloneSize", "TRA_CDR3_CloneSize_Norm", "TRA_CloneID")]
barcodes<-left_join(barcodes, x, by = c("BC", "TRA_CDR3"))
writeLines(paste0("Found ", max(barcodes$TRA_CloneID, na.rm = T), " TRA CDR3 distinct clonse \n"), logfile)

## TRB
dfb<-df %>% filter(TCR_Recovery%in%c("TRA and TRB", "TRB only"))
dfb<-dfb %>% group_by(TRB_CDR3) %>% mutate(TRB_CDR3_CloneSize = n()) %>% ungroup()
# Normalize by total cell number with TRB
dfb<-dfb %>% mutate(TRB_CDR3_CloneSize_Norm = (TRB_CDR3_CloneSize/nrow(dfb))*100)
# Rank / clone names
dfb <- dfb %>% arrange(desc(TRB_CDR3_CloneSize))
dfb$TRB_CDR3 <- factor(dfb$TRB_CDR3, levels = unique(dfb$TRB_CDR3))
dfb <- dfb %>% group_by(TRB_CDR3) %>% mutate(TRB_CloneID=cur_group_id())

# Add info to barcodes table
x<-dfb[,c("BC","TRB_CDR3", "TRB_CDR3_CloneSize", "TRB_CDR3_CloneSize_Norm", "TRB_CloneID")]
barcodes<-left_join(barcodes, x, by = c("BC", "TRB_CDR3"))
writeLines(paste0("Found ", max(barcodes$TRB_CloneID, na.rm = T), " TRB CDR3 distinct clonse \n"), logfile)

## ARI TRA and TRB clones
ari<-adjustedRandIndex(barcodes$TRA_CloneID, barcodes$TRB_CloneID)
writeLines(paste0("Adjusted Rand Index for TRA and TRB CDR3 clones: ", ari, "\n"), logfile)

## Save table
write.csv(barcodes, paste0(mysample, "_barcode_results.csv"), row.names = F)

close(logfile)
