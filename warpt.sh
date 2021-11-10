#############################################
### WARPT pipeline
#############################################


#!/bin/bash

###################################### Help ################################################
Help()
{
   # Display Help
   echo
   echo "Workflow for Association of Receptor Pairs from TREK-seq"
   echo
   echo "Help:"
   echo "Syntax: warpt.sh [--bc --tcr --sample --rna --wd]"
   echo "options:"
   echo "   -h     Print this Help."
   echo "   -b     cell barcode and UMI fastq.gz full path and name"
   echo "   -t     tcr sequence fastq.gz full path and name"
   echo "   -s     sample name"
   echo "   -r     RNA seq annotations full path and name"
   echo "   -d     working directory for analysis"
   echo "   -c     perform barcode correction"
   echo "   -u     perform UMI correction"
   echo
}

################################ Default variables ########################################
BCfastq="bc fastq"
TCRfastq="tcr fastq"
BaseFolder=$(pwd)
MySample="mysample"
scRNAannotation="scrnaanotation"
CorrectBC="false"
CorrectUMI="false"
MinProportion=0.5
MinRatio=2

####################################### Options ###########################################

while getopts ":hb:t:s:r:d:cu" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      b) 
         BCfastq=$OPTARG;;
      t) 
         TCRfastq=$OPTARG;;
      s) 
         MySample=$OPTARG;;
      r) 
         scRNAannotation=$OPTARG;;
      d) 
         BaseFolder=$OPTARG;;
      c) 
         CorrectBC="true";;
      u) 
         CorrectUMI="true";;
      
      :)    # If expected argument omitted:
         echo "Error: -${OPTARG} requires an argument."
         exit;;
      \?) # Invalid option
         echo "Error: Invalid option $OPTARG. Try warpt.sh -h"
         exit;;
   esac
done


## Load software
source /broad/software/scripts/useuse 
reuse R-4.0
export PATH=$PATH:/broad/vangalenlab/ainciburu/TCRseq_scripts/


# create folders
if [ ! -e ${BaseFolder}/fastq_processed/ ]; then
	mkdir ${BaseFolder}/fastq_processed/
fi

if [ ! -e ${BaseFolder}/warpt/ ]; then
	mkdir ${BaseFolder}/warpt/
	mkdir ${BaseFolder}/warpt/QC/
fi

if [ ! -e ${BaseFolder}/results/ ]; then
	mkdir ${BaseFolder}/results/
	mkdir ${BaseFolder}/results/plots/
fi



# Log file 
touch ${BaseFolder}/warpt/pipeline.log
echo -e "REFORMATING FASTQ FILES \n" |& tee -a ${BaseFolder}/warpt/pipeline.log

# Rearrange Fastq --------------------------------
cd ${BaseFolder}/fastq_processed/

## Get Cell barcode + UMI
zcat ${BCfastq} | awk 'NR%4==2' > "BCSeq.txt"

## Cell Barcode + UMI quality
zcat ${BCfastq} | awk 'NR%4==0' > "BCQC.txt"

## Get TCR sequence
zcat ${TCRfastq} | awk 'NR%4==2' > "TCRSeq.txt" 

## Get TCR quality 
zcat ${TCRfastq} | awk 'NR%4==0' > "TCRQC.txt"  

## Get header 
zcat ${BCfastq} | awk 'NR%4==1' | awk 'BEGIN{FS=" "};{print $1};END{}' > "seqHeaders.txt" 

## Add BC - UMI to header
paste -d ":" "seqHeaders.txt" "BCSeq.txt"  > "seqHeaders_BC.txt"

## Add qualHeader (line 3)
sed 's~@~+~' "seqHeaders.txt" > "qualHeaders.txt"

## Correct barcodes
echo -e "CORRECTING BARCODES \n" |& tee -a ${BaseFolder}/warpt/pipeline.log
python /broad/vangalenlab/ainciburu/TCRseq_scripts/CorrectBC.py -b BCSeq.txt  \
						  -w /broad/vangalenlab/ainciburu/references/3M-february-2018.txt.gz \
						  -d 1 \
						  -r 8 \
						  -o BCcorrected |& tee -a ${BaseFolder}/warpt/pipeline.log

## extract UMIs and join to corrected barcodes
cat BCSeq.txt | awk '{print substr($0,17)}' > UMI.txt
paste -d '' BCcorrected.txt UMI.txt > BCcorrectedUMI.txt

## Mask every UMI with no valid barcode
cat BCcorrectedUMI.txt | sed '/NNNNNNNNNNNNNNNN/c\NNNNNNNNNNNNNNNNNNNNNNNNNNNN' > BCcorrectedUMImasked.txt 

## Correct UMIs
echo -e "CORRECTING UMI \n" |& tee -a ${BaseFolder}/warpt/pipeline.log
python /broad/vangalenlab/ainciburu/TCRseq_scripts/CorrectUMI.py -b BCcorrectedUMImasked.txt \
							-o UMIcorrected |& tee -a ${BaseFolder}/warpt/pipeline.log

## Join corrected BC and UMI
paste -d '' BCcorrected.txt UMIcorrected.txt > BCSeq_corrected.txt

## New fastq
paste -d '' BCSeq_corrected.txt TCRSeq.txt > Read1.txt
paste -d '' BCQC.txt TCRQC.txt > Read1_Q.txt
paste -d '\n' seqHeaders.txt Read1.txt qualHeaders.txt Read1_Q.txt > sample_corrected.fastq

## Remove reads with NNNN...N in BC + UMI

cat sample_corrected.fastq | paste - - - - | awk -F '\t' '{if ($2 !~/^NNNNNNNNNNNNNNNNNNNNNNNNNNNN/){ print $0}}'| tr "\t" "\n" > sample_corrected_filtered.fastq

mv sample_corrected_filtered.fastq ${BaseFolder}/warpt/

### pRESTO / Change-O pipeline
cd ${BaseFolder}/warpt/

# 1) Find barcodes and convert them in tag
# MaskPrimers extract => look for the barcodes in a fixed sequence region
	# -s: input fastq
	# --start: starting position of the sequence to extract
	# --len: lenght of the sequence to extract
	# --pf: name for the resulting tag containing the barcode
	# --mode cut: remove barcode region from sequence
	# --failed: creates file containing records that failed

MaskPrimers.py extract -s sample_corrected_filtered.fastq --start 0 --len 28 --pf BARCODE --mode cut --failed --log MP.log --nproc 16 |& tee -a ${BaseFolder}/warpt/pipeline.log

# Count Ns per line
cat sample_corrected_filtered_primers-pass.fastq | awk 'NR%4==2' | grep -o -n 'N' | cut -d : -f 1 | uniq -c > QC/ns.txt

# Q average score per read
cat sample_corrected_filtered_primers-pass.fastq | perl -ne 'chomp;print;<STDIN>;<STDIN>;$_ = <STDIN>;map{ $score += ord($_)-33} 
split(""); print " " .($score/length($_))."\n";$score=0;' > QC/qscore.txt

## QC Plots: N number and Q score
Rscript /broad/vangalenlab/ainciburu/TCRseq_scripts/QCplots_preFiltering.r $BaseFolder

# 2) Quality filter
	# -q: quality score threshold
	# -s: input fastq	
	# --failed: creates file containing records that failed

FilterSeq.py quality -q 25 -s sample_corrected_filtered_primers-pass.fastq --failed --log FS.log --nproc 16 |& tee -a ${BaseFolder}/warpt/pipeline.log


# 3) Cluster sequences with same BC+UMI => discard products of barcode swapping
	# -s: input sequence
	# -f: annotation field used for grouping
	# -k: output field name
	# --ident: sequence identity threshold for the uclust algorithm
	##** need to download usearch binary from http://www.drive5.com/usearch/download.html, rename as "usearch" and give permissions

ClusterSets.py set -s sample_corrected_filtered_primers-pass_quality-pass.fastq -f BARCODE -k CLUSTER --log CS.log --failed --ident 0.9 --exec /home/unix/maincibu/.conda/envs/presto/usearch |& tee -a ${BaseFolder}/warpt/pipeline.log

# Join BARCODE and CLUSTER header tags
ParseHeaders.py merge -s sample_corrected_filtered_primers-pass_quality-pass_cluster-pass.fastq -f BARCODE CLUSTER -k BARCODE_CLUSTER --delim "|" "=" "_" |& tee -a ${BaseFolder}/warpt/pipeline.log

# Extract BC + UMI + cluster
cat sample_corrected_filtered_primers-pass_quality-pass_cluster-pass_reheader.fastq | awk 'NR%4==1' | awk 'BEGIN{FS="="};{print $4};END{}' > "${BaseFolder}/fastq_processed/BCSeq_corrected_filtered_qfiltered_cluster.txt" 

## QC plots: cluster proportion and ratio
Rscript /broad/vangalenlab/ainciburu/TCRseq_scripts/QCplots_clustering.r $BaseFolder

## Count reads passing clustering filters
#nreads=$(cat ${BaseFolder}/fastq_processed/BCSeq_corrected_filtered_qfiltered_cluster_masked.txt | grep -v "NNNNNNNNNNNNNNNNNNNNNNNNNNNN" | wc -l)
#echo -e "Reads passing clustering filters: ${nreads} \n" |& tee -a ${BaseFolder}/warpt/pipeline.log

# 4) Determines consensus for each BC/UMI
	# -s: input fastq	
	# --bf: tag by which to group sequences
	# -n: minimun number of sequences required to define a consensus
	# --maxerror: calculate error rate (number of missmatches) per consensus and remove groups exceding the given value
	# --maxgap: frequency of allowed gap values for each position. Positions exceeding the threshold are deleted from the consensus

BuildConsensus.py -s sample_corrected_filtered_primers-pass_quality-pass_cluster-pass_reheader.fastq  --bf BARCODE_CLUSTER -n 3 --maxerror 0.50 --maxgap 0.5 --outname consensus --log BC.log --failed --nproc 16 |& tee -a ${BaseFolder}/warpt/pipeline.log

nreads=$(cat  BC.log | grep "SEQCOUNT" | awk -F " " '{sum+=$2};END{print sum}')   ## count number of sequences reviewed

echo -e "Reads used for consensus building: ${nreads} \n" |& tee -a ${BaseFolder}/warpt/pipeline.log

# convert consensus fastq into fasta
paste - - - - < consensus_consensus-pass.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > file.fa


# 5) Align consensus with igblast => assign VDJ genes
	# --loci tr: look for T cell receptor
	# --format blast: output format
	# -b: IgBLAST database directory

AssignGenes.py igblast -s file.fa --organism human --loci tr --format blast -b /home/unix/maincibu/.conda/envs/presto/share/igblast |& tee -a ${BaseFolder}/warpt/pipeline.log

# 6) Create database file to store alignment results
	# -i: alignment output file 
	# -s: consensus file
	# -r: directory to alginment reference sequences
	# --extended: include additional aligner specific fields in the output
	# --partial: include incomplete V(D)J alignments

MakeDb.py igblast -i file_igblast.fmt7 -s file.fa -r /home/unix/maincibu/.conda/envs/presto/share/germlines/imgt/human/vdj/ --log MDB.log --extended --failed --partial |& tee -a ${BaseFolder}/warpt/pipeline.log

# sum up consensus log
ParseLog.py -l BC.log -o stats.log -f BARCODE SEQCOUNT CONSCOUNT ERROR |& tee -a ${BaseFolder}/warpt/pipeline.log


### Results analysis
cp file_igblast_db-pass.tsv ${BaseFolder}/results/
cp stats.log ${BaseFolder}/results/
cp ./QC/BC_UMI_cluster_metrics.txt ${BaseFolder}/results/

if [ -f $scRNAannotation ]; then
	cp ${scRNAannotation} ${BaseFolder}/results/
fi

Rscript sum_up_results.r BaseFolder mysample MinProportion MinRatio scRNAannotation


if [ -f $scRNAannotation ]; then
	Rscript plots.r
fi







