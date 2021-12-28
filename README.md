# WAT3R | Workflow for the Association of T-cell receptors from 3' single-cell RNA-seq

This analysis pipeline is for data from the T-cell Receptor Enrichment to linK clonotypes (TREK-seq) protocol as reported by DePasquale et al, [bioRxiv 2021](https://www.biorxiv.org/content/10.1101/2021.12.01.470599v1). The protocol uses 10x 3' v3 or v3.1 scRNA-seq cDNA as input and recovers *TRAV* and *TRBV* variable regions, which make up the &alpha; and &beta; chain of the T-cell receptor (TCR). Sequencing is performed on the Illumina MiSeq set to 28 bp for Read 1 (cell barcode + UMI) and 150 bp for Index 1 (*TRAV* or *TRBV*). Demultiplexing is done using Illumina `bcl2fastq` with the options `--use-bases-mask Y28,I150 --barcode-mismatches 0,0 --create-fastq-for-index-reads` and a SampleSheet with 150xN as the index sequence (not provided). After demultiplexing the sequencing data, this pipeline performs downstream analyses, including alignment, quality filters and generating a results table with cell barcodes and TCR assignments.

WAT3R can be run in a [Docker](https://www.docker.com) container or as a workflow on [Terra](https://app.terra.bio). Terra provides access to Google Cloud computing resources through a simple web-based user interface for non-coding biologists. Below we describe how to run WAT3R using Docker.


## Prerequisites
Since the software is run in a self-contained Docker image, you need Docker. There are no other requirements.


## Installation of Docker image
Either pull the image from Docker hub or build from source as follows.

1. Clone Github repository to your local disk.
2. From the top-level folder, build the docker image:
```
docker build -t mainciburu/wat3r .
```
3. Start the container (replace `<source>` with a local folder of choice):
```
docker run -it \
	--name WAT3R \
	--mount 'type=bind,source=<source>,target=/workdir' \
	--workdir '/workdir' \
	mainciburu/wat3r \
	bash
```
4. If you succesfully entered the container, this should show the help menu:
```
wat3r -h
```


## Test run
The analysis is split into two parts: `wat3r` and `downstream`. 

#### `wat3r`
This command requires two fastq files: one (`-b`) with a 16 bp cell barcode (CB) and 12 bp unique molecular identifier (UMI) and one (`-t`) with the TCR sequences. If you run from a folder that is mounted to your local device, the results will be more readily accessible.
```
wat3r -b /usr/local/testdata/BCseq_test.fastq.gz -t /usr/local/testdata/TCRseq_test.fastq.gz
```
The results will be written to two new folders, **fastq_processed** with intermediate analysis files and **wat3r**, which contains:
- sample_igblast_dp-pass.tsv, which is table of alignments for each transcript.
- QC/BC_UMI_cluster_metrics.txt, which is a table of metrics for each of the BC-UMI clusters.
- stats.log, which contains error rates for consensus building.- - wat3rMetrics.txt, which contains the number of reads after each filtering step.
The **watr** folder also contains two plots:
- QC/QCplots_preFiltering.pdf, which shows the read quality scores and the masked bases.
- QC/QCplot_clusters.pdf, which shows the proportion of reads assigned to the most highly ranked TCR cluster (x) vs. the ratio of the reads in the first over the second ranked TCR cluster (y) vs. the number of reads per CB-UMI (color). The number in the upper right quadrant shows the proportion of reads that is maintained with the default quality thresholds, indicated by the red lines (`-p` and `-r` parameters for `downstream`).

#### `downstream`
Continue downstream analyses as follows. To include analyses per cluster and generate additional plots, use `-a` with an txt file that contains columns (no header) with the *cell barcode* and *cluster* (e.g. cell type annotation) from accompanying scRNA-seq.
```
downstream -u /workdir/wat3r/sample_igblast_db-pass.tsv \
	-c /workdir/wat3r/QC/BC_UMI_cluster_metrics.txt \
	-f /workdir/wat3r/wat3rMetrics.txt \
	-s /workdir/wat3r/stats.log \
	-n PB01 \
	-a /usr/local/testdata/PB01_clusters.txt
```
This will create a new folder named **downstream**:
- barcode_results.csv: final results table, summarized per cell barcode.
- barcode_UMI_results.csv: final results table, summarized per UMI.

Plots without cell annotations:
- plots/db_histograms.pdf, which shows a histogram of read numbers (>=3) per consensus sequence, and a histogram of the error rate per consensus sequence.
- plots/ReadPercentage_FilteringSteps.pdf, which shows the percent reads remaining after each filtering step.

Additional plots with cell annotations:
- plots/scRNAseq_TCRrecovery_proportions.pdf, which shows the proportion of each cell annotation for which TRA and/or TRB genes were detected
- plots/valid_reads.pdf, which shows the total number of reads assigned to TRA and TRB genes, separated by cell annotation.
- plots/CDR3_UMIcount_distribution.pdf, which shows a histogram of the number of UMI (UMI counts) assigned to TRA and TRB CDR3 sequences. Frequencies are independently counted for barcodes overlapping or not the scRNAseq dataset. 
- plots/CDR3_clones_heatmap.pdf, which shows a heatmap with the cell number matching specific TRA and TRB CDR3 sequences.
- plots/TRB_TRA_correspondence.pdf, which shows individual cell correspondence between TRA and TRB CDR3 sequences, together with cell annotation. 
- plots/TRA_TRB_clone_size.pdf, which shows the ranking of clone sizes for both TRA and TRB, expressed in total cell number and normalized cell number (i.e. % of the total cell number with TRA or TRB that belong to a clone). 
- plots/trb_top_clones.pdf, which shows the size and cell composition of the 50 biggest clones, using the TRB CDR3 sequences. Size is measured as total cell number.
- plots/trb_top_clones_norm.pdf, which shows the size and cell composition of the 50 biggest clones, using the TRB CDR3 sequences. Size is measured as normalized cell number.
- plots/TRB_clone_size_celltype.pdf, which shows the size of the clone to which each of the annotated cells belongs, separated by annotation.  
- plots/TRB_distance_heatmap.pdf, which shows a heatmap with the Hamming distance between the TRB CDR3 sequences from the 50 biggest clones. 


## Analyze your own data
1. Copy fastq.gz files into the mounted working directory.
2. Run the initial steps using `wat3r` and downstream analyses using `downstream`.
