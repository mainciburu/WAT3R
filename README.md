# WARPT | Workflow for Association of Receptor Pairs from TREK-seq

This analysis pipeline is for data from the T-cell Receptor Enrichment to linK clonotypes (TREK-seq) protocol as reported by DePasquale et al, [bioRxiv 2021](https://www.biorxiv.org/content/10.1101/2021.12.01.470599v1). The protocol uses 10x 3' v3 or v3.1 scRNA-seq cDNA as input and recovers *TRAV* and *TRBV* variable regions, which make up the &alpha; and &beta; chain of the T-cell receptor (TCR). Sequencing is performed on the Illumina MiSeq set to 28 bp for Read 1 (cell barcode + UMI) and 150 bp for Index 1 (*TRAV* or *TRBV*). Demultiplexing is done using Illumina `bcl2fastq` with the options `--use-bases-mask Y28,I150 --barcode-mismatches 0,0 --create-fastq-for-index-reads` and a SampleSheet with 150xN as the index sequence (not provided). After demultiplexing the sequencing data, this pipeline performs downstream analyses, including alignment, quality filters and generating a results table with cell barcodes and TCR assignments.


### Prerequisites
Since the software is run in a self-contained Docker image, you need [Docker](https://www.docker.com). There are no other requirements.


### Installation of Docker image
Either pull the image from Docker hub or build from source as follows.

1. Clone Github repository to your local disk.
2. From the top-level folder, build the docker image:
```
docker build -t WARPT .
```
3. Start the container (replace `<source>` with a local folder of choice):
```
docker run -it \
	--name WARPT \
	--mount 'type=bind,source=<source>,target=/warpt_wd' \
	--workdir '/warpt_wd' \
	WARPT \
	bash
```


### Test run
From within the Docker container, this should show the help menu:
```
warpt -h
```
Perform test run as follows. If you run from a folder that is mounted to your local device, the results are more readily accessible.
```
warpt -b /usr/local/test/data/BCseq_sub1.fastq.gz -t /usr/local/test/data/TCRseq_sub1.fastq.gz
```
The results will be written to two new folders, **fastq_processed** with intermediate analysis files and **warpt**, which contains:
- QC/QCplots_preFiltering.pdf, which shows the read quality scores and the masked bases.
- QC/QCplot_clusters.pdf, which shows the proportion of reads assigned to the most highly ranked TCR cluster (x) vs. the ratio of the reads in the first over the second ranked TCR cluster (y) vs. the number of reads per BC-UMI (color). The number in the upper right quadrant shows the proportion of reads that is maintained with the default quality thresholds (`-p` and `-r` parameters for `analyze-results`).
- sample_igblast_dp-pass.tsv, which is table of alignments for each transcript.


Continue downstream analyses as follows. To include analyses per cluster and generate additional plots, use `-a` with an txt file that contains columns with the *cell barcode* and *cluster* (e.g. cell type annotation) from accompanying scRNA-seq.
```
analyze_results -t /warpt_wd/warpt/sample_igblast_db-pass.tsv -s test -d /warpt_wd/ -a /usr/local/test/data/BM191119_Groups.txt
```
This will create a new folder named **results**:
- barcode_results.csv: final results table, summarized per cell barcode.
- barcode_UMI_results.csv: final results table, summarized per UMI.
- plots/...
- plots/...
- plots/...
- plots/...


### Analyze your own data
1. Copy fastq.gz files into the mounted working directory.
2. Run the alignment using `warpt` and downstream analyses using `analyze_results`.
