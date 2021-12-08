# WARPT
Workflow for Association of Receptor Pairs from TREK-seq

This analysis pipeline is developed for data from the T-cell Receptor Enrichment to linK clonotypes (TREK-seq) protocol as reported by DePasquale et al, [bioRxiv 2021](https://www.biorxiv.org/content/10.1101/2021.12.01.470599v1). The protocol uses 10x 3' v3 or v3.1 scRNA-seq cDNA as input and recovers *TRAV* and *TRBV* variable regions, which make up the &alpha; and &beta; chain of the T-cell receptor (TCR). Data is generated on the Illumina MiSeq set to sequence Read 1 (28 bp, cell barcode + UMI) and Index 1 (150 bp, *TRAV* or *TRBV*). After demultiplexing the sequencing data, this pipeline performs downstream analysis, including alignment, quality filters and calling TCR variable regions. 


### Installation of Docker image
1. Clone Github repository to your local disk
2. From the top directory, build the docker image using `docker build`
3. Start the container using `docker run`


### Test run
From within the Docker container, perform test run. The results will be written to *fastq_processed* and *warpt* in your current directory (it may be convenient to execute from a mounted directory).
```
warpt -b /usr/local/test/data/BCseq_sub1.fastq.gz -t /usr/local/test/data/TCRseq_sub1.fastq.gz
```


### Test analysis
This will create a new directory *results*.
```
analyze_results -t /warpt_wd/warpt/sample_igblast_db-pass.tsv -s test -d /warpt_wd/ -a /usr/local/test/data/BM191119_Groups.txt
```