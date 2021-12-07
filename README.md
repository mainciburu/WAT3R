# WARPT
Workflow for Association of Receptor Pairs from TREK-seq


### Set up and run Docker container
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