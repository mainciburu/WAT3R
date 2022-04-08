FROM continuumio/miniconda3:4.10.3

# Set the working directory to be used when the docker gets run
WORKDIR /

# Install required software
RUN apt-get update && apt-get install -y gnupg && apt-get install -y wget

# Install R version 4.1.2
#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
RUN echo "deb http://cloud.r-project.org/bin/linux/debian buster-cran40/" >> /etc/apt/sources.list
RUN apt-get update && apt install -y -t buster-cran40 r-base-core && apt install -y -t buster-cran40 r-base

# Install python modules
RUN conda install python=3.6 \
&& pip install presto==0.6.2 changeo==1.1.0 \
&& conda install -c bioconda -c conda-forge  umi_tools=1.1.2

# Install R packages
RUN R -e 'install.packages(c("dplyr", "ggplot2", "RColorBrewer", "cowplot", "reshape2", "gridExtra", "stringr", "mclust", "pheatmap", "BiocManager", "e1071", "pals"), repos = c(CRAN = "http://cran.rstudio.org"))' \
&& R -e 'install.packages("http://cran.r-project.org/src/contrib/Archive/rjson/rjson_0.2.20.tar.gz", repo = NULL)' \
&& R -e 'BiocManager::install("ComplexHeatmap")'

# Get 10x BC whitelist
RUN wget http://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz -P /usr/local/

# Get usearch executable
RUN wget http://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz -P /usr/local \
&& gzip -d /usr/local/usearch11.0.667_i86linux32.gz \
&& mv /usr/local/usearch11.0.667_i86linux32 /usr/local/bin/usearch

#### Configure IgBLAST

RUN wget http://bitbucket.org/kleinstein/immcantation/raw/9e5f6fb95edda9901238abc28da8a29948f9de82/scripts/fetch_igblastdb.sh -P /usr/local/bin/ \
&& wget http://bitbucket.org/kleinstein/immcantation/raw/9e5f6fb95edda9901238abc28da8a29948f9de82/scripts/fetch_imgtdb.sh -P /usr/local/bin/ \
&& wget http://bitbucket.org/kleinstein/immcantation/raw/9e5f6fb95edda9901238abc28da8a29948f9de82/scripts/clean_imgtdb.py -P /usr/local/bin/ \
&& wget http://bitbucket.org/kleinstein/immcantation/raw/9e5f6fb95edda9901238abc28da8a29948f9de82/scripts/imgt2igblast.sh -P /usr/local/bin/ \
&& chmod -R 777 /usr/local/bin/ \
&& wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.17.1/ncbi-igblast-1.17.1-x64-linux.tar.gz -P /usr/local/ \
&& tar -zxf /usr/local/ncbi-igblast-1.17.1-x64-linux.tar.gz -C /usr/local/ \
&& chmod -R 777 /usr/local/ncbi-igblast-1.17.1 \
&& rm /usr/local/ncbi-igblast-1.17.1-x64-linux.tar.gz \
&& cp /usr/local/ncbi-igblast-1.17.1/bin/* /usr/local/bin/ \
&& fetch_igblastdb.sh -o /usr/igblast \
&& cp -r /usr/local/ncbi-igblast-1.17.1/internal_data /usr/igblast \
&& cp -r /usr/local/ncbi-igblast-1.17.1/optional_file /usr/igblast \
&& fetch_imgtdb.sh -o /usr/germlines/imgt \
&& imgt2igblast.sh -i /usr/germlines/imgt -o /usr/igblast \
&& chmod -R 777 /usr/igblast \
&& chmod -R 777 /usr/germlines

#### Add WAT3R scripts and test data
COPY scripts/* /usr/local/bin/
COPY testdata/* /usr/local/testdata/
RUN chmod -R 777 /usr/local/bin/ \
&& chmod -R 777 /usr/local/testdata

#USER vscode

#ENTRYPOINT ["sleep", "infinity"]







