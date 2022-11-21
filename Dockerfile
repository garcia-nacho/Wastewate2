FROM garcianacho/fhibase:v1
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"
ENV qual 15
ENV noise 0.15
USER docker
RUN cd /home/docker \
    && wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /home/docker/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \
    && conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install ivar \
    && conda install -c bioconda samtools\
    && conda install -c bioconda seqkit \
    && conda update -n base -c defaults conda \
    && conda install -c bioconda bedtools \
    && conda install -c bioconda nextalign \
    && conda install -c bioconda bowtie2 \
    && conda install -c bioconda minimap2 \
    && conda create -n nextclade \
    && ln -s /home/docker/miniconda3/lib/libcrypto.so.1.1 /home/docker/miniconda3/lib/libcrypto.so.1.0.0    
RUN /bin/bash -c ". activate nextclade && \
    conda install -c bioconda nextclade"
USER root
RUN RUN apt-get update \
    && apt get install -y --no-install-recommends \
    libcairo2-dev \
    libclang-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libgdal-dev \
    libgeos-dev \
    libharfbuzz-dev \
    libjpeg-dev \
    libproj-dev \
    libpng-dev \
    libpq-dev \
    libsodium-dev \
    libssl-dev \
    libtiff5-dev \
    libudunits2-dev \
    libx11-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*    
RUN Rscript -e "install.packages(c('doSNOW', \
'progress','foreach','parallel', 'seqinr', 'doParallel', \
 'ggplot2',  'reshape2', 'ggpubr', 'readxl','tidyverse','writexl', 'remotes', 'data.table','digest'))"
RUN Rscript -e "remotes::install_github('davidsjoberg/ggsankey')"
ENV start=1250
ENV end=2250
ENV M=1300
ENV m=500
ENV mode=i
ENV trim=0
USER docker
RUN conda create -n cutadaptenv cutadapt
USER root
RUN mkdir -p /Data /home/docker/CommonFiles
COPY CommonFiles/ /home/docker/CommonFiles/
RUN chmod -R +rwx /home/docker/CommonFiles/* \
    && chmod 777 /Data 
USER docker
WORKDIR /Data
CMD ["sh", "-c", "/home/docker/CommonFiles/WWAnalysis.sh ${qual} ${noise} ${start} ${end} ${m} ${M} ${mode} ${trim}"]
