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
RUN Rscript -e "install.packages(c('doSNOW', \
'progress','foreach','parallel', 'seqinr', 'doParallel', \
 'ggplot2',  'reshape2', 'ggpubr', 'readxl','tidyverse','writexl', 'remotes', 'data.table','digest', 'BiocManager'))"
RUN Rscript -e "remotes::install_github('davidsjoberg/ggsankey')"
RUN rm /usr/bin/gcc /usr/bin/gcc-ar /usr/bin/gcc-nm /usr/bin/gcc-ranlib \
    && ln /usr/bin/gcc-nm-9 /usr/bin/gcc-nm \
    && ln /usr/bin/gcc-ar-9 /usr/bin/gcc-ar \
    && ln /usr/bin/gcc-9 /usr/bin/gcc \
    && ln /usr/bin/gcc-ranlib-9 /usr/bin/gcc-ranlib 
RUN Rscript -e "BiocManager::install('msa')"
ENV start=1250
ENV end=2250
ENV M=1300
ENV m=500
ENV poi="0"
ENV mode=d
ENV trim=0
USER docker
RUN conda create -n cutadaptenv cutadapt
USER root
RUN Rscript -e "install.packages('phylotools')"
RUN Rscript -e "install.packages(c('umap','plotly','htmlwidgets'))"
RUN wget https://github.com/jgm/pandoc/releases/download/3.1.1/pandoc-3.1.1-1-amd64.deb && dpkg -i pandoc-3.1.1-1-amd64.deb
RUN Rscript -e "install.packages('uwot'))"
RUN mkdir -p /Data /home/docker/CommonFiles
COPY CommonFiles/ /home/docker/CommonFiles/
RUN chmod -R +rwx /home/docker/CommonFiles/* \
    && chmod 777 /Data 
USER docker
WORKDIR /Data
CMD ["sh", "-c", "/home/docker/CommonFiles/WWAnalysis.sh ${qual} ${noise} ${start} ${end} ${m} ${M} ${mode} ${trim} ${poi}"]
