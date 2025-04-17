FROM ubuntu:24.04
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get -y update && apt-get -y upgrade
RUN apt-get -y install wget
RUN apt-get -y install cutadapt
RUN apt-get -y install fastqc
RUN apt-get -y install trim-galore
RUN apt-get install xxd
RUN apt-get -y install gcc
RUN apt-get -y install g++
RUN apt-get -y install make
RUN apt-get install zlib1g-dev
RUN apt-get install unzip
RUN wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz
RUN tar -xzf 2.7.11b.tar.gz
RUN cd STAR-2.7.11b/source && make STAR
RUN cp STAR-2.7.11b/source/STAR /
RUN apt-get -y update && apt-get -y upgrade && apt-get -y install r-base
RUN apt-get -y install gdebi-core
RUN apt-get update
RUN apt-get install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg \
    lsb-release
RUN apt-get update && apt install -y libudunits2-dev libgdal-dev
RUN apt-get update
RUN apt-get -y install gfortran
RUN apt-get install -y libharfbuzz-dev
RUN apt-get install -y libfribidi-dev
RUN apt-get install -y libcairo2-dev
RUN apt-get install -y cmake
RUN apt-get install -y libgmp-dev
RUN apt-get install -y libmpfr-dev
RUN  Rscript -e 'install.packages(c("ggplot2", "ggthemes", "reshape2", "ggplotify"), dependencies = TRUE)'
RUN Rscript -e 'install.packages(c("dplyr", "tidyr", "stringr", "viridis", "tidyverse", "ggsignif", "compareGroups", "dbscan", "BiocManager", "tibble", "ggrepel", "igraph", "readxl", "pals"), dependencies = TRUE)'
RUN  Rscript -e 'BiocManager::install(c("limma", "ensembldb", "topGO", "AnnotationHub", "AnnotationFilter", "ComplexHeatmap", "ADImpute", "impute", "preprocessCore", "DESeq2", "pheatmap", "Glimma", "edgeR", "scran", "fgsea", "BiocGenerics", "clusterProfiler", "enrichplot", "pathview"), update = TRUE, ask = FALSE)'
RUN  Rscript -e 'install.packages(c("RNAseqQC", "WGCNA"), dependencies = TRUE)'
RUN Rscript -e 'BiocManager::install(c("org.Mm.eg.db", "org.Hs.eg.db"), update = TRUE, ask = FALSE)'
RUN Rscript -e 'install.packages(c("dplyr", "tidyr", "stringr", "viridis", "tidyverse", "ggsignif", "compareGroups", "dbscan", "BiocManager", "tibble", "ggrepel", "igraph", "readxl", "pals"), dependencies = TRUE)'
RUN apt-get -y install pip
RUN pip install argparse --break-system-packages
RUN pip install numpy --break-system-packages
RUN pip install gtftools --break-system-packages
COPY Rscript.R /home/Rscript.R
EXPOSE 8787
ENTRYPOINT ["tail"]
CMD ["-f","/dev/null"]