# Bulk_RNA_seq_pipeline
Automated pipeline for the analysis of bulk gene expression data

# Installation and configuration
All the packages required to run the pipeline have already been installed within a Docker image available on our repository as hedgelab/bulk_image:image4.
The use of Docker is recommended to ensure reproducibility of the results. In case you don't have Docker installed on your computer, you can install it following the instructions: https://docs.docker.com/get-started/get-docker/
Otherwise, it is possible to recreate the working environment following the installation commands in the provided Dockerfile (optimized for a bash shell).

Ensure that the RAM and CPU allocated to Docker are sufficient to run the analysis. If you wish to run only the analysis of already aligned counts, 10 GB of RAM will be enough. Else, if you need to perform either trimming or alignment, at least 30 GB is required. In case you do not have an already-indexed genome, you will need to run the indexing step, which requires at least 100 GB of RAM. To set the proper RAM limits for your analysis, some versions of Docker Desktop allow you to directly edit these parameters in the settings, while others require the use of an external settings file. In the latter case, create a .wslconfig file in the User/user_name folder of your computer. You can follow the example template provided. After making these changes, restart the Docker engine to apply them effectively. 

Once Docker is installed and properly configured, activate the Docker engine directly opening the Docker Desktop application or by running on your terminal:
start docker

# Preparation of the input files
The first step is to choose a working directory on your computer. This directory will be the one shared with the Docker container, so anything needed for the analysis needs to be located here, either directly or by organizing your input files in subdirectories. Also the final output will be written in the same directory. The container will be connected to this directory and will recognise it from a default internal path assigned when the container is created, which is /home/shared_folder. This means that to make your files recognisable within the container, all the subsequently required paths need to start with the same prefix /home/shared_folder/ folowed by the internal paths within your working directory. For instance, if a fastq file is saved in D:/Data/Working_dir/Fastq/sample1.fastq.gz, and I decide to use Working_dir as my working directory, the path of this file will be /home/shared_folder/Fasq/sample1.fastq.gz.

In your working directory you need to have:
- your input files, either zipped fastq raw files or STAR aligned counts
- genome information:
    -  the GTF (.gtf) and FASTA (.fa or .fasta) files if you need to perform the indexing. These files can be easily downloaded from the ENSEMBL site: for the human genome, the latest version released is GRCh38.113, and you can directly download the FASTA at this link https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz and the gtf at this link https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz. For mouse genome the latest version is GRCm39.113 and you can find the FASTA at https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz and the gtf at https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz.
    -  the genome folder if you already have an indexed genome and you need to perform the alignment (the use of an already indexed genome is possible only in the case it has been created with STAR 2.7.11b, which is the same version used for the subsequent alignment)
    -  the geneInfo.tab and gene_length.txt files if you only need to run the analysis from aligned counts. In this GitHub you can already find some examples of these files for human (GRCh38.113) and mouse (GRCm39.112). Else, when you have run the full pipeline once and you just need to repeat the last part of the analysis without the alignment, you will find these files in the output (gene_length.txt directly in the Output_yymmdd folder, and geneInfo.tab in the Output_yymmdd/genome folder)
- the table.xlsx file with all the information of your samples. You can find two examples of this file in the provided folders. The name of this file must always be table.xlsx and must be stored directly in the working directory (not in subfolders). This is important since the pipeline directly llok for this file as /home/shared_folder/table.xlsx. This file is devided in different columns, each with a specific name needed to recognise the corresponding information:
    - sample_ID: contains a unique identifier of each sample
    - condition: contains the information related to the biological condition of each sample (ex: control, treated). Two samples in the same condition must have the same name in this column
    - replicate: biological replicate identifier, needed to correct, if specified, for the batch effect.
    - fastq: path of the fastq files within the shared folder (only for single-end reads), if you want to start from the alignment
    - fastq1: path of fastq files read1 within the shared folder (only for paired-end reads), if you want to start from the alignment
    - fastq2: path of fastq files read2 within the shared folder (only for paired-end reads), if you want to start from the alignment
    - counts: path of the STAR aligned counts as SampleReadsPerGene.out.tab, if you want to start after the alignment
    - cfr_comparisonName: variable number of columns needed if you want to perform the complete analysis, including differential gene expression analysis. Each of these columns must start with cfr_ and is followed by a unique identifier of the comparison. In each of these columns, you must set up the comparisons among two groups of samples by indicating the name of the two groups. Two samples belonging to the same group must have the same name in this column. The name assigned could be the biological condition, but also any other name. If you have many conditions and you want to compare only two specific conditions, you just have to write the group names for the two conditions of interest and "no" in all the others.
- the settings.xlsx file, which specifies the steps of the analysis that you want to perform and the parameters needed to run each step. Two examples of this file are provided. Similarly to the table.xlsx also for this file the path must be fixed as /home/shared_folder/settings.xlsx, meaning that you cannot change the name of this file or move it in subdirectories. This file is made up of at least three columns:
    - the first contains a unique identifier of the parameters which must never be changed
    - the second contains a brief explenation about how to compile this parameter
    - the other ones contain the information needed for this specific parameter and can be edited. All the information must be added in order without leaving intermediate empty columns: for example, if only one info is needed (ex: alignment = TRUE/FALSE) the first column must be used, if more you should start top compile from the first column, then the second, and so on. The parameters which could need more columns are heatmaps, boxplots and added_genes, which require to insert a list of gene (one for each column), or GSEA and pathways which similarly require a list of GO terms or KEGG pathways of interest.




