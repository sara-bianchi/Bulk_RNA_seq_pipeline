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
    - cfr_comparisonName: variable number of columns needed if you want to perform the complete analysis, including differential gene expression analysis. Each of these columns must start with cfr_ and is followed by a unique identifier of the comparison. In each of these columns, you must set up the comparisons among two groups of samples by indicating the name of the two groups. Two samples belonging to the same group must have the same name in this column. The name assigned could be the biological condition, but also any other name. If you have multiple conditions and want to compare only two sof them, simply write the group names for the two conditions of interest and "no" in all the others.
- the settings.xlsx file, which specifies the steps of the analysis that you want to perform and the parameters needed to run each step. Two examples of this file are provided. Similarly to the table.xlsx also for this file the path must be fixed as /home/shared_folder/settings.xlsx, meaning that you cannot change the name of this file or move it in subdirectories. This file is made up of at least three columns:
    - the first contains a unique identifier of the parameters which must never be changed
    - the second contains a brief explenation about how to compile this parameter
    - the other ones contain the information needed for this specific parameter and can be edited. All the information must be added in order without leaving intermediate empty columns: for example, if only one info is needed (ex: alignment = TRUE/FALSE) the first column must be used, if more you should start top compile from the first column, then the second, and so on. The parameters which could need more columns are heatmaps, boxplots and added_genes, which require to insert a list of gene (one for each column), or GSEA and pathways which similarly require a list of GO terms or KEGG pathways of interest. In this file it is also required to add the paths of the genome information files. All the possible parameters are provided, but once you have chosen with the TRUE/FALSE parameters to selectively run only some specific steps of your analysis, all the parameters related to deselected steps can be left empty (they will just be ignored by the pipeline). For example if you do not need to run the trimming step and you set it to FALSE, all the parameters related to this step, like fastqc, adapter, adapter2, clip_5_1, clip_5_2, clip_3_1 and clip_3_2 can be left empty.

# Running the pipeline
Once all the files have been put in the working directory and all the paths have been indicated in the settings and table files, open your terminal to run the analysis.
First you have to create the connection between your folder and the container, and to do this you have to run the following command:

docker run -d -v [Path/To/Your/Folder]:/home/shared_folder --name [container_name] hedgelab/bulk_image:image4

what you need to edit is everything between squared brackets [], adding the path of your folder and assigning a unique name to your container. If you have not downloaded the image yet on your local Docker repository it will be directly downloaded from the Docker Hub. Once this download has been done the creation of futher containers for other analyses will be faster.
The next step consists in launching an Rscript which is found direclty in the downloaded image, and it is also provided in the GitHub in case you want to edit and replace it. This script can performe all the steps of the analysis indicated in the settings file using the files provided in the table. It can then save all the outputs within a folder named Output_yymmdd in your working directory. To run this script you just have to run the following command:

docker exec -it [container_name] Rscript /home/Rscript.R

The only part that you need to edit is the container_name which should match the one assigned in the previous step.

# Output interpretation
In your working directory it will be created a new folder system:
Output_yymmdd
- counts = this folder contains one directory per sample named with the corresponding sample_ID, and it will be created if alignment = TRUE. Each sample folder contains ReadsPerGene.out.tab counts file obtained from the alignment
- genome = this folder contains the indexed genome. It will be created if indexing = TRUE
- trimmed = this folder contains trimmed fastq files and the corresponding report. It will be created if trimming = TRUE
- Plots
  	- Analysis
  		- Boxplots = boxplots of the selected genes (in the settings file), normalizing counts for the sizes of the samples
  		- Heatmaps = clustered heatmaps of z-score of TPM normalized counts for: the whole dataset, the selected genes in the settings file and the significative genes for each comparison (based on the logFC and pvalue threshold set in the settings file)
  		- Modules = results of WGCNA analysis: violin plot and heatmap for the first principal component of each module in each sample, and scale free topology index plot.
  	  		- GO = GO for the selected ontology of the genes belonging to each module
  	  		- Heatmaps = clustered heatmap of the z-score of TPM normalized counts for the genes belonging to each module
  	  	- Volcano_plots = volcano plots for each comparison (significative genes for logFC and pvalue are highlighted)	
  	- GSEA
  	  	- dotplots = dotplots for the GSEA analysis of the selected ontology for each comparison
  	  	- enrichment_plots = enrichment plots and heatmap of the enrichment core of genes of the selected GO term for each comparison
  	  	- pathways = pathways plots of the selected pathways id where for each comparison, the logFC of significant genes (pvalue < selected threshold indicated in the settings file)
  	- QC
  	  	- Biotypes = percentage of counts of mitochondrial genes in each sample + biotypes distribution + biotype distribution excluding protein coding genes (plots obtained from row data before filtering)
  	  	- Chromosomes = heatmap of gene expression for each chromosome after filtering and Deseq2 vsd correction
  	  	- Correlation = heatmap of Pearson correlation generated for the top 1000 variable genes before and after batch effect correction (after filtering and Deseq2 vsd correction)
  	  	- PCA = PCA analysis was performed on filtered and Deseq2 vsd normalized counts before and after batch effect correction. For either condition is represented: the scree plot, the PC1-PC2 plot, the plots of the combination of the first n components (where n is set in the settings file), the loading plot for the first n components where the top 10 and bottom 10 genes are highlighted.
  	  	- QC = counts distribution (number of counts for each sample), library complexity, gene detection and variance stabilization plot (used for vsd normalization)
  	  	- Variability = plots of log2FC of each sample in relation to the avarge of the replicates for that condition
- RData = environment.RData file with the R data saved at the end of the analysis. Datax.Rdata files where x is the comparison name: for each comparison it contains the differential gene expression analysis object obtained from limma (diff), the counts matrix of TPM normalized counts considered significant for logFC and pvalue (counts), the gene ste enrichment analysis object obtained from ClusterProfiler for the ontology selected (gse) and keg pathways (kk) and the corresponding ordered results tables (res_gsea and resk),.
- Tables = folders with csv files
  	- Counts = raw counts filtered for protein coding genes and genes with a minimum level of expression (at least 3 counts in at least n samples where n is the minimum number of replicates per condition) with gene annotations + TPM normalized counts with the corresponding WGCNA module indicated in the “colors” column.
  	- DGEs = differential gene expression analysis for  each comparison (generated with limma)
  	- GSEA = GSEA results files for kegg pathway and for the selected ontology for each comparison (generated with ClusterProfiler)
  	- PCA = PCA loadings for the first 1000 most variable genes before and after the batch effect correction for n PCA components where n = number of samples

# Working with modified genomes
If your the genome of your biological material has undergone any gene editing process and you want to monitor transgene expression with this analysis you have to perform an additional step. Before running the indexing you have to modify the FASTA and GTF file by adding the information of your transgenes. 

- first create a FASTA file with the first line containing > gene_name (name of your transgene), and the second with the sequence of your transgene
- create a GTF file with 8 tab separated columns containing:
  - gene_name
  - exon name
  - chromosome name (for ENSEMBL annotation only the number is needed, for UCSC you need to add also the prefix chr, ex: 1, chr1)
  - length in bases of your gene
  - a separator dot (.)
  - the strand (+/-)
  - a separator dot (.)
  - a column like the following in which you have to insert the information of your gene within the "". In order to be recognised by the following steps of the pipeline gene_type must be set to protein_coding: gene_id ""; transcript_id ""; gene_name ""; gene_biotype "";
- add the FASTA file to your original FASTA file as: cat original.fa transgene.fa > final.fa
- add the GTF file to your original GTF file as: cat original.gtf transgene.gtf > final.gtf
You can either repeat these steps for all your transgenes individually or create a single FASTA and a single GTF with different lines, each one related to a specific transgene. Once you have performed these steps you can proceed with the pipeline: remember to set modified_genome to TRUE in the settings file, and to add the list of the intserted transgene names in the added_genes parameter. These names must match the gene_name parameter used in the GTF.


  

