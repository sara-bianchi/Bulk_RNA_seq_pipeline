#------------------------------------------------------------------------------
#BULK RNA SEQ
#------------------------------------------------------------------------------

#packages
library(tidyr)
library(stringr)
library(grid)
library(viridis)
library(ggthemes)
library(dplyr)
library(tidyverse)
library(edgeR)
library(compareGroups)
library(dbscan)
library(reshape2)
library(scran)
library(fgsea)
library(RNAseqQC)
library(DESeq2)
library(ensembldb)
library(tibble)
library(ggplotify)
library(limma)
library(Glimma)
library(topGO)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggrepel)
library(ggsignif)
library(igraph)
library(readxl)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(WGCNA)
library(readr)

#set up environment
theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size = 15, face="bold", margin = margin(10,0,10,0)),
                  axis.text.x = element_text(angle = 45, hjust = 1)))
options(stringsAsFactors = FALSE)
update_geom_defaults("point", aes(size = 4))
set.seed(1234)

pal = c("#3283FE", "#FA0087", "#009E73", "#FBE426", "#56B4E9", "#FEAF16", "#DEA0FD", "#1CBE4F", "#F6222E", "#1CFFCE", "#325A9B", "#AA0DFE","#D55E00", "#2ED9FF", "#f0E442", "#1C8356", "#0072B2", "#CC79A7")

current_time = format(Sys.time(), "%y%m%d")
dir.create(paste0("/home/shared_folder/Output_", current_time))


#upload settings

print("uploading settings file")

settings = as.data.frame(read_excel("/home/shared_folder/settings.xlsx", col_names = FALSE))
rownames(settings) = settings$...1
settings = settings %>% dplyr::select(-c("...1", "...2"))
settings = as.data.frame(t(settings))

if(length(unique(na.omit(settings$trimming))) != 1 | as.logical(settings$trimming[1]) == "NA"){
  print("invalid argument for trimming"); stop()
} else{trimming = as.logical(settings$trimming[1])}

if(trimming == TRUE){
  if(length(unique(na.omit(settings$fastqc))) != 1 | as.logical(settings$fastqc[1]) == "NA"){
    print("invalid argument for fastqc"); stop()
  } else{fastqc = as.logical(settings$fastqc[1]);
  if(fastqc == TRUE){fqc = " --fastqc"}
  else if(fastqc == FALSE){fqc = ""}}
  if(length(unique(na.omit(settings$mode))) != 1 | !as.character(settings$mode[1]) %in% c("single", "paired")){
    print("invalid argument for mode"); stop()
  } else{mode = as.character(settings$mode[1])};
  bases = strsplit(settings$adapter[1], "")
  invalid_bases = c()
  for(i in bases[[1]]){if(!i %in% c("A", "T", "G", "C")){invalid_bases = append(invalid_bases, i)}}
  if(length(unique(na.omit(settings$adapter))) != 1 | length(invalid_bases) != 0){
    print("invalid argument for adapter"); stop()
  } else{adapter = as.character(settings$adapter[1])};
  if(length(unique(na.omit(settings$clip_5_1))) != 1 | as.integer(settings$clip_5_1[1]) == "NA"){
    print("invalid argument for clip_5_1"); stop()
  } else{clip_5_1 = as.integer(settings$clip_5_1[1])};
  if(length(unique(na.omit(settings$clip_3_1))) != 1 | as.integer(settings$clip_3_1[1]) == "NA"){
    print("invalid argument for clip_3_1"); stop()
  } else{clip_3_1 = as.integer(settings$clip_3_1[1])};
  if(mode == "paired"){
    bases = strsplit(settings$adapter2[1], "")
    invalid_bases = c()
    for(i in bases[[1]]){if(!i %in% c("A", "T", "G", "C")){invalid_bases = append(invalid_bases, i)}}
    if(length(unique(na.omit(settings$adapter2))) != 1 | length(invalid_bases) != 0){
      print("invalid argument for adapter2"); stop()
    } else{adapter2 = as.character(settings$adapter2[1])};
    if(length(unique(na.omit(settings$clip_5_2))) != 1 | as.integer(settings$clip_5_2[1]) == "NA"){
      print("invalid argument for clip_5_2"); stop()
    } else{clip_5_2 = as.integer(settings$clip_5_2[1])};
    if(length(unique(na.omit(settings$clip_3_2))) != 1 | as.integer(settings$clip_3_2[1]) == "NA"){
      print("invalid argument for clip_3_2"); stop()
    } else{clip_3_2 = as.integer(settings$clip_3_2[1])};
  }
}

if(length(unique(na.omit(settings$indexing))) != 1 | as.logical(settings$indexing[1]) == "NA"){
  print("invalid argument for indexing"); stop()
} else{indexing = as.logical(settings$indexing[1])}

if(length(unique(na.omit(settings$alignment))) != 1 | as.logical(settings$alignment[1]) == "NA"){
  print("invalid argument for alignment"); stop()
} else{alignment = as.logical(settings$alignment[1])}

if(indexing == TRUE){
  if(length(unique(na.omit(settings$gtf))) != 1 | file.exists(settings$gtf[1]) == FALSE){
    print("invalid argument for gtf"); stop()
  } else{gtf = as.character(settings$gtf[1])};
  if(length(unique(na.omit(settings$fasta))) != 1 | file.exists(settings$fasta[1]) == FALSE){
    print("invalid argument for fasta"); stop()
  } else{fasta = as.character(settings$fasta[1])};
} else{if(alignment == TRUE){
  if(length(unique(na.omit(settings$genome_dir))) != 1 | dir.exists(settings$genome_dir[1]) == FALSE){
    print("invalid argument for genom_dir"); stop()
  } else{genome_dir = as.character(settings$genome_dir[1])}};
  if(length(unique(na.omit(settings$gene_length))) != 1 | file.exists(settings$gene_length[1]) == FALSE){
    print("invalid argument for gene_length"); stop()
  } else{gene_length = as.character(settings$gene_length[1])};
  if(length(unique(na.omit(settings$gene_info))) != 1 | file.exists(settings$gene_info[1]) == FALSE){
    print("invalid argument for gene_info"); stop()
  } else{geneInfo = as.character(settings$gene_info[1]); geneInfo = read.table(geneInfo, quote = "\"", comment.char = "", skip = 1)}
}

if(length(unique(na.omit(settings$analysis))) != 1 | as.logical(settings$analysis[1]) == "NA"){
  print("invalid argument for analysis"); stop()
} else{analysis = as.logical(settings$analysis[1])}

if(alignment == TRUE){if(length(unique(na.omit(settings$mode))) != 1 | !as.character(settings$mode[1]) %in% c("single", "paired")){
  print("invalid argument for mode"); stop()
} else{mode = as.character(settings$mode[1])}}

if(analysis == TRUE){
  if(length(unique(na.omit(settings$organism))) != 1 | !as.character(settings$organism[1]) %in% c("human", "mouse")){
    print("invalid argument organism"); stop()
  } else{organism = as.character(settings$organism[1])}
  
  if(organism == "human"){
    ah_record = "AH89426";
    chromosomes = c(1:22, "X", "Y", "MT")
    org = org.Hs.eg.db;
    library("org.Hs.eg.db", character.only = TRUE);
    kegg_organism = "hsa";
    specie = "Hs"
  } else if(organism == "mouse"){
    ah_record = "AH89211";
    chromosomes = c(1:19, "X", "Y", "MT")
    org = org.Mm.eg.db;
    library("org.Mm.eg.db", character.only = TRUE);
    kegg_organism = "mmu";
    specie = "Mm"
  }
  
  if(length(unique(na.omit(settings$batch_correction))) != 1 | as.logical(settings$batch_correction[1]) == "NA"){
    print("invalid argument for batch_correction"); stop()
  } else{batch_correction = as.logical(settings$batch_correction[1])}
  
  if(length(unique(na.omit(settings$ontology))) != 1 | !as.character(settings$ontology[1]) %in% c("BP", "CC", "MF")){
    print("invalid argument ontology"); stop()
  } else{ontology = as.character(settings$ontology[1])}
  
  if(length(unique(na.omit(settings$GO_n))) != 1 | as.integer(settings$GO_n[1]) == "NA"){
    print("invalid argument for GO_n"); stop()
  } else{GO_n = as.integer(settings$GO_n[1])}
  
  if(length(unique(na.omit(settings$PCA))) != 1 | as.integer(settings$PCA[1]) == "NA"){
    print("invalid argument for PCA"); stop()
  } else{PCA = as.integer(settings$PCA[1])}
  
  if(length(unique(na.omit(settings$logFC))) != 1 | as.numeric(settings$logFC[1]) == "NA" | as.numeric(settings$logFC[1]) < 0){
    print("invalid argument for logFC"); stop()
  } else{logFC = as.numeric(settings$logFC[1])}
  
  if(length(unique(na.omit(settings$pvalue))) != 1 | as.numeric(settings$pvalue[1]) == "NA" | as.numeric(settings$pvalue[1]) > 1 | as.numeric(settings$pvalue[1]) < 0){
    print("invalid argument for pvalue"); stop()
  } else{pvalue = as.numeric(settings$pvalue[1])}
  
  go_terms = na.omit(settings$GSEA);
  invalid_go = c()
  for(i in go_terms){i_s = strsplit(i, split = ""); if(paste0(i_s[[1]][1], i_s[[1]][2], i_s[[1]][3]) != "GO:"){invalid_go = append(invalid_go, i)}}
  if(length(invalid_go) != 0){
    print("invalid argument for GSEA"); stop()
  } else{GSEA = na.omit(settings$GSEA)}
  
  kegg_terms = na.omit(settings$pathways);
  invalid_path = c()
  for(i in kegg_terms){i_s = strsplit(i, split = ""); if(paste0(i_s[[1]][1], i_s[[1]][2], i_s[[1]][3]) != kegg_organism){invalid_path = append(invalid_path, i)}}
  if(length(invalid_path) != 0){
    print("invalid argument for pathways"); stop()
  } else{pathways = na.omit(settings$pathways)}
  
}

#table

print("uploading table file")

table = read_excel("/home/shared_folder/table.xlsx")
if(alignment == TRUE){
  if(mode == "single"){
    if(!"fastq" %in% colnames(table)){
      print("fastq column is missing"); stop()
    } else{
      for(i in table$fastq){
        if(file.exists(i) == FALSE){
          print("incorrect path in fastq column"); stop()
        }
      }
    }
  } else if(mode == "paired"){
    if(!"fastq1" %in% colnames(table) | !"fastq2" %in% colnames(table)){
      print("fastq1 or fastq2 columns are missing"); stop()
    } else{
      for(i in table$fastq1){if(file.exists(i) == FALSE){
        print("incorrect path in fastq1 column"); stop()}
      }
      for(i in table$fastq2){if(file.exists(i) == FALSE){
        print("incorrect path in fastq2 column"); stop()}
      }
    }
  }
} else if(alignment == FALSE & analysis == TRUE){
  if(!"counts" %in% colnames(table)){
    print("counts column is missing"); stop()
  } else{
    for(i in table$counts){
      if(file.exists(i) == FALSE){
        print("incorrect path in counts column"); stop()
      }
    }
  }
}

#trimming
if(trimming == TRUE){
  print("starting trimming")
  
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/trimmed"));
  out_dir = paste0("/home/shared_folder/Output_", current_time, "/trimmed");
  if(mode == "single"){
    if(clip_3_1 == 0){
      clipping_R1_3 = ""
    } else {
      clipping_R1_3 = paste0(" --three_prime_clip_R1 ", as.character(clip_3_1))
    }
    if(clip_5_1 == 0){
      clipping_R1_5 = ""
    } else {
      clipping_R1_5 = paste0(" --clip_R1 ", as.character(clip_5_1))
    }
    
    for(i in 1:length(rownames(table))){fastq = table$fastq[i]; sample_ID = table$sample_ID[i]
      system2("/usr/bin/trim_galore", paste0("--adapter ", as.character(adapter),  clipping_R1_5, clipping_R1_3, " --gzip", fqc, " --basename ", sample_ID, " --output_dir ", as.character(out_dir), " ", as.character(fastq)))
    }
  } else if(mode == "paired"){
    if(clip_3_1 == 0){
      clipping_R1_3 = ""
    } else {
      clipping_R1_3 = paste0(" --three_prime_clip_R1 ", as.character(clip_3_1))
    }
    if(clip_5_1 == 0){
      clipping_R1_5 = ""
    } else {
      clipping_R1_5 = paste0(" --clip_R1 ", as.character(clip_5_1))
    }
    if(clip_3_2 == 0){
      clipping_R2_3 = ""
    } else {
      clipping_R2_3 = paste0(" --three_prime_clip_R2 ", as.character(clip_3_2))
    }
    if(clip_5_2 == 0){
      clipping_R2_5 = ""
    } else {
      clipping_R2_5 = paste0(" --clip_R2 ", as.character(clip_5_2))
    }
    for(i in 1:length(rownames(table))){fastq1 = table$fastq1[i]; fastq2 = table$fastq2[i]; sample_ID = table$sample_ID[i]
      system2("/usr/bin/trim_galore", paste0("--adapter ", as.character(adapter), " --adapter2 ", as.character(adapter2), clipping_R1_5, clipping_R1_3, clipping_R2_5, clipping_R2_3," --gzip --paired", fqc, " --basename ", sample_ID, " --output_dir ", as.character(out_dir), " ", as.character(fastq1), " ", as.character(fastq2)))
    }
  }
  print("finishing trimming")
}

#indexing
if(indexing == TRUE){
  print("starting indexing")
  
  system2("mkdir", "/home/genome")
  system2("chmod", "777 /home/genome")
  system2("/STAR", paste0("--runThreadN 4 --runMode genomeGenerate --limitGenomeGenerateRAM=150000000000 --genomeDir /home/genome --genomeFastaFiles ", fasta, " --sjdbGTFfile ", gtf));
  system2("rm", paste0("-r /home/shared_folder/Output_", current_time, "/genome"));
  system2("mv", paste0("/home/genome /home/shared_folder/Output_", current_time));
  genome_dir = paste0("/home/shared_folder/Output_", current_time, "/genome");
  system2("gtftools", paste0("-l /home/shared_folder/Output_", current_time, "/gene_length.txt ", gtf));
  gene_length = paste0("/home/shared_folder/Output_", current_time, "/gene_length.txt");
  geneInfo = read.table(paste0(genome_dir, "/geneInfo.tab"), quote = "\"", comment.char = "", skip = 1)
  
  print("finishing indexing")
}

colnames(geneInfo) = c("gene_ID", "gene_name", "gene_type")

if(analysis == TRUE){
  if(length(unique(na.omit(settings$modified_genome))) != 1 | as.logical(settings$modified_genome[1]) == "NA"){
    print("invalid argument for modified_genome"); stop()
  } else{modified_genome = as.logical(settings$modified_genome[1])}
  
  if(modified_genome == TRUE){
    added = na.omit(settings$added_genes);
    invalid_genes = c()
    for(i in added){if(!i %in% geneInfo$gene_ID){invalid_genes = append(invalid_genes, i)}}
    if(length(invalid_genes) != 0){
      print("invalid argument for added_genes"); stop()
    } else{added_genes = na.omit(settings$added_genes)};
  } else{added_genes = c()}
  
  geneInfo_ext = geneInfo %>% dplyr::filter(gene_ID %in% added_genes)
  for(i in geneInfo_ext$gene_name){geneInfo_ext$gene_name = recode(geneInfo_ext$gene_name, i = paste0(i, "_ext"))}
  
  geneInfo_s = geneInfo %>% dplyr::filter(!gene_ID %in% added_genes)
  colnames(geneInfo_s) = c("gene_ID", "gene_name", "gene_type")
  
  if(modified_genome == TRUE){geneInfo = rbind(geneInfo_s, geneInfo_ext)}
  
  box = na.omit(settings$boxplots)
  invalid_genes = c()
  for(i in box){if(!i %in% geneInfo$gene_name){invalid_genes = append(invalid_genes, i)}}
  if(length(invalid_genes) != 0){
    print("invalid argument for boxplots"); stop()
  } else{boxplots = na.omit(settings$boxplots)}
  
  heat = na.omit(settings$heatmaps)
  invalid_genes = c()
  for(i in heat){if(!i %in% geneInfo$gene_name){invalid_genes = append(invalid_genes, i)}}
  if(length(invalid_genes) != 0){
    print("invalid argument for heatmpas"); stop()
  } else{heatmaps = na.omit(settings$heatmaps)}
}

#alignment
if(alignment == TRUE){
  print("starting alignment")
  
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/counts"));
  if(mode == "single"){
    for(i in 1:length(rownames(table))){
      sample_ID = table$sample_ID[i];
      counts_dir = paste0("/home/shared_folder/Output_", current_time, "/counts/", sample_ID);
      dir.create(counts_dir);
      system2("chmod", paste0("777 ", counts_dir));
      if(trimming == FALSE){fastq = table$fastq[i]} else if(trimming == TRUE){
        fastq = paste0(out_dir, "/", sample_ID, "_trimmed.fq.gz")};
      system2("/STAR", paste0("--runThreadN 4 --genomeDir ", genome_dir, " --readFilesIn ", fastq, " --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --readFilesCommand zcat --outFileNamePrefix ", counts_dir, "/", sample_ID))}
  } else if(mode == "paired"){
    for(i in 1:length(rownames(table))){sample_ID = table$sample_ID[i];
    counts_dir = paste0("/home/shared_folder/Output_", current_time, "/counts/", sample_ID);
    if(trimming == FALSE){fastq1 = table$fastq1[i]; fastq2 = table$fastq2[i]} else if(trimming == TRUE){
      fastq1 = paste0(out_dir, "/", sample_ID, "_val_1.fq.gz");
      fastq2 = paste0(out_dir, "/", sample_ID, "_val_2.fq.gz")};
    system2("/STAR", paste0("--runThreadN 4 --genomeDir ", genome_dir, " --readFilesIn ", fastq1, " ", fastq2, " --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --readFilesCommand zcat --outFileNamePrefix ", counts_dir, "/", sample_ID))}
  }
  print("finishing alignment")
}


if(analysis == TRUE){
  
  #folders creation
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/RData"))
  
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Tables"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Tables/DGEs"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Tables/GSEA"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Tables/PCA"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Tables/Counts"))
  
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/QC"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/QC"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/Biotypes"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/Correlation"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/PCA"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/Variability"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/Chromosomes"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Boxplots"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Heatmaps"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Modules"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Modules/GO"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Modules/Heatmaps"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Volcano_plots"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/GSEA"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/GSEA/dotplots"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/GSEA/enrichment_plots"))
  dir.create(paste0("/home/shared_folder/Output_", current_time, "/Plots/GSEA/pathways"))
  
  #ANALYSIS
  
  print("starting QC")
  
  #creation of starting datasets
  meta = table %>% dplyr::select(c("sample_ID", "condition", "replicate"))
  
  if(alignment == TRUE){counts = c();
  for(i in table$sample_ID){counts = append(counts, paste0("/home/shared_folder/Output_", current_time, "/counts/", i, "/", i, "ReadsPerGene.out.tab"))};
  table$counts = counts}
  
  counts_all = geneInfo
  for(i in 1:length(rownames(table))){
    sample_i = read.table(table$counts[i], quote = "\"", comment.char = "", skip = 4) %>% dplyr::select("V1", "V4");
    colnames(sample_i) = c("gene_ID", table$sample_ID[i]);
    counts_all = counts_all %>% dplyr::left_join(sample_i, by = "gene_ID")
  }
  
  #biotypes and counts distribution
  #create meta biotypes
  meta_biotypes = meta
  counts_matrix_all = counts_all %>% dplyr::select(- colnames(geneInfo))
  counts_matrix_all = as.matrix(counts_matrix_all)
  all_counts_sum = as.data.frame(colSums(counts_matrix_all))
  colnames(all_counts_sum) = "all_genes_sum"
  all_counts_sum$sample_ID = rownames(all_counts_sum)
  ggplot(all_counts_sum, aes(x = sample_ID, y = all_genes_sum)) + geom_bar(stat = "identity") + scale_fill_manual(values = pal) + ggtitle("total number of counts")
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/QC/counts_distribution.pdf"),
         width = 14, height = 8)
  
  meta_biotypes = meta_biotypes  %>% left_join(all_counts_sum, by = "sample_ID")
  
  #mt genes
  mt_genes = c()
  if(organism == "human"){for(i in geneInfo$gene_name) {x =  grepl("MT-", i); if(x == TRUE) {mt_genes = append(mt_genes, i)}}}
  if(organism == "mouse"){for(i in geneInfo$gene_name) {x =  grepl("mt-", i); if(x == TRUE) {mt_genes = append(mt_genes, i)}}}
  mt_counts = counts_all
  mt_counts = mt_counts %>% dplyr :: filter(gene_name %in% mt_genes)
  mt_counts = mt_counts %>% dplyr :: select(- colnames(geneInfo))
  mt_counts = as.matrix(mt_counts)
  mt_counts_sum = as.data.frame(colSums(mt_counts))
  colnames(mt_counts_sum) = "mt_genes_sum"
  mt_counts_sum$sample_ID = rownames(mt_counts_sum)
  meta_biotypes = meta_biotypes  %>% left_join(mt_counts_sum, by = "sample_ID")
  meta_biotypes = meta_biotypes %>% dplyr :: mutate(mt_percentage = mt_genes_sum/all_genes_sum*100)
  ggplot(meta_biotypes, aes(x = sample_ID, y = mt_percentage, fill = condition)) + geom_bar(stat = "identity") + scale_fill_manual(values = pal) + ggtitle("percentages of mitochondrial counts")
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/Biotypes/mitochondrial_percentages.pdf"),
         width = 14, height = 8)
  
  #pointplot
  for(i in unique(geneInfo$gene_type)){
    counts_i = as.matrix(counts_all %>% dplyr::filter(gene_type == i) %>% dplyr::select(-colnames(geneInfo)));
    sum_i = as.data.frame(colSums(counts_i));
    colnames(sum_i) = "sum"
    sum_i$sample_ID = rownames(sum_i);
    sum_i = sum_i %>% dplyr::left_join(all_counts_sum, by = "sample_ID") %>% mutate(perc = sum/all_genes_sum*100) %>% dplyr::select(-c(sum, all_genes_sum));
    colnames(sum_i) = c("sample_ID", i)
    meta_biotypes = meta_biotypes %>% dplyr::left_join(sum_i, by = c("sample_ID"))
  }
  
  meta_biotypes = meta_biotypes %>% dplyr::select(- c(mt_percentage, mt_genes_sum, all_genes_sum))
  biotypes = melt(meta_biotypes)
  colnames(biotypes) = c("sample_ID", "condition", "replicate", "biotype", "percentage")
  
  for(i in unique(geneInfo$gene_type)){
    melt_i = biotypes %>% dplyr::filter(biotype == i);
    sum_i = sum(melt_i$percentage);
    if(sum_i == 0){biotypes = biotypes %>% dplyr::filter(biotype != i)}
  }
  
  ggplot(biotypes, aes(x = sample_ID, y = percentage, group = biotype, color = biotype)) +
    geom_line() + geom_point()
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/Biotypes/biotypes_percentages.pdf"),
         width = 14, height = 8)
  
  ggplot(biotypes %>% dplyr::filter(biotype != "protein_coding"), aes(x = sample_ID, y = percentage, group = biotype, color = biotype)) +
    geom_line() + geom_point()
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/Biotypes/biotypes_percentages_noproteincoding.pdf"),
         width = 14, height = 8)
  
  #select protein coding genes
  counts_pc = counts_all %>% dplyr::filter(gene_type == "protein_coding")
  
  #RNA seq QC
  #create deseq object
  counts_std = counts_pc %>% dplyr::filter(gene_ID %in% geneInfo_s$gene_ID)
  rownames(counts_std) = counts_std$gene_ID
  counts_std = counts_std %>% dplyr::select(- colnames(geneInfo))
  system2("mkdir", "/home/rstudio/.cache")
  system2("mkdir", "/home/rstudio/.cache/R")
  system2("mkdir", "/home/rstudio/.cache/R/AnnotationHub")
  proxy = httr::use_proxy(Sys.getenv('http_proxy'))
  httr::set_config(proxy)
  AnnotationHub::setAnnotationHubOption("PROXY", proxy)
  dds = make_dds(counts = counts_std, metadata = meta, ah_record = ah_record)
  
  #library complexity
  plot_library_complexity(dds)
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/QC/library_complexity.pdf"),
         width = 14, height = 8)
  
  #gene detection
  plot_gene_detection(dds)
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/QC/gene_detection.pdf"),
         width = 14, height = 8)
  
  #gene filtering
  replicates = c()
  for(i in unique(meta$condition)){
    meta_i = meta %>% dplyr::filter(condition == i);
    replicates = append(replicates, length(rownames(meta_i)))
  }
  n = min(replicates)
  
  dds = filter_genes(dds, min_count = 3, min_rep = n)
  
  deseq_counts = dds@assays@data@listData[["counts"]]
  selected_IDs = c(rownames(deseq_counts), added_genes)
  geneInfo_f = geneInfo %>% dplyr::filter(gene_ID %in% selected_IDs)
  counts_f = counts_pc %>% dplyr::filter(gene_ID %in% geneInfo_f$gene_ID)
  rownames(counts_f) = counts_f$gene_ID
  counts_matrix_f = counts_f %>% dplyr::select(- colnames(geneInfo))
  
  rowdata = as.data.frame(dds@rowRanges@elementMetadata)
  rowdata$gene_ID = rownames(deseq_counts)
  rowdata = rowdata %>% dplyr::select(-gene_name)
  counts_f = counts_f %>% dplyr::left_join(rowdata, by = "gene_ID")
  write.csv(counts_f, file = paste0("/home/shared_folder/Output_", current_time, "/Tables/Counts/rawcounts_filtered.csv"))
  
  #Variance stabilization
  vsd = vst(dds)
  mean_sd_plot(vsd)
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/QC/variance_stabilization.pdf"),
         width = 14, height = 8)
  
  #Chromosomal expression
  for(i in chromosomes) {as.ggplot(plot_chromosome(vsd, as.character(i))); ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/Chromosomes/chr", as.character(i), ".pdf"),
                                                                                  width = 14, height = 7)}
  
  #Replicate variability
  plot_ma = plot_sample_MAs(vsd, group = "condition")
  for(i in 1:length(plot_ma)) {as.ggplot(plot_ma[[i]]); ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/Variability/", plot_ma[[i]][["labels"]][["title"]], ".pdf"),
                                                               width = 14, height = 7)}
  
  #clustering
  #Batch effect
  vsd2 = vsd
  assay(vsd2) = limma::removeBatchEffect(assay(vsd2), vsd2$replicate)
  
  #Clustering
  as.ggplot(plot_sample_clustering(vsd, n_feats = 1000, anno_vars = c("condition", "replicate"), distance = "pearson")) + ggtitle(" Pearson correlation before batch effect correction") + theme(plot.title = element_text(face = "bold"))
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/Correlation/clustering_beforecorrection.pdf"),
         width = 14, height = 8)
  as.ggplot(plot_sample_clustering(vsd2, n_feats = 1000, anno_vars = c("condition", "replicate"), distance = "pearson")) + ggtitle(" Pearson correlation after batch effect correction") + theme(plot.title = element_text(face = "bold"))
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/Correlation/clustering_aftercorrection.pdf"),
         width = 14, height = 8)
  
  #PCA
  #pca before batch correction
  plot_pca_scatters(vsd, n_PCs = PCA,  color_by = "condition", shape_by = "replicate")
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/PCA/PCA_1to", as.character(PCA), "_beforecorrection.pdf"),
         width = 14, height = 8)
  pca_res1 = plot_pca(vsd, show_plot = FALSE, n_feats = 1000)
  for(i in 1:PCA) {plot_loadings(pca_res1, PC = i, annotate_top_n = 10); ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/PCA/pca_top10_", as.character(i), "_beforecorrection.pdf"), width = 14, height = 8)}
  
  pca_var1 = as.data.frame(pca_res1[["var_exp"]])
  colnames(pca_var1) = "pca_var"
  pca_var1$pca = c(1:length(rownames(meta)))
  ggplot(pca_var1, aes(x = pca, y = pca_var)) + geom_point()
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/PCA/scree_plot_beforecorrection.pdf"),
         width = 14, height = 8)
  
  pca1 = as.data.frame(pca_res1[["data"]])
  ggplot(pca1, aes(x = PC1, y = PC2, shape = replicate)) + geom_point(aes(color = condition), size = 4) + scale_color_manual(values = pal) + ggtitle("PCA of the 1000 most variable genes before batch effect correction")
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/PCA/PCA_1_2_beforecorrection.pdf"),
         width = 14, height = 8)
  
  loadings1 = as.data.frame(pca_res1[["loadings"]])
  write.csv(loadings1, file = paste0("/home/shared_folder/Output_", current_time, "/Tables/PCA/loadings_beforecorrection.csv"))
  
  #pca after batch correction
  plot_pca_scatters(vsd2, n_PCs = PCA,  color_by = "condition", shape_by = "replicate")
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/PCA/PCA_1to", as.character(PCA), "_aftercorrection.pdf"),
         width = 14, height = 8)
  pca_res2 = plot_pca(vsd2, show_plot = FALSE, n_feats = 1000)
  for(i in 1:PCA) {plot_loadings(pca_res2, PC = i, annotate_top_n = 10); ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/PCA/pca_top10_", as.character(i), "_aftercorrection.pdf"), width = 14, height = 8)}
  
  pca_var2 = as.data.frame(pca_res2[["var_exp"]])
  colnames(pca_var2) = "pca_var"
  pca_var2$pca = c(1:length(rownames(meta)))
  ggplot(pca_var2, aes(x = pca, y = pca_var)) + geom_point()
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/PCA/scree_plot_aftercorrection.pdf"),
         width = 14, height = 8)
  
  pca2 = as.data.frame(pca_res2[["data"]])
  ggplot(pca2, aes(x = PC1, y = PC2, shape = replicate)) + geom_point(aes(color = condition), size = 4) + scale_color_manual(values = pal) + ggtitle("PCA of the 1000 most variable genes after batch effect correction")
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/QC/PCA/PCA_1_2_aftercorrection.pdf"),
         width = 14, height = 8)
  
  loadings2 = as.data.frame(pca_res2[["loadings"]])
  write.csv(loadings2, file = paste0("/home/shared_folder/Output_", current_time, "/Tables/PCA/loadings_aftercorrection.csv"))
  
  print("finishing QC")
  
  #boxplots
  print("starting boxplot creation")
  counts_norm = counts_matrix_f %>% dplyr::select(meta$sample_ID)
  counts_norm = as.matrix(counts_norm)
  counts_norm = as.data.frame(counts_norm/all_counts_sum$all_genes_sum)
  counts_norm_names = counts_norm
  counts_norm_names$gene_ID = rownames(counts_norm_names)
  counts_norm_names = counts_norm_names %>% dplyr::left_join(geneInfo_f, by = "gene_ID") %>% dplyr::select(-c("gene_ID", "gene_type"))
  counts_melt = melt(counts_norm_names)
  colnames(counts_melt) = c("gene_name", "sample_ID", "value")
  combs = as.list(as_tibble(as.data.frame(combn(unique(table$condition), 2))))
  
  if(length(boxplots) != 0) {
    for(i in boxplots){if(i %in% geneInfo$gene_name){melt_i = counts_melt %>% dplyr::filter(gene_name == i) %>% dplyr::left_join(meta, by = "sample_ID")
    ggplot(melt_i, aes(x = condition, y = value)) + geom_boxplot(aes(fill = condition)) + geom_point() + geom_signif(comparisons = combs, test = "wilcox.test", step_increase = 0.3) + scale_fill_manual(values = pal) + ggtitle(paste0("Boxplot ", i))
    ggsave(filename = paste0(paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Boxplots/boxplot_", i, ".pdf")), width = 14, height = 8)}}
  }
  
  print("finishing boxplot creation")
  
  #normalization
  #upload gene_length file
  gl = read_csv(gene_length)
  gl = as.data.frame(gl)
  colnames(gl) = c("gene_length")
  hgnc_symbol = c()
  transcript_length = c()
  for(i in gl$gene_length){row = strsplit(i, split = "\t"); hgnc_symbol = append(hgnc_symbol, row[[1]][1]); transcript_length = append(transcript_length, row[[1]][2])}
  geneLength = data.frame(hgnc_symbol, transcript_length)
  counts_tpm = as.data.frame(ADImpute::NormalizeTPM(as.matrix(counts_matrix_f), sce = NULL, log = FALSE, tr_length = geneLength, scale = 1))
  if(batch_correction == TRUE){counts_tpm = as.data.frame(limma::removeBatchEffect(counts_tpm, meta$replicate))}
  
  #heatmaps
  print("starting boxplot creation")
  
  ann_col = meta
  ann_col = as.data.frame(ann_col %>% dplyr::select(condition))
  rownames(ann_col) = meta$sample_ID
  
  as.ggplot(pheatmap(counts_tpm, scale = "row", show_rownames = FALSE, annotation_col = ann_col, color = viridis::viridis(n = 25, option = "H")))
  ggsave(filename = paste0(paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Heatmaps/heatmap_all_genes.pdf")),
         width = 14, height = 8)
  
  if(length(heatmaps) > 2){
    counts_hm = counts_tpm
    counts_hm$gene_ID = rownames(counts_hm)
    counts_hm = counts_hm %>% dplyr::left_join(geneInfo_f, by = "gene_ID") %>% dplyr::filter(gene_name %in% heatmaps)
    rownames(counts_hm) = counts_hm$gene_name
    counts_hm = as.matrix(counts_hm %>% dplyr::select(-colnames(geneInfo_f)))
    
    as.ggplot(pheatmap(counts_hm, scale = "row", show_rownames = TRUE, annotation_col = ann_col, color = viridis::viridis(n = 25, option = "H")))
    ggsave(filename = paste0(paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Heatmaps/heatmap_selected_genes.pdf")),
           width = 14, height = 8)
  }
  
  print("finishing heatmap creation")
  
  #WGCNA
  #select power
  print("starting moduls analysis")
  
  allowWGCNAThreads()
  powers = c(1:30)
  counts_w = t(counts_tpm)
  sft = pickSoftThreshold(counts_w, powerVector = powers, verbose = 5)
  indices = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
  for(i in 2:length(indices)){
    if(indices[i - 1] > 0){delta = abs(indices[i] - indices[i -1]); if(delta < 0.01){picked_power = i; break} else{picked_power = 30}
    } else{picked_power = 30}}
  indices = as.data.frame(indices)
  indices$n = powers
  picked_index = indices$indices[picked_power]
  ggplot(indices, aes(x = n, y = indices)) + geom_point() + geom_line() + geom_hline(aes(yintercept = picked_index), linetype = "dashed", color = "red") + geom_point(x = picked_power, y = picked_index, color = "red") + ggtitle("Scale indipendence") + xlab("Soft Threshold (power)") + ylab("Scale Free Topology Model Fit, signed R^2")
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Modules/scale_indipendence.pdf"),
         width = 14, height = 8)
  
  #create modules
  netwk = blockwiseModules(counts_w,                
                           power = picked_power, 
                           networkType = "signed",
                           deepSplit = 2,
                           pamRespectsDendro = F,
                           minModuleSize = 30,
                           maxBlockSize = 4000,
                           reassignThreshold = 0,
                           mergeCutHeight = 0.25,
                           saveTOMs = T,
                           saveTOMFileBase = "ER",
                           numericLabels = T,
                           verbose = 3)
  
  mergedColors = labels2colors(netwk$colors)
  
  module = data.frame(
    gene_ID = names(netwk$colors),
    colors = mergedColors)
  
  #plots of modules
  MEs0 = moduleEigengenes(counts_w, mergedColors)$eigengenes
  as.ggplot(pheatmap(MEs0, scale = "column", color = viridis::viridis(n = 25, option = "H"))) + ggtitle("Heatmap of modules")
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Modules/heatmap_modules.pdf"),
         width = 14, height = 8)
  
  module = module %>% dplyr::left_join(geneInfo_f, by = "gene_ID")
  
  cal_z_score = function(x){(x - mean(x)) / sd(x)}
  
  w_norm = as.data.frame(t(apply(as.data.frame(counts_w), 2, cal_z_score)))
  w_norm$sample_ID = rownames(w_norm)
  m_w = melt(w_norm)
  colnames(m_w) = c("gene_ID", "sample_ID", "value")
  m_w = m_w %>% dplyr::left_join(module, by = "gene_ID") %>% dplyr::left_join(meta, by = "sample_ID")
  ggplot(m_w, aes(x = sample_ID, y = value)) + geom_violin(aes(fill = condition, color = condition)) + facet_wrap(~colors) + scale_fill_manual(values = pal) + scale_color_manual(values = pal) + ggtitle("Violin plots of modules")
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Modules/violin_modules.pdf"),
         width = 14, height = 8)
  
  #covariance modules
  MEs0$sample_ID = rownames(MEs0)
  meta_modules = table %>% dplyr::select(c("sample_ID", "condition"))
  eigenvector_modules = MEs0 %>% dplyr::left_join(meta_modules, by = "sample_ID")
  melt = melt(eigenvector_modules)
  mean = melt %>% dplyr::group_by(condition, variable) %>% dplyr::summarise(mean = mean(value), sd = sd(value)) %>% dplyr::mutate(CV_rec = mean/sd)
  ggplot(mean) +
    geom_boxplot(aes(x = variable, y = mean)) +
    geom_point(aes(x = variable, y = mean, color = condition)) +
    geom_hline(aes(yintercept = 0), linetype = "dashed")
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Modules/boxplots_modules_mean.pdf"),
         width = 14, height = 8)
  
  ggplot(mean) +
    geom_boxplot(aes(x = variable, y = CV_rec)) +
    geom_point(aes(x = variable, y = CV_rec, color = condition)) +
    geom_hline(aes(yintercept = 0), linetype = "dashed")
  ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Modules/boxplots_modules_mean-sd.pdf"),
         width = 14, height = 8)
  
  colors = unique(module$colors)
  for(i in colors){
    genes_i = module %>% dplyr::filter(colors == i);
    c_i = counts_tpm;
    c_i$gene_ID = rownames(c_i);
    c_i = c_i %>% dplyr::filter(gene_ID %in% genes_i$gene_ID) %>% dplyr::select(-gene_ID);
    as.ggplot(pheatmap(c_i, scale = "row", annotation_col = ann_col, show_rownames = FALSE, color = viridis::viridis(n = 25, option = "H"))) + ggtitle(paste0("Heatmap of genes in ", as.character(i), " module"));
    ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Modules/Heatmaps/heatmap_", as.character(i), ".pdf"), width = 14, height = 8);
    
    GO = ensembldb::select(org, keys = genes_i$gene_name, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL");
    go = goana(GO$ENTREZID, species = specie, convert = T)
    res = topGO(go, ontology = c(ontology), number = GO_n, truncate.term = 50)
    ggplot(res, aes(x = reorder(Term, +DE), fill = - log10(res$P.DE), y = DE)) + geom_bar(stat = 'identity', width = 0.8) + coord_flip() +
      labs (title = "GO", colour = "-Log10Pval") + xlab("") + ylab("n of genes") + ggtitle(paste0("GO of ", as.character(ontology), " of genes in ", as.character(i), " module"));
    ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Modules/GO/GO_", as.character(ontology), "_", as.character(i), ".pdf"), width = 14, height = 14)
  }
  
  #save file
  counts_modules = counts_tpm
  counts_modules$gene_ID = rownames(counts_modules)
  counts_modules = counts_modules %>% dplyr::left_join(module, by = "gene_ID")
  write.csv(counts_modules, file = paste0("/home/shared_folder/Output_", current_time, "/Tables/Counts/TPMnormalized_counts_modules.csv"))
  write.csv(MEs0, file = paste0("/home/shared_folder/Output_", current_time, "/Tables/Counts/eigenvector_modules.csv"))
  
  print("finishing moduls analysis")
  
  #differential gene expression analysis
  print("starting differential gene expression analysis")
  
  cfr = c()
  for(i in colnames(table)){x =  grepl("cfr_", i); if(x == TRUE) {
    cfr = append(cfr, i);
    cfr_name = strsplit(i, "_")[[1]][2]
    meta_cfr = table %>% dplyr::select(c("sample_ID", "condition", "replicate", i));
    colnames(meta_cfr) = c("sample_ID", "condition", "replicate", "cfr");
    meta_cfr = meta_cfr %>% dplyr::filter(cfr != "no");
    counts_cfr = as.data.frame(counts_matrix_f) %>% dplyr::select(meta_cfr$sample_ID);
    comparison_names = meta_cfr$cfr;
    elements = unique(comparison_names);
    comp = c()
    for(j in comparison_names){
      if(j == elements[1]){comp = append(comp, "up")} else if(j == elements[2]){comp = append(comp, "down")}
    };
    design = model.matrix(~0+comp);
    contrasts = makeContrasts(cfr = compup - compdown, levels = colnames(design));
    v = voom(counts_cfr, design, plot = FALSE);
    if(batch_correction == TRUE){v[["E"]] = limma::removeBatchEffect(v[["E"]], meta_cfr$replicate)};
    l = lmFit(v, design);
    vfit = contrasts.fit(l, contrasts = contrasts);
    efit = eBayes(vfit);
    diff = topTreat(efit, coef = 1, n = Inf);
    diff$gene_ID = rownames(diff);
    diff = diff %>% dplyr::left_join(geneInfo_f, by = "gene_ID");
    col_to_plot = c();
    genes_to_plot = c();
    for(j in 1:length(rownames(diff))){
      if(diff$logFC[j] > logFC & diff$adj.P.Val[j] < pvalue){col_to_plot = append(col_to_plot, paste0("up_", elements[1])); genes_to_plot = append(genes_to_plot, diff$gene_name[j])}
      else if(diff$logFC[j] < - logFC & diff$adj.P.Val[j] < pvalue){col_to_plot = append(col_to_plot, paste0("down_", elements[2])); genes_to_plot = append(genes_to_plot, diff$gene_name[j])}
      else{col_to_plot = append(col_to_plot, "other"); genes_to_plot = append(genes_to_plot, "")}
    };
    diff$col_to_plot = col_to_plot;
    diff$genes_to_plot = genes_to_plot;
    if(length(unique(diff$col_to_plot)) == 3){
      ggplot(diff, aes(x = logFC, y = -log10(adj.P.Val))) + geom_point(aes(color = col_to_plot, size = 0.5, alpha = 0.5)) +
        scale_colour_manual(values = c("#3283FE", "grey", "#FA0087")) +
      ggrepel::geom_text_repel(aes(label = genes_to_plot), max.overlaps = 50,	size = 2, hjust = 1.2) + ggtitle(cfr_name);
      ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Volcano_plots/", cfr_name, ".pdf"), width = 14, height = 8)}
    if(length(unique(diff$col_to_plot)) == 1){
      ggplot(diff, aes(x = logFC, y = -log10(adj.P.Val))) + geom_point(aes(color = col_to_plot, size = 0.5, alpha = 0.5)) +
        scale_colour_manual(values = c("grey")) +
      ggrepel::geom_text_repel(aes(label = genes_to_plot), max.overlaps = 50,	size = 2, hjust = 1.2) + ggtitle(cfr_name);
      ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Volcano_plots/", cfr_name, ".pdf"), width = 14, height = 8)}
    if(length(unique(diff$col_to_plot)) == 2){if(sort(unique(diff$col_to_plot))[2] == "other"){
      ggplot(diff, aes(x = logFC, y = -log10(adj.P.Val))) + geom_point(aes(color = col_to_plot, size = 0.5, alpha = 0.5)) +
        scale_colour_manual(values = c("#3283FE", "grey")) +
      ggrepel::geom_text_repel(aes(label = genes_to_plot), max.overlaps = 50,	size = 2, hjust = 1.2) + ggtitle(cfr_name);
      ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Volcano_plots/", cfr_name, ".pdf"), width = 14, height = 8)}}
    if(length(unique(diff$col_to_plot)) == 2){if(sort(unique(diff$col_to_plot))[1] == "other"){
      ggplot(diff, aes(x = logFC, y = -log10(adj.P.Val))) + geom_point(aes(color = col_to_plot, size = 0.5, alpha = 0.5)) +
        scale_colour_manual(values = c("grey", "#FA0087")) +
      ggrepel::geom_text_repel(aes(label = genes_to_plot), max.overlaps = 50,	size = 2, hjust = 1.2) + ggtitle(cfr_name);
      ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Volcano_plots/", cfr_name, ".pdf"), width = 14, height = 8)}}
    write.csv(diff, file = paste0("/home/shared_folder/Output_", current_time, "/Tables/DGEs/", cfr_name, ".csv"));
    diff_sign = diff %>% dplyr::filter(col_to_plot != "_other");
    module_sign = module %>% dplyr::filter(gene_ID %in% diff_sign$gene_ID);
    module_sign = module_sign %>% arrange(colors);
    counts = counts_tpm;
    counts$gene_ID = rownames(counts);
    counts = counts %>% dplyr::filter(gene_ID %in% diff_sign$gene_ID) %>% dplyr::left_join(module_sign, by = "gene_ID") %>% dplyr::arrange(colors);
    rownames(counts) = counts$gene_ID;
    counts = counts %>% dplyr::select(meta$sample_ID);
    ann_row = module_sign;
    rownames(ann_row) = ann_row$gene_ID;
    ann_row = ann_row %>% dplyr::select(colors);
    colors = unique(module_sign$colors);
    names(colors) = colors;
    ann_colors = list(colors = colors);
    as.ggplot(pheatmap(as.matrix(counts), scale = "row", cluster_rows = TRUE, show_rownames = FALSE, annotation_row = ann_row, annotation_col = ann_col, annotation_legend = TRUE, annotation_names_col = TRUE, annotation_names_row = TRUE, annotation_colors = ann_colors, color = viridis::viridis(n = 25, option = "H"))) + ggtitle(paste0("Differentially expressed genes ", cfr_name));
    ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/Analysis/Heatmaps/", cfr_name, ".pdf"), width = 14, height = 8);
    genes = diff$logFC;
    names(genes) = diff$gene_ID;
    genes = sort(genes, decreasing = TRUE);
    gse = gseGO(geneList = genes, 
                ont = ontology, 
                keyType = "ENSEMBL",
                pvalueCutoff = 1, 
                OrgDb = org, 
                pAdjustMethod = "BH",
                by = "fgsea");
    res_gsea = gse@result;
    write.csv(res_gsea, file = paste0("/home/shared_folder/Output_", current_time, "/Tables/GSEA/", cfr_name, ".csv"));
    dotplot(gse, showCategory = GO_n, split = ".sign") + facet_grid(.~.sign);
    ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/GSEA/dotplots/", ontology, "_", cfr_name, ".pdf"), width = 14, height = 14);
    res_gsea$n = c(1:length(rownames(res_gsea)))
    for(j in GSEA){if(j %in% res_gsea$ID){
      res_j = res_gsea %>% dplyr::filter(ID == j);
      n = res_j$n;
      gseaplot(gse, by ="runningScore", color.line = "#009E73", title = paste0(gse$Description[n], " p = ", gse$p.adjust[n]), geneSetID = n) + theme(plot.title = element_text(face = "bold")) + scale_y_continuous(limits = c(-0.7, 0.7));
      ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/GSEA/enrichment_plots/runningScore_", cfr_name, "_", gse$Description[n], ".pdf"), width = 14, height = 8);
      core_genes = res_j$core_enrichment;
      core_genes = strsplit(core_genes, "/")[[1]];
      counts_gsea = counts_tpm;
      counts_gsea$gene_ID = rownames(counts_gsea);
      counts_gsea = counts_gsea %>% dplyr::filter(gene_ID %in% core_genes) %>% dplyr::select(-gene_ID);
      if(length(rownames(counts_gsea)) >= 2){as.ggplot(pheatmap(as.matrix(counts_gsea), scale = "row", cluster_rows = TRUE, show_rownames = TRUE, annotation_col = ann_col, annotation_legend = TRUE, annotation_names_col = TRUE, annotation_colors = ann_colors, color = viridis::viridis(n = 25, option = "H"))) + ggtitle(paste0("Enrichment core genes ", j, " ", cfr_name));
        ggsave(filename = paste0("/home/shared_folder/Output_", current_time, "/Plots/GSEA/enrichment_plots/core_", cfr_name, "_", gse$Description[n], ".pdf"), width = 14, height = 8)
      }}};
    diff_k = diff %>% dplyr::filter(adj.P.Val < pvalue);
    if(length(rownames(diff_k)) != 0){
      genes_k = diff_k$logFC;
      names(genes_k) = diff_k$gene_ID;
      ids = bitr(names(genes_k), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org);
      dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),];
      diff_k2 = diff_k[diff_k$gene_ID %in% dedup_ids$ENSEMBL,];
      diff_k2$K_id = dedup_ids$ENTREZID;
      kegg = diff_k2$logFC;
      names(kegg) = diff_k2$K_id;
      kegg = na.omit(kegg);
      kegg = sort(kegg, decreasing = TRUE);
      if(length(names(kegg)) > 2){
        kk = gseKEGG(geneList = kegg,
                     organism = kegg_organism,
                     pvalueCutoff = 1,
                     pAdjustMethod = "BH",
                     keyType = "kegg");
        resk = kk@result;
        write.csv(resk, file = paste0("/home/shared_folder/Output_", current_time, "/Tables/GSEA/Kegg_", cfr_name, ".csv"));
        for(j in pathways){if(j %in% resk$ID){
          pathview(gene.data = kegg, pathway.id = j, species = kegg_organism, low = list(gene = "#3283FE", cpd = "#3283FE"), mid = list(gene = "gray", cpd = "gray"), high = list(gene =  "#FA0087", cpd =  "#FA0087"), na.col = "transparent", kegg.dir = paste0("/home/shared_folder/Output_", current_time, "/Plots/GSEA/pathways/"));
          system2("mv", paste0(j, ".pathview.png /home/shared_folder/Output_", current_time, "/Plots/GSEA/pathways/", cfr_name, "_", j, ".png"))
        }}} else{kk = 0; resk = 0}} else{kk = 0; resk = 0}
    save(diff, counts, gse, res_gsea, kk, resk, file = paste0("/home/shared_folder/Output_", current_time, "/RData/data_", cfr_name, ".RData"))
  }}
  
  print("finishing differential gene expression analysis")
  

  #save Rfile
  print("saving R environment")
  
  save.image(paste0("/home/shared_folder/Output_", current_time, "/RData/environment.RData"))
  
}

print("finished succesfully")