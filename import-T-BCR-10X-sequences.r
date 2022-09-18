# Loading the necessary libraries.
library(Seurat)
library(tidyverse)

# Loading the Seurat object.
BCR <- readRDS("BCR-file.rds")

# Here, I seek to extract CDR1, CDR2, and CDR3, as well as all four FWR, amino acid and nucleotide sequences.

# First, we need the CTaa fields from the Seurat object ("BCR"), which correspond to the CDR3.
# This will be used to extract the other two CDR sequences, and FWR sequences for each chain.
# The Seurat object contains sequences combined into it using the scRepertoire package, to include (in this case) BCR data.


sequences <- BCR@meta.data %>%
  group_by(CTaa, sample, cell,  v_gene, j_gene, c_gene, lkv_gene, lkj_gene, lkc_gene) %>%
  as.data.frame() # Where "lk" means lambda/kappa chain, and other genes are for the heavy chains.
sequences <- mutate(sequences, cdr3 = CTaa, CTaa = NULL)
# Extracting the CDR3 region for each (heavy/light) chain for each cell.
sequences$hv3 <- vapply(str_split(sequences$cdr3, "[_]"), "[", "", 1)
sequences$lt3 <- vapply(str_split(sequences$cdr3, "[_]"), "[", "", 2)
# Earlier, using scRepertoire and Seurat, I pasted the sample name, which is the same as the folder containing the sample data, into each cell's barcode.
# Here I create a "folder" column to make it easier to handle each sample separately.
sequences$folder <- vapply(str_split(sequences$sample, "[_]"), "[", "", 1)

# Creating the data frame that will hold all the required sequences.
CDR123 <- data.frame()

# Looping through the cells in each folder to match cell barcodes in Seurat object to those in BCR filtere_contig_annotations.csv file for each sample.
for (i in levels(factor(sequences$folder))) {
  
    cdr <- read.table(paste0("from_cellranger/", i, "/vdj_b/filtered_contig_annotations.csv"), sep = ",", header = T)
    # The airr is needed to get the nucleotide sequences for whole BCR, but we used it to extract constant region (C gene) sequences.
    airr <- read.table(paste0("from_cellranger/", i, "/vdj_b/airr_rearrangement.tsv"), sep = "\t", header = T)
    
  
    for (j in 1:nrow(sequences[sequences$folder == i,])) {
      df <- data.frame(cell = sequences[sequences$folder == i,][j,]$cell,
                      sample = as.character(sequences[sequences$folder == i,][j,]$sample),
                      cdr3 = sequences[sequences$folder == i,][j,]$cdr3,
      # The following few lines should seek to match every cell barcode with a heavy and a light chain in the filtered contig file. Failing to find either, it will give the missing chaines NA values. "h" or "hv" denote heavy chains, and "l" or "lt" imply, you know, light chains; "nt" denotes nucleotide(s).
      hv1 = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$cdr1)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$cdr1),
      hv2 = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$cdr2)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$cdr2),
      hv3 = sequences[sequences$folder == i,][j,]$hv3,
      
      fwrh1 = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$fwr1)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$fwr1),
      fwrh2 = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$fwr2)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$fwr2),
      fwrh3 = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$fwr3)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$fwr3),
      fwrh4 = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$fwr4)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$fwr4),
      
      fwrl1 = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$fwr1)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$fwr1),
      fwrl2 = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$fwr2)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$fwr2),
      fwrl3 = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$fwr3)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$fwr3),
      fwrl4 = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$fwr4)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$fwr4),
      
      lt1 = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$cdr1)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$cdr1),
      lt2 = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$cdr2)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$cdr2),
      lt3 = sequences[sequences$folder == i,][j,]$lt3,
     
      hv1_nt = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$cdr1)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$cdr1_nt),
      hv2_nt = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$cdr2)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$cdr2_nt),
      hv3_nt = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$cdr2)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$cdr3_nt),
      
      fwrh1_nt = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$fwr1)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$fwr1_nt),
      fwrh2_nt = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$fwr2)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$fwr2_nt),
      fwrh3_nt = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$fwr3)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$fwr3_nt),
      fwrh4_nt = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$fwr4)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & cdr$chain == "IGH",]$fwr4_nt),
      
      fwrl1_nt = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$fwr1)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$fwr1_nt),
      fwrl2_nt = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$fwr2)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$fwr2_nt),
      fwrl3_nt = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$fwr3)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$fwr3_nt),
      fwrl4_nt = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$fwr4)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$fwr4_nt),
      
      lt1_nt = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$cdr1)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$cdr1_nt),
      lt2_nt = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$cdr2)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$cdr2_nt),
      lt3_nt = ifelse(length(cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$cdr2)==0,NA_character_,cdr[cdr$barcode == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & (cdr$chain == "IGL" | cdr$chain == "IGK"),]$cdr3_nt),
      
      HV = sequences[sequences$folder == i,][j,]$v_gene, HJ = sequences[sequences$folder == i,][j,]$j_gene, Ig = sequences[sequences$folder == i,][j,]$c_gene, 
      LV = sequences[sequences$folder == i,][j,]$lkv_gene, LJ = sequences[sequences$folder == i,][j,]$lkj_gene, LC = sequences[sequences$folder == i,][j,]$lkc_gene,
      
      # In the airr-formatted file, we can find the whole sequence, with each gene/region having a start location and an end location. I use these to exctract constant region sequences, but the same approach can be used to extract any other regions of interest. The function substr() extracts the sequences by using the start and end locations for each chain. It's worthy of note here that this file shows sequences start and end locations for nucleotides only, not amino acids.
      constant_lt_nt = ifelse(length(airr[airr$cell_id == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & grepl(airr$v_call, pattern = "(L)|(K)"),]$sequence_aa)==0,NA_character_, substr(x = airr[airr$cell_id == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & grepl(airr$v_call, pattern = "(L)|(K)"),]$sequence, 
                                                                                                                                                                                                                           start = airr[airr$cell_id == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & grepl(airr$v_call, pattern = "(L)|(K)"),]$c_sequence_start,
                                                                                                                                                                                                                           stop = airr[airr$cell_id == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & grepl(airr$v_call, pattern = "(L)|(K)"),]$c_sequence_end)),
      
     
      constant_hv_nt = ifelse(length(airr[airr$cell_id == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & grepl(airr$v_call, pattern = "(H)"),]$sequence_aa)==0,NA_character_, substr(x = airr[airr$cell_id == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & grepl(airr$v_call, pattern = "(H)"),]$sequence, 
                                                                                                                                                                                                                           start = airr[airr$cell_id == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & grepl(airr$v_call, pattern = "(H)"),]$c_sequence_start,
                                                                                                                                                                                                                           stop = airr[airr$cell_id == vapply(strsplit(sequences[sequences$folder == i,][j,]$cell, "[_]"), "[", "", 3) & grepl(airr$v_call, pattern = "(H)"),]$c_sequence_end))
      
      
      
      
    )
     # Combining the dataframes for the individual samples.
     CDR123 <- rbind(CDR123, df)
    }
       

  }
  
# Saving the data frame as a table, including all V(D)J genes and all fwrs of both chains, and constant regions, and also combining duplicate chains, giving us an indication of clonotype abundance/expansion for each sample.
CDRabundance <- CDR123 %>%
  group_by(sample, fwrh1, hv1, fwrh2, hv2, fwrh3, hv3, fwrh4, cdr3, fwrl1, lt1, fwrl2, lt2, fwrl3, lt3, fwrl4,
           fwrh1_nt, hv1_nt, fwrh2_nt, hv2_nt, fwrh3_nt, hv3_nt, fwrh4_nt, constant_hv_nt, 
           fwrl1_nt, lt1_nt, fwrl2_nt, lt2_nt, fwrl3_nt, lt3_nt, fwrl4_nt, constant_lt_nt,
           HV, HJ, Ig, LV, LJ, LC) %>% 
  count() %>% arrange(desc(n)) %>% 
  mutate(full = paste0(fwrh1, hv1, fwrh2, hv2, fwrh3, hv3, fwrh4, fwrl1, lt1, fwrl2, lt2, fwrl3, lt3, fwrl4),
         full.length = sapply(full, nchar),
         full_nt = paste0(fwrh1_nt, hv1_nt, fwrh2_nt, hv2_nt, fwrh3_nt, hv3_nt, fwrh4_nt, constant_hv_nt, fwrl1_nt, lt1_nt, fwrl2_nt, lt2_nt, fwrl3_nt, lt3_nt, fwrl4_nt, constant_lt_nt),
         full_nt.length = sapply(full_nt, nchar),
         ) %>%
  relocate(n, .after = full_nt.length)
write.table(CDRabundance, "sequences-by-sample.tsv", sep = "\t", col.names = NA)

