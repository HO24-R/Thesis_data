# DE Seq analysis
rm(list = ls())
setwd("add the data path here")
getwd() # to confirm your path is valid

library(DESeq2)

# Loading the .out.tab files into your variable
file.list <- list.files( path = "./", pattern = "*ReadsPerGene.out.tab$")

# Forming a combined counts table
counts.files <- lapply(file.list, read.table, skip = 4)

# Selecting counts column
counts <- as.data.frame( sapply( counts.files, function(x) x[ ,2]))

# Setting columns and rows for the data frame
colnames(counts) <- file.list
rownames(counts) <- counts.files[[1]]$V1

# See what I got
print(colnames(counts))
print(rownames(counts))

# Mine was 4 biological replicates; adjust according to your number of biological replicates
condition <- c(rep("WT_Treated",4), rep("WT_Untreated",4))

print(condition) # to confirm that your conditions were captured effectively

# Forming a meta data matrix
sampleTable <- data.frame(sampleName = file.list, condition = condition)

print(sampleTable) # to confirm that your conditions were captured effectively

# DESeq dataset formation
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sampleTable, design = ~ condition)

# Running DESeq
output <- DESeq(dds)

# Obtaining results
results_WT_tr_vs_WT_unt <- results(output, contrast=c("condition","WT_Treated","WT_Untreated"))

# Sorting results by Adjusted P-value
results_WT_tr_vs_WT_unt_PValue <- results_WT_tr_vs_WT_unt[order(results_WT_tr_vs_WT_unt$padj),]

# Heading of results and write to file
head(results_WT_tr_vs_WT_unt_PValue)

write.csv(as.data.frame(results_WT_tr_vs_WT_unt_PValue), file="WT_exp_tr_vs_WT_exp_unt_DE.csv") # To save the output; change title as you deem fit

# End


# Bacterial gene IDs do not always show like eukaryotes where we can map it using ENSEMBL ID, so I used this process

# Mapping of gene names to gene ID in DESeq files
rm(list = ls())

library(ggplot2)
library(dplyr)


DExp <- read.csv("WT_exp_tr_vs_WT_exp_unt_DE.csv")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rtracklayer")

library(rtracklayer)


gff_file <- "/Users/homeoga/Documents/NGS_data/Annot.gff3" # Download the "Annot.gff3" file and give the path

gff_data <- import.gff3(gff_file, format = "GFF3")

View(gff_data)
gene_ids <- DExp$gene


# Extracting Gene ID and Gene Name Info from gff annotation file
gene_IDs <- mcols(gff_data)$gene_id
gene_names <- mcols(gff_data)$Name

gene_mapping <- data.frame(genes = gene_IDs, GeneName = gene_names, stringsAsFactors = FALSE)
merged_data <- left_join(DExp, gene_mapping, by = "genes")

write.csv(merged_data, "WT_exp_tr_vs_WT_exp_unt_DE_with_Gene_Names.csv", row.names=FALSE)