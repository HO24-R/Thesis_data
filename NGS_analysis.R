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


#Heatmap using pheatmap

setwd("/Volumes/begleylab/Humphrey/PhD_project/NGS_2/30-990495916/00_fastq/fastq_files/STARoutput_NGS2_2/DE_files/Stationary_phase/dmnmA_vs_WT_analysis")
files <- list.files(path = "/Volumes/begleylab/Humphrey/PhD_project/NGS_2/30-990495916/00_fastq/fastq_files/STARoutput_NGS2_2/DE_files/Stationary_phase/dmnmA_vs_WT_analysis", full.names = FALSE)
print(files)

# Step 1: Load Data
DEGs <- read.csv("Complete_DEGs_sta_to_WT_untreated.csv")

# Step 2: Rename Columns for Clarity
colnames(DEGs) <- c(
  "GeneName",
  "∆mnmA tr vs WT untr",  # Replace dmnmA_tr_vs_WT_unt
  "∆mnmA untr vs WT untr", # Replace dmnmA_unt_vs_WT_unt
  "WT tr vs WT untr"       # Replace WT_tr_vs_WT_unt
)

print(colnames(DEGs))  # Check before conversion to matrix

# Step 3: Reshape to Long Format to Handle Duplicates
library(tidyverse)
DEGs_long <- DEGs %>%
  pivot_longer(
    cols = -GeneName,            # Exclude GeneName from reshaping
    names_to = "Condition",      # New column for condition names
    values_to = "Log2FC"         # New column for Log2FC values
  )

print(colnames(DEGs_long))  # Check before conversion to matrix

# Step 4: Aggregate Data to Handle Duplicate GeneNames
DEGs_aggregated <- DEGs_long %>%
  group_by(GeneName, Condition) %>%        # Group by GeneName and Condition
  summarise(
    Log2FC = mean(Log2FC, na.rm = TRUE),   # Calculate mean Log2FC for duplicates
    .groups = 'drop'
  )

print(colnames(DEGs_aggregated))  # Check before conversion to matrix

# Step 5: Pivot Back to Wide Format for Heatmap Clustering
DEGs_wide <- DEGs_aggregated %>%
  pivot_wider(
    names_from = Condition,      # Use Condition as column names
    values_from = Log2FC         # Use Log2FC values to fill cells
  )

print(colnames(DEGs_wide))  # Check before conversion to matrix

# Step 6: Handle Missing or Invalid Values

# First, ensure that the dataframe is in the correct format
# Check if GeneName is indeed a column
if ("GeneName" %in% colnames(DEGs_wide)) {
  
  # Ensure GeneName is not NA or empty
  DEGs_wide <- DEGs_wide[!is.na(DEGs_wide$GeneName) & DEGs_wide$GeneName != "", ]
  
  # Replace NaN, NA, and Inf in numeric columns only
  numeric_cols <- sapply(DEGs_wide, is.numeric)  # Identify numeric columns
  DEGs_wide[, numeric_cols] <- lapply(DEGs_wide[, numeric_cols], function(x) {
    x[is.na(x) | is.nan(x) | is.infinite(x)] <- 0  # Replace invalid values with 0
    return(x)
  })
  
} else {
  stop("GeneName column is missing from DEGs_wide.")
}

print(colnames(DEGs_wide))  # Check before conversion to matrix

# Ensure unique GeneName by appending suffixes if necessary
DEGs_wide$GeneName <- make.unique(as.character(DEGs_wide$GeneName))

print(colnames(DEGs_wide))  # Check before conversion to matrix

# Convert to matrix
DEGs_matrix <- DEGs_wide %>%
  column_to_rownames(var = "GeneName") %>%
  as.matrix()

print(colnames(DEGs_matrix))  # Check before conversion to matrix


### if scaling is not needed, proceed to step 8; otherwise continue with step 7 for scaling

# Step 7 Apply scaling to each row (not column)
DEGs_matrix_scaled <- t(apply(DEGs_matrix, 1, function(x) {
  scaled <- scale(x)  # Standardize to mean = 0 and SD = 1
  scaled[scaled > 1.5] <- 1.5  # Cap max at 1.5
  scaled[scaled < -1.5] <- -1.5  # Cap min at -1.5
  return(scaled)
}))

# Restore row and column names
rownames(DEGs_matrix_scaled) <- rownames(DEGs_matrix)
colnames(DEGs_matrix_scaled) <- colnames(DEGs_matrix)  # Ensure column names are retained

# Check column names again
print(colnames(DEGs_matrix_scaled))  # Should now correctly show column names


DEGs_matrix <- t(DEGs_matrix)  # Transpose back to original dimensions


# Step 8: Create Clustered Heatmap with Labeled Columns
library(pheatmap)
pheatmap(
  DEGs_matrix,
  scale = "none",               # No additional scaling; already normalized to -1.5 to +1.5
  clustering_distance_rows = "euclidean",  # Distance metric for rows
  clustering_distance_cols = "euclidean",  # Distance metric for columns
  clustering_method = "complete",         # Hierarchical clustering method
  color = colorRampPalette(c("green", "white", "red"))(50),  # Custom color gradient
  breaks = seq(-1.5, 1.5, length.out = 51),  # Force color scale to span -2.0 to +2.0
  legend_breaks = seq(-1.5, 1.5, by = 0.5),  # Set legend ticks at internvals of 0.5
  main = "Clustered Heatmap of mnmA WT all_genes_scaled",
  fontsize_row = 8,             # Font size for row labels
  fontsize_col = 10,            # Font size for column labels
  labels_col = colnames(DEGs_matrix)  # Add descriptive column names
)
