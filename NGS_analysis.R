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

#load l=ibraries

library(ggplot2)
library(dplyr)
library(rtracklayer)


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


# Getting upregulated and downregulated genes that do align with volcano plot number of genes

rm (list = ls())
setwd("path_to_directory_or_folder_of_interest")

# Getting upregulated and downregulated genes to align with volcano plot

DE_csv <- "name_of_file.csv"

dmnmH_tr_vs_dmnmH_untr <- read.csv(DE_csv)

#Modified to align with volcano plot

# Define thresholds
fold_change_threshold <- 1.5
significance_threshold <- 0.05

# Add classification column using pvalue
dmnmH_tr_vs_dmnmH_untr$diffexpressed <- with(dmnmH_tr_vs_dmnmH_untr,
                                             ifelse(log2FoldChange < -fold_change_threshold & pvalue < significance_threshold, "Downregulated",
                                                    ifelse(log2FoldChange > fold_change_threshold & pvalue < significance_threshold, "Upregulated", "NS"))
)


# Filter genes
upregulated_genes <- subset(dmnmH_tr_vs_dmnmH_untr, diffexpressed == "Upregulated")
downregulated_genes <- subset(dmnmH_tr_vs_dmnmH_untr, diffexpressed == "Downregulated")

# Save results
write.csv(upregulated_genes, "Upregulated_dmnmH_treated_vs_dmnmH_untreated.csv", row.names=FALSE)
write.csv(downregulated_genes, "Downregulated_dmnmH_treated_vs_dmnmH_untreated.csv", row.names=FALSE)


################### VOLCANO PLOT INCREASE THE X-AXIS BY 1 ###################################
rm(list = ls())

setwd("path_to_directory_or_folder_of_interest")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
install.packages("installr")
library(installr)
updateR()

library(readr)

DExp <- read_csv("name_of_file.csv")

# Assuming DExp is your original data frame

# Add -log10(padj) column
DExp$log_padj <- -log10(DExp$padj)

# Define thresholds
fold_change_threshold <- 1.5
significance_threshold <- 0.05

# add a column of NAs
DExp$diffexpressed <- "NS"
# if log2Foldchange < 1.5, >-1.5
DExp$diffexpressed[DExp$log2FoldChange < 1.5 & DExp$log2FoldChange > -1.5 & DExp$pvalue > 0.05] <- "NS"
# if log2Foldchange > 1.5 and pvalue > 0.05, set as "log2FoldChange" 
DExp$diffexpressed[DExp$log2FoldChange > 1.5 & DExp$pvalue > 0.05] <- "log2FoldChange"
# if log2Foldchange < -1.5 and pvalue > 0.05, set as "log2FoldChange" 
DExp$diffexpressed[DExp$log2FoldChange < -1.5 & DExp$pvalue > 0.05] <- "log2FoldChange"
# if log2Foldchange > -1.5, < 1.5, and pvalue < 0.05, set as "pvalue"
DExp$diffexpressed[DExp$log2FoldChange > -1.5 & DExp$log2FoldChange < 1.5 & DExp$pvalue < 0.05] <- "pvalue"
# if log2Foldchange > -1.5, < 1.5, and pvalue < 0.05, set as "p-value and log2 FC"
DExp$diffexpressed[DExp$log2FoldChange < -1.5 & DExp$pvalue < 0.05] <- "p-value and log2 FC"
DExp$diffexpressed[DExp$log2FoldChange > 1.5 & DExp$pvalue < 0.05] <- "p-value and log2 FC"


# Create volcano plot with x-axis increment of 1
p <- ggplot(data = DExp, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed)) +
  geom_point(size = 2.5) +
  theme_minimal() +
  scale_color_manual(values = c("green", "grey", "red", "blue")) +
  geom_vline(xintercept = c(-fold_change_threshold, fold_change_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(significance_threshold), linetype = "dashed", color = "black") +
  scale_x_continuous(
    breaks = seq(floor(min(DExp$log2FoldChange, na.rm = TRUE)), 
                 ceiling(max(DExp$log2FoldChange, na.rm = TRUE)), 
                 by = 1)  # Increment by 1
  ) +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-log10(p-value)"
  ) +
  theme(
    axis.text = element_text(size = 10)
  )

# Print the plot
print(p)

############ PCA PLOT WITH ELLIPSE ##############################

rm(list = ls())
setwd("path_to_directory_or_folder_of_interest")

BiocManager::install("DESeq2") #if not previously installed
library(DESeq2)
library(ggplot2)
library(ggforce)

path <- "/Volumes/begleylab/Humphrey/PhD_project/NGS_2/30-990495916/00_fastq/fastq_files/STARoutput_NGS2_2/DE_files/PCA/dmnmA"
file_names <- list.files(path, pattern="ReadsPerGene.out.tab$", full.names = TRUE)
conditions <- rep(c("dmnmA_exp_treated", "dmnmA_exp_untreated", "dmnmA_sta_treated", "dmnmA_sta_untreated"), each = 4)

# Function to read counts and gene names from each file
read_data <- function(file) {
  dat <- read.delim(file, header=FALSE, stringsAsFactors=FALSE)
  gene_names <- dat[-c(1:4), 1]  # Assuming gene IDs start from the 5th row
  counts <- as.numeric(dat[-c(1:4), 2])  # Assuming counts are in the second column
  list(counts = counts, gene_names = gene_names)
}

# Read data from all files, extract counts and gene names
all_data <- lapply(file_names, read_data)
count_matrix <- do.call(cbind, lapply(all_data, `[[`, "counts"))
gene_names <- all_data[[1]]$gene_names  # Assuming all files have the same gene order

# Set gene names as row names in the count matrix
rownames(count_matrix) <- gene_names
colnames(count_matrix) <- conditions

# Create a DESeq2 dataset
col_data <- DataFrame(condition = factor(conditions))
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ condition)

# Normalize data and perform variance stabilizing transformation
dds <- DESeq(dds)
vst_data <- vst(dds, blind = FALSE)

# Assuming 'vst_data' is the variance-stabilized data from DESeq2
# Ensure to perform PCA on the transposed data to have samples as observations
pca_results <- prcomp(t(assay(vst_data)))

# Create a data frame for PCA results
# Ensure 'conditions' matches the number of principal component scores
pca_data <- data.frame(
  PC1 = pca_results$x[, 1], 
  PC2 = pca_results$x[, 2], 
  Condition = conditions  # Assuming 'conditions' has a length equal to the number of columns in count_matrix
)
# Check dimensions to ensure they match
print(nrow(pca_results$x))  # Should print 16
print(length(conditions))   # Should also print 16

if (!require("ggforce")) install.packages("ggforce")
library(ggforce)


# Create a PCA plot with dotted ellipses around each cluster for 4 groups
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) +  # Plot points
  geom_mark_ellipse(aes(fill = Condition, group = Condition), 
                    colour = "black", linetype = "dotted", 
                    show.legend = FALSE, alpha = 0.3) +  # Add dotted ellipses
  theme_minimal() +
  labs(title = "PCA of RNA-seq Data",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  scale_fill_manual(values = c("red", "white", "yellow", "purple"))  # Customize colors if needed


## for 6 groups
# Use coord_cartesian for adjusting view without clipping data
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) +
  geom_mark_ellipse(aes(fill = Condition, group = Condition), 
                    colour = "black", linetype = "dotted", 
                    show.legend = FALSE, alpha = 0.3) +
  theme_minimal() +
  coord_cartesian(xlim = c(min(pca_data$PC1) - 5, max(pca_data$PC1) + 5),
                  ylim = c(min(pca_data$PC2) - 5, max(pca_data$PC2) + 5)) +
  labs(title = "PCA of RNA-seq Data",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  scale_fill_manual(values = c("red", "blue", "green", "yellow", "white", "purple"))  # Customize colors if needed



#Heatmap using pheatmap

# Load necessary libraries
library(tidyr)
library(dplyr)
library(pheatmap)
library(tidyverse)

setwd("path_to_directory_or_folder_of_interest")

# Step 1: Load Data
DEGs <- read.csv("name_of_file.csv")

# Step 2: Rename Columns for Clarity
colnames(DEGs) <- c(
  "GeneName",
  "∆mnmA tr vs WT untr",  # Replace dmnmA_tr_vs_WT_unt
  "∆mnmA untr vs WT untr", # Replace dmnmA_unt_vs_WT_unt
  "WT tr vs WT untr"       # Replace WT_tr_vs_WT_unt
)

print(colnames(DEGs))  # Check before conversion to matrix

# Step 3: Reshape to Long Format to Handle Duplicates using tidyverse
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


# Step 8: Create Clustered Heatmap with Labeled Columns using pheatmap
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


# To identify a group of genes within the larger set of differentially expressed genes

############################# Collating genes from a Master List ###########################################

setwd("")

# Load the necessary library for working with CSV files
library(readr)

# Load the master list CSV file into a data frame
master_list <- read.csv("list.csv") #change name according to yours


# rpoS genes
# List of genes of interest

genes_of_interest <- c("yahO", "bolA", "ybaY", "ybaS", "ybaT", "ybgS", 
                       "ybhE", "dps", "ycaC", "cbpA", "yccJ", "wrbA", "msyB", "ycfH", "ycgZ", 
                       "ymgA", "ymgB", "ymgC", "ycgB", "ymgE", "treA", "osmB", "ydcS", "ydcV", 
                       "narU", "sra", "osmC", "xasA", "gadB", "ydeI", "sodC", "ynhG", "sufS", 
                       "sufD", "sufB", "sufA", "ydiZ", "osmE", "yeaG", "yeaR", "otsA", "otsB", 
                       "yodD", "yodC", "hchA", "yegP", "fbaB", "yohC", "elaB", "talA", "tktB", 
                       "ygaU", "ygiW", "yqjC", "yqjD", "yqjE", "yhcO", "bfr", "fic", "yhfG", 
                       "slp", "hdeB", "hdeA", "hdeD", "yhiE", "gadW", "gadX", "gadA", "yiaG", 
                       "yjbJ", "yjdI", "yjdJ", "aidB", "ytfQ", "osmY", "yahO", "xthA", "gadC", 
                       "crp", "nmpC", "yliH", "ompF", "serT", "flgM", "flgA", "flgB", "flgC", 
                       "flgD", "flgE", "flgF", "flgG", "flgH", "flgI", "ompW", "ynaJ", "uspE", 
                       "fnr", "abgT", "abgB", "abgR", "ydaM", "ydaN", "ydaO", "intR", "ydaC", 
                       "lar", "recT", "racC", "kil", "ydaG", "ydaR", "ynaA", "insH5", "ynfE", 
                       "ynfF", "fliZ", "fliA", "fliD", "fliM", "fliN", "yfiD", "ansB", "tnaA", 
                       "glyY", "pyrB", "fimA", "yjbJ", "yjbE", "aldB", "gabP", "ygaU", "ygdI", 
                       "yqhE", "mltB", "narY", "ygaF", "ybaY", "yjgR", "ydaM", "ydcS", "ydcK", 
                       "gabD”, “gabT", "yhjY", "yhiN", "yeaG", "phnP", "yebF", "yodC", "gabD", 
                       "talA", "aroM", "ugpE", "nlpA”, “yicM", "argH”,   ybiO", "yhjG", "argH", 
                       "yliI", "ugpC", "yhiV", "uspB", "ldcC", "appB", "yhiU", "aidB", "yehX", 
                       "yphA", "osmY", "yfcG", "otsA", "ecnB", "yhjD", "yjiN", "ybiO", "katE", "oxyR")

# Filter the master list to get details of the genes of interest

genes_of_interest_data <- master_list[master_list$Gene %in% genes_of_interest, ]

# Print the result
print(genes_of_interest_data)

write.csv(genes_of_interest_data, "rpoS_regulated_genes_dmnmH_exp_tr_vs_dmnmH_exp_unt.csv")
