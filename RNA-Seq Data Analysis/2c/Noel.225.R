# set working directory to final project folder
setwd("/fs/ess/PAS3124/noel.225/MG5795-2025/Final5795Project/RNA-Seq Data Analysis/2c")

# creating tables of the data using the .tab files
c1 <- read.table("c1_zm_ReadsPerGene.out.tab",
                 sep = "\t",
                 skip = 4,
                 col.names = c("gene_id", "unstranded", "forward", "reverse"))
c2 <- read.table("c2_zm_ReadsPerGene.out.tab",
                 sep = "\t",
                 skip = 4,
                 col.names = c("gene_id", "unstranded", "forward", "reverse"))
m1 <- read.table("m1_zm_ReadsPerGene.out.tab",
                 sep = "\t",
                 skip = 4,
                 col.names = c("gene_id", "unstranded", "forward", "reverse"))
m2 <- read.table("m2_zm_ReadsPerGene.out.tab",
                 sep = "\t",
                 skip = 4,
                 col.names = c("gene_id", "unstranded", "forward", "reverse"))

# check to ensure gene order matches in the tables
stopifnot(all(c1$gene_id == c2$gene_id))
stopifnot(all(c2$gene_id == m1$gene_id))
stopifnot(all(m1$gene_id == m2$gene_id))
stopifnot(all(m2$gene_id == c1$gene_id))

# creating the count table using the unstranded counts
count_table <- data.frame(
  c1 = c1$unstranded,
  c2 = c2$unstranded,
  m1 = m1$unstranded,
  m2 = m2$unstranded)

# making gene_id the row names
rownames(count_table) <- c1$gene_id

# making the table into a file 
write.csv(count_table, "/fs/ess/PAS3124/noel.225/MG5795-2025/Final5795Project/RNA-Seq Data Analysis/2c/unstranded_count_table_4maize.csv")

# packages required to compare datasets and visualize them
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# identify genes with sufficient counts to analyse
count_at_least_10 = count_table >= 10
num_samples_with_enough_counts_each_gene = rowSums(count_at_least_10)
gene_passes_filter = num_samples_with_enough_counts_each_gene >= 3
count_table_filtered = count_table[gene_passes_filter, ]

# compare control replicates - scatter plot 
ggplot(count_table_filtered, aes(x = log10(c1+1), y = log10(c2+1))) + # map the x and y variables to ggplot aesthetics
  geom_point(color = "blue", alpha = 0.6) +
  labs(
    title = "Maize Control 1 vs 2",
    x = "Control 1",
    y = "Control 2"
  )

# compare LARP4B kd replicates - scatter plot 
ggplot(count_table_filtered, aes(x = log10(m1+1), y = log10(m2+1))) + # map the x and y variables to ggplot aesthetics
  geom_point(color = "red", alpha = 0.6) +
  labs(
    title = "Maize Mutant 1 vs 2",
    x = "Mutant 1",
    y = "Mutant 2"
  )

# DESeq2 analysis
# gathering necessary packages
library("DESeq2")
library("biomaRt")

# read in count table
count_table = read.csv("unstranded_count_table_4maize.csv",
                       row.names = 1)

#identify genes with sufficient counts to analyse
count_at_least_10 = count_table >= 10
num_samples_with_enough_counts_each_gene = rowSums(count_at_least_10)
gene_passes_filter = num_samples_with_enough_counts_each_gene >= 3
count_table_filtered = count_table[gene_passes_filter, ]

#define metadata
metadata = data.frame("condition" = c("ctrl", "ctrl", "mut", "mut"))
rownames(metadata) = colnames(count_table_filtered)

# turn condition column of metadata into a factor
metadata$condition = factor(metadata$condition)

#create DESeqDataSet Object
DE_dataset = DESeqDataSetFromMatrix(countData = count_table_filtered,
                                    colData = metadata,
                                    design = ~ condition)

# Make a PCA plot
# first do variance stabilizing transformation
vst = varianceStabilizingTransformation(DE_dataset)

# Now plot PCA
plotPCA(vst)

# Set control knockdown as reference level
DE_dataset$condition <- relevel(DE_dataset$condition, "ctrl")

# run the differential expression analysis
DE_dataset = DESeq(DE_dataset)

# plot MA plot
plotMA(DE_dataset,
       alpha = 0.1,
       main = "",
       xlab = "mean of normalized counts",
       colNonSig = "gray60",
       colSig = "blue",
       colLine = "grey40",
       returnData = FALSE,
       MLE = FALSE
)

# get the results for control vs mutant comparison
condition_results = results(DE_dataset, 
                            contrast = c("condition","mut", "ctrl"),
                            independentFiltering = FALSE)

# connect to Ensembl Plants, maize dataset
ensembl_zm <- useMart(
  biomart = "plants_mart",
  dataset = "zmays_eg_gene",
  host = "plants.ensembl.org"
)

# remove version suffix from gene IDs if present (safe even if none)
rownames(condition_results) <- gsub("\\.[0-9]{1,3}$", "", rownames(condition_results))

# get a mapping table from Ensembl gene ID -> external gene name (symbol-like)
gene_id_table <- getBM(
  mart = ensembl_zm,
  attributes = c("ensembl_gene_id", "external_gene_name")
)

# remove duplicates (rare)
gene_id_table <- gene_id_table[!duplicated(gene_id_table$ensembl_gene_id), ]

# merge names onto DESeq2 results
condition_results_with_names <- merge(
  gene_id_table,
  as.data.frame(condition_results),
  all.y = TRUE,
  by.x = "ensembl_gene_id",
  by.y = "row.names"
)

# save results
write.csv(condition_results_with_names,
          "maize_mut_vs_ctrl_DESeq2_results_with_names.csv",
          row.names = FALSE)

# Analysis of DESeq2 results
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

Maize_DEseq <- read_csv("maize_mut_vs_ctrl_DESeq2_results_with_names.csv")

glimpse(Maize_DEseq)

# Select certain columns and get summary
log2fc <- dplyr::select(Maize_DEseq, log2FoldChange) 
summary(log2fc)

# Create a new column for 2-fold significance
Maize_DEseq$sig_2f <- with(Maize_DEseq, ifelse(log2FoldChange > 1 & padj < 0.05, "Up",
                                                 ifelse(log2FoldChange < -1 & padj < 0.05, "Down", "Not significant")))

# Counts significant genes
table(Maize_DEseq$sig_2f)

# Volcano Plot 1 (2-fold change)
ggplot(Maize_DEseq, aes(x = log2FoldChange, y = -log10(padj), color = sig_2f)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not significant" = "black")) +
  labs(title = "Mutant vs Control KD (2-fold change)", x = "Log2 Fold Change", y = "-log10(Adjusted P-value)") +
  theme_minimal()

# Next we will make a list of genes that go significantly up and down. But first remove rows where hgnc_symbol is NA
Maize_DEseq <- subset(Maize_DEseq, !is.na(external_gene_name) & external_gene_name != "")

# Subset for Up-regulated genes
up_genes <- subset(Maize_DEseq, sig_2f == "Up")

# Subset for Down-regulated genes
down_genes <- subset(Maize_DEseq, sig_2f == "Down")

# Write to CSV (only hgnc_symbol column)
write.csv(up_genes$external_gene_name, "up_genes.csv", row.names = FALSE)
write.csv(down_genes$external_gene_name, "down_genes.csv", row.names = FALSE)
write.csv(Maize_DEseq$external_gene_name, "all_genes.csv", row.names = FALSE)