# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

# Read the differentially expressed genes (Lung and Blood)
de_genes_lung <- read.csv("significant_genes_lung.csv")$ENSEMBL
de_genes_blood <- read.csv("significant_genes_blood.csv")$ENSEMBL

# Get the list of all genes (Ensembl IDs)
all_genes <- keys(org.Mm.eg.db, keytype = "ENSEMBL")

# Perform GO enrichment analysis for Lung
ego_lung <- enrichGO(
  gene = de_genes_lung,
  universe = all_genes,
  OrgDb = org.Mm.eg.db,
  ont = "BP",         # You can change this to "MF", "CC", or "ALL"
  keyType = "ENSEMBL",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Perform GO enrichment analysis for Blood
ego_blood <- enrichGO(
  gene = de_genes_blood,
  universe = all_genes,
  OrgDb = org.Mm.eg.db,
  ont = "BP",         # You can change this to "MF", "CC", or "ALL"
  keyType = "ENSEMBL",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Create a Barplot for Lung (top 20 GO terms)
barplot(ego_lung, showCategory = 20, main = "GO Terms for Lung DE Genes")

# Create a Dotplot for Lung (top 20 GO terms)
dotplot(ego_lung, showCategory = 20, main = "GO Terms for Lung DE Genes")

# Create a Barplot for Blood (top 20 GO terms)
barplot(ego_blood, showCategory = 20, main = "GO Terms for Blood DE Genes")

# Create a Dotplot for Blood (top 20 GO terms)
dotplot(ego_blood, showCategory = 20, main = "GO Terms for Blood DE Genes")
