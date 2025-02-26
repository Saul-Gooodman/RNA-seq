# Load necessary libraries
library(DESeq2)

# Extract results for a pairwise comparison
results_lung <- results(dds, contrast = c("condition", "Lung_WT_Case", "Lung_WT_Control"))

# Filter significant genes (padj < 0.05 and absolute log2FoldChange greater than 1)
significant_genes_lung <- results_lung[
  !is.na(results_lung$padj) & results_lung$padj < 0.05 & 
    !is.na(results_lung$log2FoldChange) & abs(results_lung$log2FoldChange) > 1, 
]

# Count the number of upregulated and downregulated genes
upregulated_genes <- subset(significant_genes_lung, log2FoldChange > 0)
downregulated_genes <- subset(significant_genes_lung, log2FoldChange < 0)

# Print the number of upregulated genes
cat("Number of upregulated genes:", nrow(upregulated_genes), "\n")
# Print the number of downregulated genes
cat("Number of downregulated genes:", nrow(downregulated_genes), "\n")

# Extract normalized expression data for genes of interest 
genes_of_interest <- c("Gbp2", "Fcer1g", "Tgfp1")
vsd_subset <- vsd[genes_of_interest, ]

# View the normalized expression data for the genes of interest
print(vsd_subset)
