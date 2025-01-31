# Load necessary libraries
library(DESeq2)
library(ggplot2)

# Load the count matrix and sample information
counts_file <- "~/Desktop/counts_matrix.txt"
counts_data <- read.table(counts_file, header = TRUE, row.names = 1)

# Clean column names by removing the path part
colnames(counts_data) <- gsub(".*(SRR\\d+_sorted\\.bam)$", "\\1", colnames(counts_data))

# Remove columns related to the genome (Chr, Start, End, Strand, Length)
counts_data <- counts_data[, !(colnames(counts_data) %in% c("Chr", "Start", "End", "Strand", "Length"))]

# Create a col_data data frame containing the experimental group information for each sample
col_data <- data.frame(
  condition = factor(c(
    "Lung_WT_Case", "Lung_WT_Case", "Lung_WT_Case", "Lung_WT_Case", "Lung_WT_Case", 
    "Lung_WT_Control", "Lung_WT_Control", "Lung_WT_Control", 
    "Blood_WT_Case", "Blood_WT_Case", "Blood_WT_Case", "Blood_WT_Case", "Blood_WT_Case", 
    "Blood_WT_Control", "Blood_WT_Control", "Blood_WT_Control"
  )),
  row.names = c(
    "SRR7821918", "SRR7821919", "SRR7821920", "SRR7821921", "SRR7821922", 
    "SRR7821937", "SRR7821938", "SRR7821939", 
    "SRR7821949", "SRR7821950", "SRR7821951", "SRR7821952", "SRR7821953", 
    "SRR7821968", "SRR7821969", "SRR7821970"
  )
)

# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = col_data, design = ~ condition)

# Run the DESeq2 analysis
dds <- DESeq(dds)

# Get the results of differential expression analysis
res <- results(dds)

# Normalize the data using vst() to remove the dependence of variance on the mean
vsd <- vst(dds, blind = TRUE)

# Generate PCA plot
pcaData <- plotPCA(vsd, returnData = TRUE)  # Fix undefined pcaData
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  ggtitle("PCA of Gene Expression Data")

# Create the volcano plot
plot(res$log2FoldChange, -log10(res$padj), pch = 20, 
     col = ifelse(res$padj < 0.05, "red", "black"), 
     xlab = "Log2 Fold Change", ylab = "-log10(padj)")
