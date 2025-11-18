###################################################################
### GSE Title: Gene expression data from healthy, adjacent normal
###            and tumor colon cells
### GSE ID: GSE44076
### Comparison: Colorectal Tumor vs. Adjacent Normal Tissue
###################################################################

# --- 0. Setup: Define Output Path ---

output_path <- "D:/DOWNLOADS"
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}
cat("All output files will be saved to:", output_path, "\n")


# --- 0. Setup: Install and Load All Packages ---

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# List of Bioconductor packages
bioc_packages <- c("GEOquery", "DESeq2", "EnhancedVolcano")

# Install Bioconductor packages
cat("Checking and installing Bioconductor packages...\n")
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# List of CRAN packages
cran_packages <- c("pheatmap", "RColorBrewer", "ggplot2", "ggrepel", "magrittr", "dplyr")

# Install CRAN packages
cat("Checking and installing CRAN packages...\n")
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load all required libraries
cat("Loading all required libraries...\n")
library(GEOquery)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(magrittr)
library(dplyr)
library(EnhancedVolcano)


###################################################################
### PART 1: Data Loading and Normalization (GSE44076)
###################################################################

cat("--- Starting Part 1: Data Loading ---\n")

# --- 1.1. Download Data from GEO ---
gse_id <- "GSE44076"
gse <- getGEO(gse_id, GSEMatrix = TRUE)
data <- gse[[1]]

# --- 1.2. Extract and Filter Data ---
expression_data <- exprs(data)
pheno_data <- pData(data)
cat("Using column 'source_name_ch1' for sample groups.\n")
subtype <- pheno_data$`source_name_ch1`

# --- *** CRITICAL ERROR FIX *** ---
# Using the *exact* labels you found with the debug script.
correct_labels_to_keep <- c("Primary colon adenocarcinoma cells", "Normal distant colon mucosa cells")
keep_samples <- which(subtype %in% correct_labels_to_keep)

# Check if any samples were found.
if (length(keep_samples) == 0) {
  stop("FATAL ERROR: Could not find 'Primary colon adenocarcinoma cells' or 'Normal distant colon mucosa cells'.")
}

# Filter expression data and subtype labels
expression_data <- expression_data[, keep_samples]
subtype <- subtype[keep_samples]

# Rename labels to "Tumor" and "Normal" for simplicity in our plots
cat("Renaming labels for analysis...\n")
subtype[subtype == "Primary colon adenocarcinoma cells"] <- "Tumor"
subtype[subtype == "Normal distant colon mucosa cells"] <- "Normal"
# --- *** END OF FIX *** ---

cat("Sample subtypes identified and filtered:\n")
print(table(subtype))
cat("Data loaded. Normalization will be done by DESeq2.\n")


###################################################################
### PART 2: Run DESeq2 and Create QC Plots
###################################################################

cat("--- Starting Part 2: DESeq2 and QC Plots ---\n")

# --- 2.1. Prepare DESeq2 Object ---
cat("Preparing DESeq2 data object...\n")
# Ensure 'subtype' is a factor and set 'Normal' as the base level
col_data <- data.frame(condition = factor(subtype))
col_data$condition <- relevel(col_data$condition, ref = "Normal")
dds <- DESeqDataSetFromMatrix(countData = round(expression_data), 
                              colData = col_data, 
                              design = ~ condition)
# Run the main DESeq analysis
cat("Running DESeq analysis...\n")
dds <- DESeq(dds)

# --- 2.2. NEW PLOT 1: Dispersion Plot ---
cat("Generating and saving Dispersion Plot...\n")
png(file.path(output_path, "deseq2_dispersion_plot.png"), width=800, height=600)
plotDispEsts(dds)
dev.off()

# --- 2.3. Get Normalized Data for QC Plots ---
cat("Applying variance stabilizing transformation...\n")
vsd <- vst(dds, blind = TRUE)

# Create annotation data frame for all plots
annotation_df <- data.frame(Subtype = col_data$condition)
rownames(annotation_df) <- colnames(vsd)
annotation_colors <- list(Subtype = c(Tumor = "red", Normal = "blue"))


# --- 2.4. NEW PLOT 2: PCA Plot ---
cat("Generating and saving PCA Plot...\n")
pca_plot <- plotPCA(vsd, intgroup = "condition") + 
  labs(title = "PCA of Tumor vs. Normal Samples (GSE44076)") +
  theme_minimal()
ggsave(file.path(output_path, "pca_plot.png"), plot = pca_plot)


# --- 2.5. NEW PLOT 3: Sample-to-Sample Distance Heatmap ---
cat("Generating and saving Sample Distance Heatmap...\n")
sample_dist_matrix <- dist(t(assay(vsd)))
sample_dist_matrix_plot <- as.matrix(sample_dist_matrix)
png(file.path(output_path, "sample_distance_heatmap.png"), width=1000, height=1000)
pheatmap(sample_dist_matrix_plot,
         clustering_distance_rows = sample_dist_matrix,
         clustering_distance_cols = sample_dist_matrix,
         main = "Sample-to-Sample Distance Heatmap (GSE44076)",
         annotation_col = annotation_df,
         annotation_colors = annotation_colors,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()


# --- 2.6. PLOT 4: Top 50 Variable Genes Heatmap ---
cat("Generating and saving Top 50 Variable Gene Heatmap...\n")
top_50_genes_indices <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
heatmap_data <- assay(vsd)[top_50_genes_indices, ]

pheatmap(heatmap_data, 
         annotation_col = annotation_df, 
         scale = "row",
         show_rownames = FALSE,
         cluster_rows = TRUE, # Explicitly cluster genes
         main = "Heatmap of Top 50 Variable Genes (GSE44076)",
         annotation_colors = annotation_colors,
         filename = file.path(output_path, "heatmap_top50_variable_genes.png"))


###################################################################
### PART 3: DGE Results and Visualization
###################################################################

cat("--- Starting Part 3: DGE Results and Visualization ---\n")

# --- 3.1. Extract DGE Results ---
results <- results(dds, contrast = c("condition", "Tumor", "Normal"))


# --- 3.2. Save Significant Gene Results (CSV) ---
results_df <- as.data.frame(results)
significant_genes <- results_df[!is.na(results_df$padj) & results_df$padj < 0.05, ]
cat("Number of significant DEGs found:", nrow(significant_genes), "\n")
cat("Saving significant gene list to significant_genes_tumor_vs_normal.csv...\n")
write.csv(significant_genes, 
          file = file.path(output_path, "significant_genes_tumor_vs_normal.csv"))


# --- 3.3. NEW PLOT 5: Heatmap of Top 50 Significant Genes ---
cat("Generating and saving Top 50 Significant Gene Heatmap...\n")
# Sort by adjusted p-value
top_50_sig_genes <- significant_genes %>%
  arrange(padj) %>%
  head(50)
# Get the normalized expression data for these genes
top_50_sig_heatmap_data <- assay(vsd)[rownames(top_50_sig_genes), ]

pheatmap(top_50_sig_heatmap_data, 
         annotation_col = annotation_df, 
         scale = "row", # Scale genes
         cluster_rows = TRUE,
         show_rownames = TRUE, # Show the gene names
         fontsize_row = 8,
         main = "Heatmap of Top 50 Significant DEGs (GSE44076)",
         annotation_colors = annotation_colors,
         filename = file.path(output_path, "heatmap_top50_significant_genes.png"))


# --- 3.4. PLOT 6: Enhanced Volcano Plot ---
cat("Generating and saving Enhanced Volcano Plot...\n")
ev_plot <- EnhancedVolcano(results_df,
                           lab = rownames(results_df),
                           x = 'log2FoldChange',
                           y = 'padj',
                           title = 'Volcano Plot (Tumor vs. Normal)',
                           subtitle = 'GSE44076',
                           pCutoff = 0.05,
                           FCcutoff = 1.0,
                           pointSize = 2.0,
                           labSize = 4.0,
                           legendPosition = 'bottom',
                           colAlpha = 0.5)
ggsave(file.path(output_path, "enhanced_volcano_plot.png"), plot = ev_plot, width = 12, height = 10)


# --- 3.5. PLOT 7: MA Plot ---
cat("Generating and saving MA Plot...\n")
volcano_df <- as.data.frame(results)
volcano_df$significant <- ifelse(volcano_df$padj < 0.05, "Significant", "Not significant")

ma_plot <- ggplot(volcano_df, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(alpha = 0.7, color = "gray") +
  geom_point(data = subset(volcano_df, padj < 0.05), aes(color = significant), size = 2) +
  scale_color_manual(values = c("red")) +
  labs(title = "MA Plot (GSE44076)", x = "Log2 Mean Expression", y = "Log2 Fold Change") +
  theme_minimal()

ggsave(file.path(output_path, "ma_plot.png"), plot = ma_plot)


# --- 3.6. NEW PLOT 8: Boxplot of Top Gene ---
cat("Generating and saving boxplot for the #1 significant gene...\n")
# Find the #1 gene (smallest padj)
top_gene_id <- rownames(significant_genes)[1]
# Get the normalized count data for this gene
top_gene_counts <- plotCounts(dds, gene = top_gene_id, intgroup = "condition", returnData = TRUE)

# Create the boxplot
top_gene_boxplot <- ggplot(top_gene_counts, aes(x = condition, y = count, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.7) +
  scale_fill_manual(values = c("Normal" = "blue", "Tumor" = "red")) +
  labs(title = paste("Expression of Top Gene:", top_gene_id),
       x = "Sample Group",
       y = "Normalized Count") +
  theme_minimal()
ggsave(file.path(output_path, "top_gene_boxplot.png"), plot = top_gene_boxplot)


cat("\n#################################################\n")
cat("Analysis complete. All 8 plots and CSV saved to:", output_path, "\n")
cat("#################################################\n")
