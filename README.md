# DESeq2-DGE-Analysis-of-GSE44076-with-QC-PCA-Heatmaps-and-EnhancedVolcano-Visualization
DESeq2 DGE Analysis of GSE44076 (Tumor vs. Normal) with QC (PCA, Heatmaps) and EnhancedVolcano Visualization
# 🦠 RNA-Seq DGE Pipeline: Colorectal Tumor vs. Normal (GSE44076)

This R script automates a comprehensive bioinformatics pipeline for **Differential Gene Expression (DGE)** analysis of RNA-Seq count data, specifically comparing **Colorectal Tumor** samples against **Adjacent Normal Tissue** using the publicly available **GSE44076** dataset from the Gene Expression Omnibus (GEO).

The pipeline uses the robust **`DESeq2`** package for DGE and generates an extensive set of eight diagnostic (QC) and results plots, ensuring a complete and reproducible analysis workflow.

## 🚀 Key Features

* **RNA-Seq Specific DGE:** Uses the **`DESeq2`** package, which is optimized for count data and provides accurate normalization and dispersion estimation.
* **Automated Data Retrieval:** Downloads raw count data and phenotype information directly from GEO using **`GEOquery`**.
* **Targeted Filtering:** Filters the dataset to compare `Primary colon adenocarcinoma cells` (renamed to **Tumor**) against `Normal distant colon mucosa cells` (renamed to **Normal**).
* **Comprehensive QC & Visualization:** Generates a suite of eight plots covering quality control and final results, including **PCA**, **Heatmaps**, **Volcano Plots**, and **Dispersion Estimates**.
* **Reproducible Output:** All results (CSV and PNG plots) are automatically saved to a dedicated output directory.

---

## 🔬 Analysis Overview

| Component | Method / Test | Purpose |
| :--- | :--- | :--- |
| **DGE** | `DESeq2` (Negative Binomial GLMs) | Identifies genes significantly up or down-regulated between Tumor and Normal tissues. |
| **Normalization** | **Median of Ratios** (within `DESeq2`) | Corrects for sequencing depth and RNA composition differences. |
| **QC Transformation** | **Variance Stabilizing Transformation (VST)** | Prepares data for global QC plots (PCA, Heatmaps) by stabilizing variance across the mean. |
| **Comparison** | Tumor vs. Normal | Identifies the unique gene signature associated with colorectal tumor pathology. |

---

## 🛠️ Prerequisites and Setup

### 📦 Packages

This script requires several Bioconductor and CRAN packages. The script includes an automated check and installation routine for all necessary libraries:
* **Bioconductor:** `GEOquery`, `DESeq2`, `EnhancedVolcano`
* **CRAN:** `pheatmap`, `RColorBrewer`, `ggplot2`, `ggrepel`, `magrittr`, `dplyr`

### ⚙️ Execution

1.  **Download** the `DESeq2 DGE Analysis of GSE44076 (Tumor vs. Normal) with QC (PCA, Heatmaps) and EnhancedVolcano Visualization.R` file.
2.  **Optional:** Modify the `output_path` variable at the beginning of the script to your desired saving location.
    ```R
    output_path <- "D:/DOWNLOADS" # Change this path
    ```
3.  **Execute** the script in your R environment:
    ```R
    source("DESeq2 DGE Analysis of GSE44076 (Tumor vs. Normal) with QC (PCA, Heatmaps) and EnhancedVolcano Visualization.R")
    ```

---

## 📁 Output Files (8 Plots + 1 CSV)

All output files are saved to the specified `output_path`.

### Statistical Results

| Filename | Type | Description |
| :--- | :--- | :--- |
| `significant_genes_tumor_vs_normal.csv` | CSV | Table containing all genes with an adjusted p-value (padj) $< 0.05$. Includes log2FoldChange, padj, and baseMean. |

### Visualization and QC Plots

| Filename | Analysis Stage | Description |
| :--- | :--- | :--- |
| `deseq2_dispersion_plot.png` | QC | **Dispersion Plot** showing estimated, fitted, and final dispersion values. |
| `pca_plot.png` | QC | **Principal Component Analysis (PCA)** plot for visual separation of Tumor vs. Normal samples. |
| `sample_distance_heatmap.png` | QC | **Sample-to-Sample Distance Heatmap** to assess sample relatedness and check for batch effects. |
| `heatmap_top50_variable_genes.png` | QC | **Heatmap of Top 50 Most Variable Genes** (VSD-transformed data). |
| `heatmap_top50_significant_genes.png` | Results | **Heatmap of Top 50 Significant DEGs** (ranked by padj). |
| `enhanced_volcano_plot.png` | Results | **Volcano Plot** illustrating log2FC vs. $\log_{10}(P_{\text{adj}})$, highlighting significant DEGs. |
| `ma_plot.png` | Results | **MA Plot** (Mean-Difference Plot) showing log2FC vs. $\log_2(\text{Mean Expression})$, colored by significance. |
| `top_gene_boxplot.png` | Results | **Boxplot** illustrating the normalized count distribution of the single most significant DEG. |
