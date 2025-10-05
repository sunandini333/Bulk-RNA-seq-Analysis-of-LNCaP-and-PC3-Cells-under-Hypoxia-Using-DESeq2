# BulkRNAseq_DESeq2_LNCaP_PC3_Hypoxia

## 📘 Project Overview
This repository presents a **comprehensive workflow for bulk RNA-seq data analysis** performed on **LNCaP and PC3 prostate cancer cell lines** under **hypoxic and normoxic** conditions.  
It integrates preprocessing in **Bash/Linux** and downstream statistical analysis in **R**, using open-source tools such as **Trimmomatic, FastQC, MultiQC, HISAT2, and DESeq2**.

The workflow covers all stages from **raw FASTQ reads to differential expression and pathway enrichment analysis**, ensuring full reproducibility and clarity.

---

## 🧬 Workflow Summary
1. **Quality Trimming:** Adapter and low-quality read removal using Trimmomatic.  
2. **Quality Control:** FastQC for individual sample QC and MultiQC for summary reports.  
3. **Alignment:** HISAT2 alignment to the human genome (GRCh38).  
4. **Quantification:** SAMtools to sort, index, and count aligned reads.  
5. **Differential Expression:** DESeq2 normalization and DE analysis in R.  
6. **Functional Analysis:** ReactomePA, fgsea, and clusterProfiler for pathway enrichment.  
7. **Visualization:** PCA, heatmaps, and volcano plots for exploration.

---

## 💻 Bash Preprocessing Workflow

```bash
#!/bin/bash

FASTQ_DIR="./data"
GENOME_INDEX="./genome_index/grch38"
LOGFILE="./results/logfile.txt"

mkdir -p $(dirname "$LOGFILE")

FILES=(
  "LNCAP_Hypoxia_S1.fastq.gz"
  "LNCAP_Hypoxia_S2.fastq.gz"
  "LNCAP_Normoxia_S1.fastq.gz"
  "LNCAP_Normoxia_S2.fastq.gz"
  "PC3_Hypoxia_S1.fastq.gz"
  "PC3_Hypoxia_S2.fastq.gz"
  "PC3_Normoxia_S1.fastq.gz"
  "PC3_Normoxia_S2.fastq.gz"
)

for f in "${FILES[@]}"; do
  SAMPLE=$(basename "$f" .fastq.gz)
  echo "Processing $SAMPLE at $(date)" | tee -a $LOGFILE
  START=$(date +%s)

  hisat2 -q -x $GENOME_INDEX -U $FASTQ_DIR/$f | samtools sort -o results/${SAMPLE}.bam
  samtools index results/${SAMPLE}.bam

  END=$(date +%s)
  echo "Finished $SAMPLE in $((END-START)) seconds at $(date)" | tee -a $LOGFILE
  echo "-----------------------------------" | tee -a $LOGFILE
done
```

---

## 🧠 R Analysis Pipeline (DESeq2, Annotation & Visualization)

```r
# Load required libraries
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

# Load counts & metadata
counts <- read.csv("results/GSE106305_counts_matrix.csv", row.names = 1)
condition <- c(rep("LNCAP_Hypoxia", 2), rep("LNCAP_Normoxia", 2),
               rep("PC3_Hypoxia", 2), rep("PC3_Normoxia", 2))
colData <- data.frame(condition)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)

# DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds)
write.csv(as.data.frame(res), "results/DEG_results.csv")
```

---

## 📊 Visualization Outputs

| Type | Description | Example |
|------|--------------|----------|
| PCA Plot | Sample clustering | ![PCA](images/pca_plot.png) |
| Volcano Plot | Differential expression | ![Volcano](images/volcano_plot.png) |
| Heatmap | Top 500 variable genes | ![Heatmap](images/heatmap.png) |
| MultiQC Report | QC summary | ![MultiQC](images/multiqc_report.png) |

---

## 🧩 Functional Enrichment Analysis

```r
library(ReactomePA)
res_df <- read.csv("results/DEG_results.csv", row.names = 1)
gene_list <- res_df$log2FoldChange
names(gene_list) <- rownames(res_df)
gene_list <- sort(gene_list, decreasing = TRUE)

reactome_results <- gsePathway(gene_list, organism = "human", pvalueCutoff = 0.05)
dotplot(reactome_results, showCategory = 20)
```

---

## 🧾 License
This project is released under the MIT License. You are free to reuse and adapt the workflow with proper attribution.

---

## 🙌 Author
**Sunandini Chowdhury**  
BTech Biotechnology | RNA-seq Analyst  
📍 India  
🔗 [LinkedIn](https://linkedin.com) | [GitHub](https://github.com)
