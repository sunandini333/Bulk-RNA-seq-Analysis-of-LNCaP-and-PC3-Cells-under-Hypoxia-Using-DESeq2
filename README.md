# Introduction
a comprehensive workflow for analyzing bulk RNA-seq data from LNCaP and PC3 prostate cancer cell lines under hypoxic and normoxic conditions. The pipeline begins with raw FASTQ files, which are preprocessed using Trimmomatic to remove adapters and low-quality bases. Quality assessment is performed with FastQC, and MultiQC provides an aggregated overview of all samples, ensuring high-quality data for downstream analysis. High-quality reads are aligned to the human genome using HISAT2, generating count matrices for each sample. In R, the count data are normalized, and differential gene expression is analyzed using DESeq2. Genes are annotated with biotypes, and low-count or zero-count genes are filtered to improve statistical robustness. Principal component analysis (PCA) and heatmaps are generated to assess sample clustering and variation, while volcano plots visualize differentially expressed genes between conditions. The workflow includes detailed gene-level analyses, allowing visualization of normalized expression for individual genes across experimental groups. Functional enrichment and pathway analyses are conducted using ReactomePA and fgsea, highlighting pathways significantly altered under hypoxia. Genes are mapped from Ensembl IDs to Entrez IDs to enable pathway and gene set enrichment analyses, with visualization of top enriched pathways and hallmark gene sets. This repository integrates Bash/Linux commands for preprocessing steps—trimming, quality control, and alignment—with R-based analyses for statistical modeling, visualization, and pathway enrichment. Example data, annotated counts, and scripts are included to facilitate reproducibility and hands-on learning. This pipeline can be easily adapted to other RNA-seq datasets, making it a practical guide for researchers aiming to perform end-to-end RNA-seq analysis from raw reads to biologically meaningful insights. By combining widely used tools such as Trimmomatic, FastQC, MultiQC, HISAT2, and DESeq2, this workflow ensures high-quality, reproducible results. It is ideal for anyone seeking a structured, transparent, and detailed approach to bulk RNA-seq analysis in human cell lines under different experimental conditions.
# Obtaining raw data from GEO
The dataset that we will be working with comes from Guo et al. Nature Communications 2019. To find the raw sequencing data, we can navigate through the Gene Expression Omnibus (GEO) using the accession number provided in the publication: GSE106305. The GEO page corresponding to this accession number is https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106305. On this webpage, there is information about where the sequencing data for each sample in the study is stored.
Our goal is to find differentially expressed genes in response to hypoxia for the LNCaP and PC3 cell lines. Therefore, we will select the control samples for both cell lines (Empty_Vector for LNCaP and siCtrl for PC3) in conditions of either normoxia or hypoxia. The specific samples we need to download are outlined in the table below:
Sample Name	GSM Identifier	SRA Identifier (SRX)	SRA Runs (SRR, download these)
LNCaP_RNA-Seq_Empty_Vector_Normoxia_rep1	GSM3145509	SRX4096735	SRR7179504, SRR7179505, SRR7179506, and SRR7179507
LNCaP_RNA-Seq_Empty_Vector_Normoxia_rep2	GSM3145510	SRX4096736	SRR7179508, SRR7179509, SRR7179510, and SRR7179511
LNCaP_RNA-Seq_Empty_Vector_Hypoxia_rep1	GSM3145513	SRX4096739	SRR7179520, SRR7179521, SRR7179522, and SRR7179523
LNCaP_RNA-Seq_Empty_Vector_Hypoxia_rep2	GSM3145514	SRX4096740	SRR7179524, SRR7179525, SRR7179526, and SRR7179527
PC3_RNA-Seq_siCtrl_Normoxia_rep1	GSM3145517	SRX4096743	SRR7179536
PC3_RNA-Seq_siCtrl_Normoxia_rep2	GSM3145518	SRX4096744	SRR7179537
PC3_RNA-Seq_siCtrl_Hypoxia_rep1	GSM3145521	SRX4096747	SRR7179540
PC3_RNA-Seq_siCtrl_Hypoxia_rep2	GSM3145522	SRX4096748	SRR7179541
The NCBI stores sequencing data in the form of Sequence Read Archive (SRA) files. Each of the samples is associated with a set of SRA accession numbers, indicated above. First, we need to download the SRA runs for each sample. Then, we will use the SRA files to generate FASTQ files. The FASTQ file is the data format required for bulk RNA-sequencing analysis.
# Download SRA files using SRA tools
To download the SRA files onto our machine, we will use the NCBI’s SRA toolkit. In linux, you can type: sudo apt install sra-toolkit in your command line to install the toolkit. You can read more about SRA toolkit here: https://www.ncbi.nlm.nih.gov/books/NBK242621/.
The toolkit works by first using the prefetch command to download the SRA file associated with the specified SRA ID. The SRA file contains a set of “instructions” for downloading the sequencing data associated with the SRA ID from NCBI.

 
We downloaded all the SRA files by making a new directory rnaseq_project. We can see by giving ls command.
 






