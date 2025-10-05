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
 
# Converting SRA files to Fastq file:
We can now use fastq-dump to extract the contents of it into a FASTQ file. The Edwards lab at SDSU provides a nice tutorial for how to use fastq-dump here: https://edwards.sdsu.edu/research/fastq-dump/. This will create a file called SRRXXXXXXX_pass.fastq.gz inside the directory where fastq-dump was called. The “.gz” just means that the file is compressed. To perform it in 20 SRA files, we wrote a Python script which will generate a fastq file for each of the 20 SRA files in a subdirectory called tempfa.
 
# REAL TIME UNDERSTANDING OF CONVERTING SRA FILE TO FASTQ FILE:
By putting command
watch -n 1 ls -la tempfa/
we can easily understand the real time conversion of SRA to Fastq in the subdirectory called tempfa to crosscheck our progress.
 
 

# Concatenation of files into total 8 files 
Each of the LNCaP samples was associated with four SRA runs, which means that we obtained four resulting FASTQ files for each sample after running fastq_download.py. For each sample, we should concatenate the four files into a single FASTQ file by using the command cat. Below, I perform the concatenation for each of the LNCaP samples:
cat SRR7179504_pass.fastq.gz SRR7179505_pass.fastq.gz SRR7179506_pass.fastq.gz SRR7179507_pass.fastq.gz  > LNCAP_Normoxia_S1.fastq.gz
cat SRR7179508_pass.fastq.gz SRR7179509_pass.fastq.gz SRR7179510_pass.fastq.gz SRR7179511_pass.fastq.gz  > LNCAP_Normoxia_S2.fastq.gz
cat SRR7179520_pass.fastq.gz SRR7179521_pass.fastq.gz SRR7179522_pass.fastq.gz SRR7179523_pass.fastq.gz  > LNCAP_Hypoxia_S1.fastq.gz
cat SRR7179524_pass.fastq.gz SRR7179525_pass.fastq.gz SRR7179526_pass.fastq.gz SRR7179527_pass.fastq.gz  > LNCAP_Hypoxia_S2.fastq.gz
In contrast, there is only one FASTQ file for each of the PC3 samples. We can just rename them from their SRR identifiers to their real sample names using mv:
mv SRR7179536_pass.fastq.gz PC3_Normoxia_S1.fastq.gz
mv SRR7179537_pass.fastq.gz PC3_Normoxia_S2.fastq.gz
mv SRR7179540_pass.fastq.gz PC3_Hypoxia_S1.fastq.gz
mv SRR7179541_pass.fastq.gz PC3_Hypoxia_S2.fastq.gz
We won’t need the individual SRA runs anymore, so we can remove them all using the command rm SRR*, which removes all the files in the folder that begin with “SRR”. Now, the folder should contain a total of 8 FASTQ files: 4 for LNCaP and 4 for PC3. We are ready to begin aligning these FASTQ files to the reference genome!

 

 

# Performing FastQC:
FastQC is a quality control tool for next-generation sequencing (NGS) data, especially FASTQ files.When you run it (like in your command), it checks the raw sequencing reads and produces reports that help you decide if your data is good enough to use for downstream analysis (alignment, quantification, etc.) or if it needs trimming/filtering.
🔹 What does FastQC check?
It generates both HTML reports and detailed stats. The main checks are:
1.	Per base sequence quality → Phred quality score for each base across all reads.
o	High = good sequencing.
o	Low = errors → may need trimming.
2.	Per sequence quality scores → Are most reads of high quality or not?
3.	Per base sequence content → Balance of A, T, G, C at each position.
o	Should be roughly even (except in special cases like RNA-seq with biases).
4.	Per base GC content → Detects unusual GC patterns.
5.	Sequence length distribution → Are all reads the same length? (important for alignment).
6.	Overrepresented sequences → Adapters or contamination.
7.	Adapter content → Shows if sequencing adapters are still present.
8.	K-mer content → Detects unexpected sequence motifs (bias/contamination).

🔹 Why is FastQC important?
•	First step in RNA-seq or DNA-seq analysis → ensures raw reads are usable.
•	Detects problems early → e.g., adapter contamination, poor quality, low diversity.
•	Guides preprocessing → tells you if you need trimming, filtering, or removing adapters before alignment.

 

 
 By installing the fastqc at tempfa directory you can perform fastqc in the files having extension .fastq.gz in tempfa directory and saving the result in a new subdirectory names fastqc_results.


# Performing MultiQC:
running a summary report that combines the results from multiple bioinformatics tools (like FastQC, HISAT2, featureCounts, etc.) into one single HTML file that’s easy to read and interpret.
Here, When you analyze RNA-seq data, you usually generate many output files  for example:
20 FastQC reports (.html + .zip files), so Reading each one manually is painful. So, MultiQC scans all these outputs and automatically creates one combined report.




 


 
After installing multiqc we need to run all the 8 fastqc results stored in the subdirectory called fastqc_results into a output directory called multiqc_report.




# Trimming:
By observing the fastqc and multiqc results we can understand that the data doesn’t need any trimming after noticing the basic statistics of each cell lines itself but still we performed trimming and compared the two multiqc result before trimming and after trimming to understand with data should we proceed for downstream processing.
After installing Trimmomatic we get Trimmomatic—0.39 zip file grom the web and then after unzipping it we run the code for trimming for total 8 times on 8 different files.







 
 
# Fastqc  and multiqc on trimmed files:
Same as the step number 6 and 7 but this time it is perfomed on the trimmed files having extension .fastq and saved in a output subdirectory names fastqc_results_trimmed and for multiqc also likewise we did.
 

 

# Quality Control and Trimming Decision
After performing quality assessment using FastQC and aggregating the results with MultiQC, we compared the reports generated from both the raw FASTQ files and the trimmed FASTQ files.
The overall quality metrics—including per-base sequence quality, GC content, adapter content, and sequence duplication levels—showed no significant improvement following trimming.
Therefore, to preserve read length and maintain data integrity, we proceeded with the untrimmed (raw) FASTQ files for all downstream analyses.
However, if in your dataset trimming results in improved read quality (e.g., removal of adapter contamination or low-quality bases), you should proceed with the trimmed FASTQ files for subsequent analysis steps.


<img width="940" height="352" alt="image" src="https://github.com/user-attachments/assets/fb011dd3-5f84-4ad5-b6a9-a7d75714810f" />





