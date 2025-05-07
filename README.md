# CRI Bioinformatics Core Workshop Series - 10x Genomics scRNA-seq data analysis workshop

**This bioinformatics learning series, hosted by the CRI Bioinformatics Core, provides a comprehensive, 
hands-on introduction to single-cell RNA sequencing (scRNA-seq) data analysis. 
Participants will be guided through each step of the typical scRNA-seq workflow using widely adopted tools 
such as 10x Genomics cellranger, seurat, and DoubletDecon.**

**This workshop will cover the topics listed below. 
You can access the presentation slides by [downloading Workshop Slides (PDF)](./docs/Presentation_June2025_scRNAseq.pdf) **

Please share your feedback on topics you'd like to see covered, or any questions you have 
for our workshop learning series, either through [this survey link](https://mycri.cri.uchicago.edu/educations/trainings/77/survey/) 
or by emailing me at yli22@bsd.uchicago.edu.

### 1. High Performance Computing (HPC)

CRI provides access to its high-performance computing (HPC) system, randi, for all registered participants in today's workshop training session. 
The complete computational environment has been pre-configured on randi to support scRNA-seq data pre-processing, quality control, 
and downstream analysis.

To access randi, use the SSH command below with your assigned BSDID and password:

```
ssh your_username@randi.cri.uchicago.edu
```

After logging in successfully, you should see a screen similar to the following:
![Login Screenshot](./images/randi-login-screenshot.png)

**💬 CRI HPC Randi Support**

You are welcome to join the **CRI HPC Randi Slack user group** to ask questions or get support related to Randi access:

🔗 [Join the Slack Group](https://join.slack.com/t/criscientific-dzi9891/shared_invite/zt-2kghy4392-1ELPfgn8pL5BcXk4oF9D4g)

---

🕒 Additionally we are offering **joint in-person office hours every Tuesday afternoon** at **Peck Pavilion N161**:

  - **HPC Team:** 12:00 PM – 2:30 PM  
  - **Bioinformatics Team:** 12:30 PM – 3:30 PM  

Feel free to stop by with any questions — we are happy to help!

### 2. Introduction to Single-Cell RNA Sequencing

#### 📊 Overview of 10x Genomics scRNA-seq Technology

10x Genomics Chromium platform enables **high-throughput single-cell RNA sequencing (scRNA-seq)**. 
It captures thousands to tens of thousands of individual cells per run, isolating and barcoding transcripts from each cell.

**Key Concepts:**

- 🧪 **GEMs (Gel Beads-in-Emulsion):** Each droplet contains one cell and one barcoded bead.
- 🧬 **Cell Barcodes:** Identify which cell each transcript came from.
- 🔁 **UMIs (Unique Molecular Identifiers):** Help quantify gene expression by correcting for PCR duplicates.
- 🎯 Supports both **3' and 5' gene expression profiling**.

---

#### 🧬 10x Genomics scRNA-seq Workflow

1. **Sample Preparation:** Generate a single-cell suspension from tissue or culture.
2. **Chromium Partitioning:** Cells are encapsulated in droplets with barcoded beads.
3. **Reverse Transcription:** Transcripts are barcoded and converted to cDNA inside droplets.
4. **Library Construction:** Amplification and preparation for Illumina


### 3. **[biocore scRNA-seq Data Analysis Workflow]**

This section outlines the complete computational pipeline used by the CRI Bioinformatics Core 
for single-cell RNA-seq analysis—from raw sequencing data to downstream interpretation and visualization.

![](./images/scRNA-biocore-workflow.png)


### 4. **[Data Preprocessing with cellRanger](./docs/03-cellranger.md)**
   + Aligning reads, demultiplexing barcodes, and generating feature-barcode matrices from raw sequencing data.

### 5. **[Quality Control and Filtering Using Seurat](./docs/04-qc-seurat.md)**  
   + Identifying low-quality cells and applying standard filtering criteria

### 6. **[Detecting Doublets with DoubletDecon](./docs/05-doublets.md)** 
   + Identifying and removing artificial cell doublets/multiplets to improve data fidelity

### 7. **[Data Normalization and Integration](./docs/06-normalization-integration.md)**  
   + Seurat-based normalization and integration across samples or conditions

### 8. **[Clustering and Cell Type Identification](./docs/07-clustering.md)**  
   + Performing unsupervised clustering and assigning putative cell identities based on marker gene expression

### 9. **[Visualization with UMAP/TSNE](./docs/08-umap.md)**  
   + Projecting high-dimensional scRNA-seq data into two dimensions for visual interpretation

### 10. **[Differential Expression Analysis](./docs/09-de-analysis.md)**  
   + Identifying marker genes and performing comparisons across clusters or conditions

