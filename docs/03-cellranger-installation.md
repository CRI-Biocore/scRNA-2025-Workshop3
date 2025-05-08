# Cell Ranger Installation Guide

This guide provides step-by-step instructions for installing the **Cell Ranger** software suite, 
developed by 10x Genomics, which is required for processing single-cell RNA-seq data.

> 💡 **Note:** The Cell Ranger environment has already been set up on the **Randi HPC**. 
Users on this system can skip installation and proceed directly to data processing on **Randi HPC**. 
Please load the module or activate the appropriate conda environment as configured on the cluster.

---

## 1. System Requirements

Before installation, ensure your system meets the following requirements:

- **Operating System**: Linux (Ubuntu 16.04+, CentOS 7+)
- **CPU**: 8 or more cores recommended
- **RAM**: At least 32 GB (higher for large datasets)
- **Disk Space**: 1 TB or more recommended
- **Python**: Comes bundled with Cell Ranger (no separate installation required)

---

## 2. Download Cell Ranger

1. Visit the official [10x Genomics Cell Ranger download page](https://www.10xgenomics.com/support/software/cell-ranger).
2. Choose the desired version (e.g., `cellranger-7.2.0`).
3. Download and extract the tarball using the following commands:

```
wget -O cellranger-7.2.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.2.0.tar.gz"
tar -xzvf cellranger-7.2.0.tar.gz
```

## 3: Set Up Environment for cellRanger

After downloading and extracting Cell Ranger, you need to configure your system environment 
so the `cellranger` command is recognized in any directory.

---

## Add Cell Ranger to Your PATH

Temporarily set the PATH in your current terminal session:

```
export PATH=$PATH:/your/path/to/cellranger-7.2.0
```

## Verify cellRanger Installation

After configuring your environment, verify that **cellRanger** is correctly installed and 
accessible from your terminal.

e.g. you can check cellranger version via running the following command:

```
cellranger --version
```

If the installation is successful, you should see output similar to:

```
cellranger 7.2.0
```

For our Randi HPC, the cellRanger path has been integrated into the module system. 
You can access the configured computational environment by loading the module using the following command:

```
module load workshop3
cellranger --version
```
