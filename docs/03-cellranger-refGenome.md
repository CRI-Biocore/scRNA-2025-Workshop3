# Cell Ranger Reference Genome Installation Guide

To enable read alignment and quantification, `cellranger` requires an appropriate reference genome file. 
10x Genomics provides pre-built reference packages for common species such as human and mouse.

In this tutorial, all test data come from **human samples**, so we will demonstrate how to install 
the human reference genome for use with `cellranger` only.

> **Note:** For multiome analysis using `cellranger-arc`, 
a different reference package is required. That installation is not covered in this document.

---

## 📦 Cell Ranger Reference for `count` / `multi`

The current human genome reference version recommended for use with `cellranger` is:

- **Human reference (GRCh38) – 2020-A**

### 🔽 Download

Run the following command to download the reference package:

```bash
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
```

### 📂 Extract

After downloading, extract the contents:

```bash
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz
```

You will then have a directory (e.g., `refdata-gex-GRCh38-2020-A`) containing the pre-built genome index and gene annotation. This directory can be passed to `cellranger` commands using the `--transcriptome` option.

---

✅ After downloading and extracting the reference files, you are ready to proceed with data processing using `cellranger`.

➡️ For step-by-step instructions on running `cellranger count`, refer to [03-cellranger-count.md](./03-cellranger-count.md).

---

### Mouse Reference Genome

The current mouse genome reference version recommended for use with `cellranger` is:

- **Mouse reference (GRCm38) – 2020-A**

#### 🔽 Download

Run the following command to download the reference package:

```bash
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm38-2020-A.tar.gz
```

#### 📂 Extract

After downloading, extract the contents:

```bash
tar -xzvf refdata-gex-GRCm38-2020-A.tar.gz
```

As with the human reference, you will have a directory (e.g., `refdata-gex-GRCm38-2020-A`) that can be used in `cellranger count` or `cellranger multi`.


