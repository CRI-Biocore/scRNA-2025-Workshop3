## **`cellranger count`** execution guide

This pipeline aligns sequencing reads from FASTQ files to a specified reference transcriptome. 
It produces a .cloupe file for visualization and analysis in the 10x Genomics Loupe Browser. 
Additionally, it generates various other outputs that are compatible with other publicly available tools,
some of which will be explored in later sections of this workshop.

The execution code is located at `/gpfs/data/cri-training/Nov_scRNAseq/analysis_demo/testData1_count.slurm`, 
and you can copy it to your home directory for execution. 
Please be aware of the execution directory and create the corresponding directory before submitting the job for execution.

### Copy the SLURM File to Workshop Lab-share Directory

```bash
cd ~
cp /gpfs/data/cri-training/Nov_scRNAseq/analysis_demo/testData1_count.slurm ~
```

`testData1_count.slurm` is the SLURM file to be executed on Randi. Depending on your operating system, you can use different tools to check and edit this file. One way to do so is via the `vi` or `vim` command to check and edit this file. For example, to open the file in `vim`, run:

```bash
vim testData1_count.slurm
```

The slurm contents are seen as below
![](./images/slurm.png)

Then press 'Esc' and type `:q!` to exit this file.

When you open this file, as shown in the presentation slide, the code is executed under the home directory `~/analysis_results1`. Therefore, we need to make sure the corresponding CSV config file is copied over to the correct directory.

### Explanation of Options

- `--id=sampl1_ln1`: This specifies the output directory for the results. The directory will be named `sampl1_ln1`.
- `--transcriptome=/path/to/refdata-gex-GRCh38-2020-A`: Specifies the path to the reference transcriptome to use. In this case, it's pointing to the human reference genome version GRCh38-2020-A.
- `--fastqs=/path/to/count`: Specifies the directory where the FASTQ files are located. These files contain the raw sequencing reads.
- `--sample=AB-EL-HE-samp1`: Indicates the specific sample(s) to process. In this case, it will process the sample named `AB-EL-HE-samp1`.
- `> testdata1_count.log 2>&1`: Redirects the output and errors to a log file (`testdata1_count.log`) for tracking and debugging.

> **Note**: There are many more options you can use with `cellranger count`. To explore all available options, run the following command:

```bash
cellranger count --help
```

### Execution Process

The execution code is located at `/gpfs/data/cri-training/Nov_scRNAseq/analysis_demo/testData1_count.slurm`, 
and you can copy it to your home directory for execution. 
Please be aware of the execution directory and create the corresponding directory before submitting the job for execution.

### Create Directory for Results

```bash
cd ~
mkdir analysis_results1
```

### Submit the SLURM Job

Once the directory is set up, you can submit the SLURM job on Randi to run with the following command:

```bash
sbatch testData1_count.slurm
```

### Check Submitted Jobs

On Randi, you can check your submitted jobs with the following command:

```bash
squeue -u your/username
```
