# Assignment 3: RNA-Seq Analysis Pipeline

## Overview
In this assignment, you will build a **Snakemake workflow** to analyze RNA-Seq data. You will perform quality control, trim reads, align them to a reference genome, quantify gene expression, and perform exploratory analysis. This assignment builds on the skills you learned in Assignment 2 and introduces new Snakemake features.

## Learning Objectives
By completing this assignment, you will:
- Build an RNA-Seq analysis pipeline from raw reads to gene counts
- Use external Python and R scripts in your workflow
- Generate gene expression matrices and visualizations
- Prepare data for differential expression analysis

## Prerequisites
Before starting this assignment, ensure you have the following installed:
- Python
- [Snakemake](https://snakemake.readthedocs.io/)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://multiqc.info/)
- [Trim Galore](https://github.com/FelixKrueger/TrimGalore)
- [STAR](https://github.com/alexdobin/STAR)
- [SAMtools](http://www.htslib.org/)
- [Subread](http://subread.sourceforge.net/) (for featureCounts)
- R with ggplot2, DESeq2, and pheatmap packages

You can install all tools using conda/mamba:
```bash
mamba create -n rnaseq-env -c bioconda -c conda-forge \
    snakemake fastqc multiqc trim-galore star samtools \
    subread r-base r-ggplot2 r-pheatmap bioconductor-deseq2
conda activate rnaseq-env
```

## Data
For this assignment, you will need RNA-Seq data with the following structure:
- Paired-end FASTQ files in `data/fastq/` (`.fastq.gz` format)
  - Naming convention: `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`
  - 6 samples total: 3 control and 3 treatment replicates
- Reference genome in `data/reference/genome.fa`
- Gene annotation in `data/reference/genes.gtf`
- Sample metadata in `config/samples.tsv`

### Obtaining Data
You have several options for obtaining data:

1. **Use your own RNA-Seq data** - If you have experimental data, organize it according to the structure above

2. **Download public data from SRA** - Use the SRA toolkit to download datasets:
   ```bash
   # Example: Download a small dataset
   prefetch SRR_ID
   fasterq-dump --split-files --gzip SRR_ID
   ```

3. **Create test data** - For testing your workflow, create small subsampled files:
   ```bash
   # Subsample existing FASTQ files (first 100,000 reads)
   zcat original_R1.fastq.gz | head -400000 | gzip > control_1_R1.fastq.gz
   ```

4. **Download reference files** - Get reference genome and annotation from Ensembl:
   ```bash
   # Example for human chromosome 22 (for testing)
   wget http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
   wget http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
   ```

**Note:** Start with a small test dataset (e.g., one chromosome) to verify your workflow works before running on full-scale data.

## Assignment Tasks

### Task 1: Quality Control and Trimming
Implement Snakemake rules to:
1. Run **FastQC** on raw FASTQ files (both R1 and R2)
2. Trim adapters and low-quality bases using **Trim Galore**
3. Run **FastQC** on trimmed reads
4. Aggregate all QC reports using **MultiQC**

#### Expected outputs:
- `results/qc/raw_fastqc/{sample}_R{read}_fastqc.html`
- `results/trimmed/{sample}_R{read}_trimmed.fastq.gz`
- `results/qc/trimmed_fastqc/{sample}_R{read}_trimmed_fastqc.html`
- `results/qc/multiqc_report.html`

### Task 2: Read Alignment
Implement Snakemake rules to:
1. Index the reference genome using **STAR**
2. Align trimmed reads to the reference genome
3. Sort BAM files using **SAMtools**
4. Index the sorted BAM files
5. Generate alignment statistics

#### Expected outputs:
- `results/star_index/` (genome index files)
- `results/aligned/{sample}.sorted.bam`
- `results/aligned/{sample}.sorted.bam.bai`
- `results/stats/{sample}.alignment_stats.txt`

### Task 3: Read Quantification
Implement Snakemake rules to:
1. Count reads per gene using **featureCounts**
2. Save individual sample counts to separate files

#### Expected outputs:
- `results/counts/{sample}.counts.txt`
- `results/counts/{sample}.counts.txt.summary`

### Task 4: Count Matrix Generation
Use the provided Python script (`workflow/scripts/aggregate_counts.py`) and create a Snakemake rule to:
1. Read all individual count files
2. Combine them into a single count matrix
3. Save the matrix with genes as rows and samples as columns

**Note:** The Python script is provided for you - focus on writing the Snakemake rule.

#### Expected outputs:
- `results/counts/count_matrix.tsv`

### Task 5: Exploratory Analysis
Use the provided R scripts and create Snakemake rules to:
1. Generate a PCA plot showing sample clustering (use `workflow/scripts/plot_pca.R`)
2. Create a sample correlation heatmap (use `workflow/scripts/plot_correlation.R`)
3. Pass the count matrix and sample metadata as inputs

**Note:** The R scripts are provided for you - focus on writing the Snakemake rules.

#### Expected outputs:
- `results/analysis/pca_plot.pdf`
- `results/analysis/correlation_heatmap.pdf`

Use the provided R script and create a Snakemake rule to:
1. Load the count matrix and sample metadata
2. Create a DESeq2 dataset object (use `workflow/scripts/prepare_deseq2.R`)
3. Calculate size factors for normalization
4. Save the DESeq2 object and normalized counts
5. Calculate size factors for normalization
6. Save the DESeq2 object and normalized counts

**Note:** The R script is provided for you - focus on writing the Snakemake rule.

#### Expected outputs:
- `results/deseq2/deseq2_dataset.rds`
- `results/deseq2/normalized_counts.tsv`

### Task 7: Workflow Configuration
Create a proper Snakemake workflow structure:
1. Use a `Snakefile` with clear rule definitions
2. Create a `config.yaml` file for configurable parameters
3. Use wildcards to handle multiple samples
4. Define an `all` rule that specifies all final outputs
5. Add comments to your Snakefile explaining each rule

### Task 8: Documentation
Document your workflow:
1. Create a brief report (`report.md`) with 2-3 paragraphs summarizing:
   - Overall read quality and trimming results
   - Alignment rates across samples
   - General patterns observed in PCA
   - Any quality concerns or interesting observations
2. Generate a workflow diagram (DAG)

## Workflow Structure
Your final directory structure should look like:
```
.
├── workflow/
│   ├── Snakefile
│   └── scripts/
│       ├── aggregate_counts.py
│       ├── plot_pca.R
│       ├── plot_correlation.R
│       └── prepare_deseq2.R
├── config/
│   ├── config.yaml
├── data/
│   ├── fastq/
│   │   ├── control_1_R1.fastq.gz
│   │   ├── control_1_R2.fastq.gz
│   │   └── ...
│   └── reference/
│       ├── genome.fa
│       └── genes.gtf
├── results/
│   ├── qc/
│   ├── trimmed/
│   ├── aligned/
│   ├── counts/
│   ├── analysis/
│   └── deseq2/
├── logs/
├── report.md
└── dag.png
```

## Running Your Workflow
Execute your workflow with:
```bash
# Dry-run to check workflow
snakemake --dry-run

# Run the full workflow
snakemake --cores 8

# Generate workflow diagram
snakemake --dag | dot -Tpng > dag.png
```

## Deliverables
Submit the following via GitHub by merging your working branch into the main branch via a pull request:
1. `workflow/Snakefile` - Your complete workflow
2. `config/config.yaml` - Configuration file
3. `report.md` - Brief analysis report
4. `dag.png` - DAG visualization of your workflow
5. `results/analysis/pca_plot.pdf` - PCA plot
6. `results/analysis/correlation_heatmap.pdf` - Sample correlation heatmap

**Do not commit large data files (FASTQ, BAM files) to GitHub!** Use `.gitignore` appropriately.

## Tips
- Start small: test your workflow one rule at a time
- Use `snakemake --dry-run` frequently to check for errors
- All analysis scripts (Python and R) are provided in `workflow/scripts/` - you just need to call them from Snakemake. Use the `script:` directive to run these scripts rather than the `shell:` directive.
- Use the `log:` directive in rules to capture error messages
- Make sure to add results and .snakemake directories to your `.gitignore` file

## Resources
- [Snakemake Documentation](https://snakemake.readthedocs.io/)
- [STAR Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
- [featureCounts Documentation](https://subread.sourceforge.net/featureCounts.html)
- [DESeq2 Vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [Trim Galore User Guide](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)

## Due Date
**Tuesday, March 10th at Noon**
