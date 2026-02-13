# Bulk RNA-seq Workshop: Pseudoalignment with Kallisto

A hands-on tutorial for running the **cribbslab `pipeline_pseudobulk`** pipeline on a small example human RNA-seq dataset. This pipeline uses [Kallisto](https://pachterlab.github.io/kallisto/) for transcript-level pseudoalignment and quantification, coordinated by the [CGAT-core / Ruffus](https://github.com/cgat-developers/cgat-core) workflow framework.

---

## Overview

This tutorial walks through every step of a bulk RNA-seq pseudoalignment analysis, from downloading raw data to interpreting output files. By the end you will have produced per-sample transcript abundance estimates (`abundance.tsv`) and quality-control reports via FastQC and MultiQC.

**Pipeline steps:**

1. **FastQC** — per-sample read quality assessment
2. **MultiQC** — aggregate QC report across all samples
3. **Kallisto index** — build a pseudoalignment index from the reference transcriptome
4. **Kallisto quant** — pseudoalign reads and estimate transcript abundances
5. **Kallisto MultiQC** — summarise alignment statistics across samples

**Estimated runtime:** ~15–30 minutes on a standard laptop (depending on CPU cores and download speed).

---

## Example Dataset

This tutorial uses the **Griffith Lab RNA-seq tutorial dataset**, a widely used teaching resource in bioinformatics. The data compares two commercially available reference RNA samples:

| Condition | Description | Replicates |
|-----------|-------------|------------|
| **UHR** | Universal Human Reference — pooled total RNA from 10 human cancer cell lines (Agilent) | 3 |
| **HBR** | Human Brain Reference — pooled total RNA from several adult human brains (Ambion) | 3 |

**Key properties:**

- **Organism:** *Homo sapiens*
- **Sequencing:** Illumina HiSeq, paired-end, 101 bp reads
- **Library prep:** TruSeq Stranded Total RNA with Ribo-Zero Gold
- **Subsampled:** Pre-filtered for reads mapping to **chromosome 22 only**, keeping the total download size around ~200 MB. This makes the dataset fast to process while still being biologically meaningful.
- **Spike-ins:** ERCC spike-in controls (Mix 1 for UHR, Mix 2 for HBR)

After downloading and renaming, you will have 12 FASTQ files (6 samples × 2 paired-end reads):

```
UHR_Rep1.fastq.1.gz    UHR_Rep1.fastq.2.gz
UHR_Rep2.fastq.1.gz    UHR_Rep2.fastq.2.gz
UHR_Rep3.fastq.1.gz    UHR_Rep3.fastq.2.gz
HBR_Rep1.fastq.1.gz    HBR_Rep1.fastq.2.gz
HBR_Rep2.fastq.1.gz    HBR_Rep2.fastq.2.gz
HBR_Rep3.fastq.1.gz    HBR_Rep3.fastq.2.gz
```

**Source:** [Griffith Lab RNA-seq tutorial](https://rnabio.org/) — original data from Illumina BaseSpace.

---

## Prerequisites

### Software

You will need the following tools installed. The recommended approach is to use [Miniconda](https://docs.conda.io/en/latest/miniconda.html) with the `mamba` package manager.

- **Python ≥ 3.6**
- **cribbslab** (the pipeline package)
- **cgat-core** (workflow engine)
- **ruffus** (pipeline framework)
- **Kallisto** (v0.48.0 recommended)
- **FastQC**
- **MultiQC**
- **wget** (for downloading data)

### Install the conda environment

```bash
# Install mamba if you haven't already
conda install -n base -c conda-forge mamba

# Clone the cribbslab repository
git clone https://github.com/cribbslab/cribbslab.git
cd cribbslab

# Create the environment from the provided YAML
mamba env create -f conda/environments/cribbslab.yml

# Activate the environment
conda activate cribbslab

# Install cribbslab in development mode
pip install -e .
```

Verify the installation:

```bash
cribbslab --help
```

You should see a list of available pipelines, including `pseudobulk`.

---

## Step-by-Step Tutorial

### Step 1: Create a working directory

Choose a location on your machine for this tutorial. All pipeline commands should be run from within this directory.

```bash
mkdir -p ~/bulkRNA_workshop
cd ~/bulkRNA_workshop
```

### Step 2: Download the data

Copy the provided download script into your working directory and run it. The script will download the example FASTQ files, rename them to the format expected by the pipeline, and download the Ensembl human reference transcriptome.

```bash
# Copy the download script (adjust path as needed)
cp /path/to/download_data.sh .
chmod +x download_data.sh
./download_data.sh
```

Alternatively, you can download everything manually:

```bash
# Create directories
mkdir -p data reference

# Download and extract FASTQ files
cd data
wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
tar -xf HBR_UHR_ERCC_ds_5pc.tar
rm HBR_UHR_ERCC_ds_5pc.tar

# Rename files to pipeline format: <sample>.fastq.<read>.gz
# The exact filenames in the archive may vary; adapt as needed:
for f in *read1.fastq.gz; do
    sample=$(echo "$f" | sed -E 's/^(UHR_Rep[0-9]+|HBR_Rep[0-9]+)_.*/\1/')
    mv "$f" "${sample}.fastq.1.gz"
done
for f in *read2.fastq.gz; do
    sample=$(echo "$f" | sed -E 's/^(UHR_Rep[0-9]+|HBR_Rep[0-9]+)_.*/\1/')
    mv "$f" "${sample}.fastq.2.gz"
done

cd ..

# Download the human reference transcriptome (Ensembl release 115, GRCh38)
cd reference
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
cd ..
```

After this step your directory should look like:

```
~/bulkRNA_workshop/
├── data/
│   ├── UHR_Rep1.fastq.1.gz
│   ├── UHR_Rep1.fastq.2.gz
│   ├── UHR_Rep2.fastq.1.gz
│   ├── UHR_Rep2.fastq.2.gz
│   ├── UHR_Rep3.fastq.1.gz
│   ├── UHR_Rep3.fastq.2.gz
│   ├── HBR_Rep1.fastq.1.gz
│   ├── HBR_Rep1.fastq.2.gz
│   ├── HBR_Rep2.fastq.1.gz
│   ├── HBR_Rep2.fastq.2.gz
│   ├── HBR_Rep3.fastq.1.gz
│   └── HBR_Rep3.fastq.2.gz
└── reference/
    └── Homo_sapiens.GRCh38.cdna.all.fa.gz
```

### Step 3: Set up the pipeline working directory

The pipeline expects FASTQ files and the reference transcriptome to be accessible from the directory where you run it. Copy (or symlink) files into a single run directory:

```bash
mkdir -p run_pseudobulk
cd run_pseudobulk

# Symlink FASTQ files into the run directory
ln -s ../data/*.fastq.*.gz .

# Symlink the reference transcriptome
ln -s ../reference/Homo_sapiens.GRCh38.cdna.all.fa.gz .
```

### Step 4: Generate the pipeline configuration

```bash
cribbslab pseudobulk config
```

This creates a `pipeline.yml` file in your current directory with default settings. Open it in a text editor to customise for this dataset.

### Step 5: Edit the configuration file

Open `pipeline.yml` and update the following parameters:

```yaml
# Point to the human reference transcriptome
cdna_fasta: Homo_sapiens.GRCh38.cdna.all.fa.gz

# Our data is paired-end, so keep this as 0
kallisto_single: 0

kal_quant:
    options: --rf-stranded   # Strand-specific: first read reverse
    threads: 4               # Adjust based on your CPU cores
    bootstraps: 100          # Number of bootstrap samples for uncertainty estimation
    pseudobam: False         # Set to True if you want pseudoBAM output
```

**Configuration notes:**

- `cdna_fasta`: Must match the exact filename of the reference transcriptome FASTA in your run directory.
- `kallisto_single`: Set to `0` for paired-end data, `1` for single-end.
- `kal_quant.options`: Use `--rf-stranded` for TruSeq Stranded libraries (first read maps to the reverse strand). Use `--fr-stranded` if the first read maps forward. Omit entirely for unstranded libraries.
- `kal_quant.threads`: Number of CPU threads per sample. 4–8 is typical.
- `kal_quant.bootstraps`: Controls the number of bootstrap resamples Kallisto uses to estimate technical variance. 100 is the standard default.

### Step 6: Run the pipeline

To run all steps locally:

```bash
cribbslab pseudobulk make full -v5 --local
```

The `-v5` flag sets verbosity to maximum, so you can see detailed progress. The `--local` flag runs all tasks on your local machine (as opposed to submitting to a cluster scheduler).

You can also run individual tasks:

```bash
# Run only FastQC
cribbslab pseudobulk make fastqc -v5 --local

# Run only MultiQC (after FastQC has completed)
cribbslab pseudobulk make multiqc -v5 --local

# Run only the Kallisto index (requires reference transcriptome)
cribbslab pseudobulk make kallisto_index -v5 --local

# Run only Kallisto quantification (after index is built)
cribbslab pseudobulk make kal_quant -v5 --local
```

### Step 7: Examine the outputs

After a successful run, your directory will contain several new folders:

```
run_pseudobulk/
├── fastqc/                          # FastQC reports (per file)
│   ├── UHR_Rep1.fastq.1_fastqc.html
│   ├── UHR_Rep1.fastq.1_fastqc.zip
│   ├── ...
├── reports/
│   └── multiqc_report.html          # Aggregated QC report
├── kalindex/
│   └── Homo_sapiens.GRCh38.cdna.all.idx  # Kallisto index
├── quant/
│   ├── UHR_Rep1/
│   │   ├── abundance.tsv            # Main output: transcript counts
│   │   ├── abundance.h5             # HDF5 format with bootstraps
│   │   └── run_info.json            # Run metadata
│   ├── UHR_Rep2/
│   ├── UHR_Rep3/
│   ├── HBR_Rep1/
│   ├── HBR_Rep2/
│   ├── HBR_Rep3/
│   └── kallisto_multiqc.html        # Kallisto QC summary
└── pipeline.yml
```

### Understanding the key outputs

**`abundance.tsv`** — This is the primary output. Each sample gets its own file with four columns:

| Column | Description |
|--------|-------------|
| `target_id` | Ensembl transcript ID |
| `length` | Effective transcript length |
| `eff_length` | Fragment-length-corrected effective length |
| `est_counts` | Estimated number of reads from this transcript |
| `tpm` | Transcripts Per Million — normalised abundance |

Example (first few lines):

```
target_id                  length  eff_length  est_counts  tpm
ENST00000456328.2          1657    1470.03     12.5432     1.23456
ENST00000450305.2          632     445.03      0           0
...
```

**`reports/multiqc_report.html`** — Open this in a browser to see an interactive summary of read quality across all samples. Check for adapter contamination, GC bias, and per-base quality scores.

**`quant/kallisto_multiqc.html`** — Open this to review alignment rates and fragment length distributions. A typical pseudoalignment rate for chr22-subsampled data will be lower than a full-genome run, since many reads from other chromosomes were filtered out during subsampling.

---

## Pipeline Architecture

The pipeline is built on the [CGAT-core](https://github.com/cgat-developers/cgat-core) framework, which uses [Ruffus](http://www.ruffus.org.uk/) for task management. Each step is a Python function decorated with Ruffus decorators that define input/output relationships and execution dependencies.

```
Input FASTQs ──> fastqc ──> multiqc ──────────────────────> full
                                                              ^
Reference FA ──> kallisto_index ──> kal_quant ──> kallisto_multiqc
```

The `full` task acts as a sentinel that depends on all other tasks. When you run `make full`, Ruffus determines which tasks need to run based on file timestamps and executes them in the correct order.

Key design decisions in this pipeline:

- **Pseudoalignment over traditional alignment:** Kallisto uses an approach based on pseudoalignment to rapidly determine the compatibility of reads with transcripts, without performing a base-level alignment. This makes it significantly faster than STAR or HISAT2 while producing comparable quantification results for transcript-level analysis.
- **Bootstrap resampling:** The `-b 100` flag tells Kallisto to perform 100 bootstrap resamplings, which can later be used by tools like [Sleuth](https://pachterlab.github.io/sleuth/) to model technical variance during differential expression analysis.
- **Strand-specificity:** The `--rf-stranded` flag matches TruSeq Stranded library preparation, where the first read corresponds to the reverse complement of the RNA fragment.

---

## Troubleshooting

**"No files matching pattern *.fastq.*.gz"**
The pipeline expects files named `<sample>.fastq.1.gz` and `<sample>.fastq.2.gz` in the current working directory. Check that your FASTQ files follow this naming convention and that you are running the pipeline from the correct directory.

**Low pseudoalignment rates**
Since this dataset is subsampled to chromosome 22, many reads will not map to the full human transcriptome. This is expected for this tutorial dataset. A full-genome dataset would typically show 60–90% pseudoalignment rates.

**Kallisto index fails**
Ensure the `cdna_fasta` path in `pipeline.yml` matches the exact filename of the gzipped FASTA file in your working directory. The file must be readable and not corrupted (check with `zcat <file> | head`).

**MultiQC locale errors**
If you see encoding-related errors, the pipeline already handles this by exporting `LC_ALL=en_US.UTF-8` and `LANG=en_US.UTF-8`. If problems persist, ensure these locale settings are available on your system (`locale -a | grep en_US`).

**Memory issues**
The Kallisto index step for the full human transcriptome requires approximately 4 GB of RAM. If you are running on a memory-constrained machine, close other applications during the indexing step.

---

## Next Steps

After completing this tutorial, you might consider:

- **Differential expression analysis:** Use the `abundance.h5` files with [Sleuth](https://pachterlab.github.io/sleuth/) (R/Bioconductor) for transcript-level differential expression, or aggregate to gene level and use [DESeq2](https://bioconductor.org/packages/DESeq2/) or [edgeR](https://bioconductor.org/packages/edgeR/).
- **Running on your own data:** Replace the example FASTQ files with your own paired-end RNA-seq data, update the sample naming to match `<sample>.fastq.1.gz` / `<sample>.fastq.2.gz`, and adjust `pipeline.yml` accordingly.
- **Cluster execution:** Remove the `--local` flag and configure your CGAT-core settings (`.cgat.yml`) to submit jobs to a DRMAA-compatible cluster scheduler (SGE, SLURM, PBS).
- **Using a different organism:** Download the appropriate cDNA FASTA from [Ensembl FTP](https://www.ensembl.org/info/data/ftp/index.html) and update `cdna_fasta` in the configuration.

---

## References and Resources

- **Cribbslab pipelines:** https://github.com/cribbslab/cribbslab
- **CGAT-core documentation:** https://cgat-core.readthedocs.io/
- **Kallisto paper:** Bray, N. L. et al. "Near-optimal probabilistic RNA-seq quantification." *Nature Biotechnology* 34, 525–527 (2016). https://doi.org/10.1038/nbt.3519
- **Kallisto manual:** https://pachterlab.github.io/kallisto/manual
- **Griffith Lab RNA-seq tutorial:** https://rnabio.org/
- **Ensembl FTP downloads:** https://www.ensembl.org/info/data/ftp/index.html
- **FastQC documentation:** https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- **MultiQC documentation:** https://multiqc.info/
