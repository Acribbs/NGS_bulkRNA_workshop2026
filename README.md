# NGS Bulk RNA-seq Workshop 2026

A comprehensive workshop covering the complete workflow for bulk RNA-seq analysis, from raw sequencing data processing to differential expression analysis.

---

## Workshop Overview

This workshop is divided into **two main stages**:

1. **Stage 1: Upstream Pipeline** — Process raw FASTQ files to generate transcript abundance estimates
2. **Stage 2: Downstream Analysis** — Perform differential expression analysis using DESeq2

**Prerequisites:** Basic familiarity with command-line tools and R is recommended.

---

## Stage 1: Upstream Pipeline Tutorial

**Location:** `Upstream_pipeline/` folder

The first stage of this workshop focuses on processing raw RNA-seq data from FASTQ files to transcript-level abundance estimates using Kallisto pseudoalignment. This is a critical preprocessing step that must be completed before performing downstream differential expression analysis.

### What You'll Learn

- How to set up and run the cribbslab `pipeline_pseudobulk` pipeline
- Quality control assessment using FastQC and MultiQC
- Building Kallisto indices from reference transcriptomes
- Quantifying transcript abundances using Kallisto
- Interpreting pipeline outputs and QC reports

### Quick Start

1. **Navigate to the upstream pipeline directory:**
   ```bash
   cd Upstream_pipeline
   ```

2. **Follow the detailed tutorial:**
   ```bash
   # Open the README for complete instructions
   cat README.md
   # Or view it in your preferred text editor
   ```

3. **Key steps:**
   - Install the conda environment (see `Upstream_pipeline/README.md` for details)
   - Download the example dataset
   - Configure and run the pipeline
   - Examine the outputs (especially the `abundance.tsv` files in the `quant/` directory)

### Expected Outputs

After completing Stage 1, you should have:
- Quality control reports (`reports/multiqc_report.html`)
- Kallisto index (`kalindex/`)
- Per-sample transcript abundance estimates (`quant/<sample>/abundance.tsv`)
- Kallisto bootstrap files (`quant/<sample>/abundance.h5`) for downstream analysis

**Important:** The upstream pipeline must be completed successfully before proceeding to Stage 2. The downstream analysis workflows use pre-processed count data, but understanding the upstream pipeline is essential for interpreting results.

**For detailed instructions, see:** [`Upstream_pipeline/README.md`](Upstream_pipeline/README.md)

---

## Stage 2: Downstream Analysis Workflows

**Location:** Root directory

After completing the upstream pipeline, you can proceed to differential expression analysis using DESeq2. This workshop includes two comprehensive tutorials covering different statistical approaches:

### Tutorial 1: DESeq2 Wald Test Analysis

**File:** `deseq2_wald_test_tutorial.Rmd`

**Use Case:** Pairwise comparisons between two conditions (e.g., treated vs. untreated)

**Dataset:** Airway smooth muscle cells treated with dexamethasone
- **Experimental Design:** Paired design with 4 cell lines × 2 treatments = 8 samples
- **Biological Question:** Which genes are differentially expressed following dexamethasone treatment?

**When to Use Wald Tests:**
- Comparing two conditions directly
- Testing specific contrasts or coefficients
- When you have a clear reference condition
- Simple experimental designs with one or few comparisons

**Key Features:**
- Log2 fold change estimates
- Standard errors for fold changes
- P-values for each gene
- LFC shrinkage for improved estimates
- Visualization with volcano plots, heatmaps, and PCA

#### How to Run the Wald Test Tutorial

1. **Install Required R Packages:**

   Open R or RStudio and install the necessary packages:
   ```r
   # Install Bioconductor packages
   if (!require("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   
   BiocManager::install(c("DESeq2", "airway", "org.Hs.eg.db", "EnhancedVolcano"))
   
   # Install CRAN packages
   install.packages(c("ggplot2", "pheatmap", "dplyr", "tidyr", "gridExtra"))
   ```

2. **Prepare the Data:**

   The tutorial uses the `airway` dataset from Bioconductor. To prepare the data files:
   ```r
   # Run the data preparation script
   source("Data_download.R")
   ```
   
   This will create three files in your working directory:
   - `airway_metadata_clean.txt` — Sample metadata with experimental design
   - `airway_counts_ensembl.txt` — Gene expression counts (Ensembl IDs)
   - `airway_gene_annotation.txt` — Gene symbols and names

3. **Run the Tutorial:**

   Open `deseq2_wald_test_tutorial.Rmd` in RStudio and:
   - Click "Knit" to render the HTML report, OR
   - Run chunks interactively by clicking "Run" on each code chunk
   
   Alternatively, render from the command line:
   ```bash
   Rscript -e "rmarkdown::render('deseq2_wald_test_tutorial.Rmd')"
   ```

4. **Expected Workflow:**
   - Load libraries and data files
   - Create DESeq2 object with design formula `~ CellLine + Treatment`
   - Estimate size factors and dispersions
   - Perform Wald test for differential expression
   - Apply log2 fold change shrinkage
   - Generate visualizations (PCA, volcano plots, heatmaps)
   - Export results tables

---

### Tutorial 2: DESeq2 LRT (Likelihood Ratio Test) Analysis

**File:** `deseq2_lrt_test_tutorial.Rmd`

**Use Case:** Testing for changes across multiple conditions or timepoints

**Dataset:** Mouse mammary gland development
- **Experimental Design:** Factorial design with 2 cell types × 3 developmental stages × 2 replicates = 12 samples
- **Biological Question:** Which genes show temporal changes during mammary gland development, and do these changes differ between cell types?

**When to Use LRT Tests:**
- Testing for ANY change across multiple conditions
- Comparing nested models (full vs. reduced)
- Identifying genes with complex expression patterns
- Testing interaction effects (e.g., different temporal patterns in different cell types)

**Key Features:**
- Tests whether genes change at all (not specific comparisons)
- Model comparison approach (full model vs. reduced model)
- Interaction term testing (`CellType:Status`)
- Often followed by pairwise Wald tests for specific contrasts
- Advanced visualizations with ComplexHeatmap

#### How to Run the LRT Test Tutorial

1. **Install Required R Packages:**

   ```r
   # Install Bioconductor packages
   if (!require("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   
   BiocManager::install(c("DESeq2", "GEOquery", "org.Mm.eg.db", "DEGreport", "ComplexHeatmap"))
   
   # Install CRAN packages
   install.packages(c("ggplot2", "pheatmap", "dplyr", "tidyr", "tibble", "circlize"))
   ```

2. **Prepare the Data:**

   The tutorial uses data from GEO (GSE60450). To prepare the data files:
   ```r
   # Run the data preparation script
   source("Data_download.R")
   ```
   
   **Note:** The `Data_download.R` script contains code for both datasets. The GSE60450 section will:
   - Download data from GEO
   - Convert Entrez IDs to Ensembl IDs
   - Handle duplicate Ensembl IDs by aggregating counts
   - Create cleaned metadata and count files
   
   This will create three files in your working directory:
   - `GSE60450_metadata_clean.txt` — Sample metadata with experimental design
   - `GSE60450_counts_ensembl.txt` — Gene expression counts (Ensembl IDs)
   - `GSE60450_gene_annotation.txt` — Gene symbols and names

3. **Run the Tutorial:**

   Open `deseq2_lrt_test_tutorial.Rmd` in RStudio and:
   - Click "Knit" to render the HTML report, OR
   - Run chunks interactively by clicking "Run" on each code chunk
   
   Alternatively, render from the command line:
   ```bash
   Rscript -e "rmarkdown::render('deseq2_lrt_test_tutorial.Rmd')"
   ```

4. **Expected Workflow:**
   - Load libraries and data files
   - Create DESeq2 object with design formula `~ CellType + Status + CellType:Status`
   - Estimate size factors and dispersions
   - Perform LRT test (comparing full model vs. reduced model)
   - Identify genes with significant changes
   - Follow-up with pairwise Wald tests for specific contrasts
   - Generate visualizations (PCA, heatmaps, expression plots)
   - Export results tables

---

## Complete Workflow Summary

### Recommended Learning Path

1. **Start with Stage 1 (Upstream Pipeline):**
   - Complete the Kallisto pipeline tutorial
   - Understand quality control metrics
   - Familiarize yourself with transcript abundance outputs

2. **Proceed to Stage 2 (Downstream Analysis):**
   - Begin with the **Wald Test Tutorial** (simpler, pairwise comparison)
   - Then try the **LRT Test Tutorial** (more complex, multiple conditions)
   - Compare the approaches and understand when to use each

### Data Flow

```
Raw FASTQ files
    ↓
[Stage 1: Upstream Pipeline]
    ↓
Transcript abundance estimates (abundance.tsv)
    ↓
[Stage 2: Downstream Analysis]
    ↓
Differential expression results
```

**Note:** The downstream tutorials use pre-processed example datasets (airway and GSE60450) for teaching purposes. In a real analysis, you would use the count data generated from Stage 1.

---

## Troubleshooting

### Common Issues

**R Package Installation Errors:**
- Ensure you have the latest version of R (≥ 4.0 recommended)
- For Bioconductor packages, use `BiocManager::install()` rather than `install.packages()`
- Some packages may require system dependencies (e.g., `libcurl`, `libxml2`)

**Data File Not Found Errors:**
- Ensure you've run `source("Data_download.R")` before running the tutorials
- Check that data files are in the same directory as the R Markdown files
- Verify file names match exactly (case-sensitive)

**Memory Issues:**
- DESeq2 analysis can be memory-intensive for large datasets
- Close other applications if running on a memory-constrained machine
- Consider filtering low-count genes before analysis

**Rendering R Markdown Errors:**
- Ensure all required packages are installed
- Check that all data files are present
- Review error messages in the R console for specific issues

---

## Additional Resources

### Documentation
- **Upstream Pipeline:** See [`Upstream_pipeline/README.md`](Upstream_pipeline/README.md) for detailed instructions
- **DESeq2 Manual:** https://bioconductor.org/packages/DESeq2/
- **Kallisto Manual:** https://pachterlab.github.io/kallisto/manual

### Key Papers
- **Kallisto:** Bray, N. L. et al. "Near-optimal probabilistic RNA-seq quantification." *Nature Biotechnology* 34, 525–527 (2016)
- **DESeq2:** Love, M. I. et al. "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." *Genome Biology* 15, 550 (2014)

### Workshop Materials
- Presentation slides are included in the repository:
  - `NGS_experimental_design.pptx` — Experimental design principles
  - `NGS_data_analysis_20210304.pptx` — Data analysis overview
  - `Deseq2_modelling.pptx` — DESeq2 statistical modeling
  - `SequencingMethods.pptx` — Sequencing technologies overview

---

## Workshop Structure

```
NGS_bulkRNA_workshop2026/
├── README.md                          # This file
├── Upstream_pipeline/                 # Stage 1: Data processing
│   ├── README.md                     # Detailed pipeline instructions
│   └── download_data.sh              # Data download script
├── Data_download.R                   # Data preparation for downstream tutorials
├── deseq2_wald_test_tutorial.Rmd     # Stage 2: Wald test tutorial
├── deseq2_lrt_test_tutorial.Rmd      # Stage 2: LRT test tutorial
└── [Presentation slides and other materials]
```

---

## Getting Help

If you encounter issues:
1. Check the troubleshooting section above
2. Review the detailed README in `Upstream_pipeline/`
3. Consult the R Markdown files for inline comments and explanations
4. Refer to the official documentation for DESeq2 and Kallisto

---

**Workshop Date:** 2026  
**Last Updated:** February 2026
