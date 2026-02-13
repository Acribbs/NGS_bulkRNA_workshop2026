#!/usr/bin/env bash
# ============================================================================
# download_data.sh
# Downloads example FASTQ files and human reference transcriptome for the
# cribbslab pipeline_pseudobulk tutorial.
#
# Dataset: Griffith Lab RNA-seq tutorial data (UHR vs HBR, chr22-subsampled)
#   - Universal Human Reference (UHR): pooled RNA from 10 cancer cell lines
#   - Human Brain Reference (HBR): pooled RNA from multiple brain donors
#   - 3 biological replicates per condition, paired-end 101 bp reads
#   - Pre-filtered to chromosome 22 for fast processing (~200 MB total)
#
# Reference: Ensembl release 115, Homo sapiens GRCh38 cDNA transcriptome
#
# Usage:
#   chmod +x download_data.sh
#   ./download_data.sh
# ============================================================================

set -euo pipefail

WORK_DIR="$(pwd)"
DATA_DIR="${WORK_DIR}/data"
REF_DIR="${WORK_DIR}/reference"

echo "============================================================"
echo " Bulk RNA-seq Workshop: Data Download Script"
echo "============================================================"
echo ""

# ── 1. Create directory structure ──────────────────────────────────────────
echo "[1/4] Creating directory structure..."
mkdir -p "${DATA_DIR}"
mkdir -p "${REF_DIR}"

# ── 2. Download example FASTQ files ───────────────────────────────────────
echo "[2/4] Downloading example FASTQ dataset (chr22-subsampled, ~200 MB)..."
echo "       Source: Griffith Lab RNA-seq tutorial (UHR vs HBR)"
echo ""

cd "${DATA_DIR}"

if [ -f "HBR_UHR_ERCC_ds_5pc.tar" ] || [ -d "UHR_Rep1" ] || [ -f "UHR_Rep1.fastq.1.gz" ]; then
    echo "       FASTQ data appears to already exist. Skipping download."
else
    wget -q --show-progress \
        http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar \
        -O HBR_UHR_ERCC_ds_5pc.tar

    echo "       Extracting archive..."
    tar -xf HBR_UHR_ERCC_ds_5pc.tar
    rm -f HBR_UHR_ERCC_ds_5pc.tar
fi

echo "       Done."
echo ""

# ── 3. Rename FASTQ files to pipeline-compatible format ───────────────────
# The pipeline expects: sampleName.fastq.1.gz and sampleName.fastq.2.gz
echo "[3/4] Renaming FASTQ files for pipeline_pseudobulk compatibility..."
echo "       Expected format: <sample>.fastq.1.gz / <sample>.fastq.2.gz"
echo ""

cd "${DATA_DIR}"

# The Griffith Lab archive typically contains files in the format:
#   UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz
#   UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz
# We need to rename them to:
#   UHR_Rep1.fastq.1.gz and UHR_Rep1.fastq.2.gz

shopt -s nullglob
for f in *read1.fastq.gz; do
    if [ -f "$f" ]; then
        # Extract the sample name (e.g., UHR_Rep1 from UHR_Rep1_ERCC-Mix1_...)
        sample=$(echo "$f" | sed -E 's/^(UHR_Rep[0-9]+|HBR_Rep[0-9]+)_.*/\1/')
        new_name="${sample}.fastq.1.gz"
        if [ ! -f "$new_name" ]; then
            mv "$f" "$new_name"
            echo "       Renamed: $f -> $new_name"
        fi
    fi
done

for f in *read2.fastq.gz; do
    if [ -f "$f" ]; then
        sample=$(echo "$f" | sed -E 's/^(UHR_Rep[0-9]+|HBR_Rep[0-9]+)_.*/\1/')
        new_name="${sample}.fastq.2.gz"
        if [ ! -f "$new_name" ]; then
            mv "$f" "$new_name"
            echo "       Renamed: $f -> $new_name"
        fi
    fi
done
shopt -u nullglob

echo "       Done."
echo ""

# ── 4. Download human reference transcriptome ─────────────────────────────
echo "[4/4] Downloading Ensembl human reference transcriptome (cDNA)..."
echo "       Release: 115 | Assembly: GRCh38 | Size: ~100 MB"
echo ""

cd "${REF_DIR}"

TRANSCRIPTOME="Homo_sapiens.GRCh38.cdna.all.fa.gz"

if [ -f "${TRANSCRIPTOME}" ]; then
    echo "       Reference transcriptome already exists. Skipping download."
else
    wget -q --show-progress \
        "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/${TRANSCRIPTOME}"
fi

echo "       Done."
echo ""

# ── Summary ────────────────────────────────────────────────────────────────
echo "============================================================"
echo " Download Complete!"
echo "============================================================"
echo ""
echo " FASTQ files:    ${DATA_DIR}/"

cd "${DATA_DIR}"
echo ""
echo "   Condition 1 - UHR (Universal Human Reference):"
shopt -s nullglob
for f in UHR_*.fastq.*.gz; do
    echo "     $f"
done
echo ""
echo "   Condition 2 - HBR (Human Brain Reference):"
for f in HBR_*.fastq.*.gz; do
    echo "     $f"
done
shopt -u nullglob

echo ""
echo " Reference:       ${REF_DIR}/${TRANSCRIPTOME}"
echo ""
echo " Next step: Follow the README.md to configure and run the pipeline."
echo "============================================================"
