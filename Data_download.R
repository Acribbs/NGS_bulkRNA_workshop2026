

# ============================================
# Airway Dataset Download and Preparation
# For DESeq2 Wald Test Tutorial
# ============================================

# Load required libraries
library(airway)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ============================================
# 1. LOAD AIRWAY DATA
# ============================================

cat("Loading airway dataset...\n")

data("airway")

# Extract counts and metadata
counts <- assay(airway)
metadata <- as.data.frame(colData(airway))

cat("Data loaded successfully!\n")

# ============================================
# 2. CLEAN METADATA
# ============================================

cat("\nCleaning metadata...\n")

metadata_clean <- data.frame(
  Sample = rownames(metadata),
  SampleName = metadata$Run,
  CellLine = metadata$cell,
  Treatment = metadata$dex,
  stringsAsFactors = FALSE
)

# Set factors with untreated as reference
metadata_clean$Treatment <- factor(metadata_clean$Treatment,
                                   levels = c("untrt", "trt"))
metadata_clean$CellLine <- factor(metadata_clean$CellLine)

rownames(metadata_clean) <- metadata_clean$SampleName

# ============================================
# 3. REORDER COUNTS TO MATCH METADATA
# ============================================

cat("\nReordering counts matrix...\n")

colnames(counts) <- metadata_clean$SampleName
counts_ordered <- counts[, rownames(metadata_clean)]

# ============================================
# 4. ADD GENE SYMBOLS
# ============================================

cat("\nAdding gene symbols...\n")

ensembl_ids <- rownames(counts_ordered)

gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

gene_names <- mapIds(org.Hs.eg.db,
                     keys = ensembl_ids,
                     column = "GENENAME",
                     keytype = "ENSEMBL",
                     multiVals = "first")

entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = ensembl_ids,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")

gene_annotation <- data.frame(
  EnsemblID = ensembl_ids,
  GeneSymbol = gene_symbols,
  GeneName = gene_names,
  EntrezID = entrez_ids,
  stringsAsFactors = FALSE
)

# ============================================
# 5. CREATE OUTPUT FILES
# ============================================

cat("\nCreating output files...\n")

counts_with_id <- data.frame(
  EnsemblID = rownames(counts_ordered),
  GeneSymbol = gene_annotation$GeneSymbol,
  counts_ordered,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ============================================
# 6. SAVE FILES
# ============================================

cat("\nSaving files...\n")

write.table(metadata_clean,
            file = "airway_metadata_clean.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

write.table(counts_with_id,
            file = "airway_counts_ensembl.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

write.table(gene_annotation,
            file = "airway_gene_annotation.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)





# ============================================
# GSE60450 Data Download and Preparation
# Converting to Ensembl IDs using org.Mm.eg.db
# Handling duplicate Ensembl IDs
# ============================================

# Load required libraries
library(GEOquery)
library(dplyr)
library(stringr)
library(org.Mm.eg.db)
library(AnnotationDbi)

# ============================================
# 1. DOWNLOAD DATA FROM GEO
# ============================================

cat("Downloading data from GEO...\n")

# Download metadata
gse <- getGEO("GSE60450", GSEMatrix = TRUE)
gse <- gse[[1]]
sample_info <- pData(gse)

# Download counts matrix
FileURL <- paste(
  "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE60450",
  "format=file",
  "file=GSE60450_Lactation-GenewiseCounts.txt.gz",
  sep = "&"
)
download.file(FileURL, "GSE60450_Lactation-GenewiseCounts.txt.gz")
counts <- read.delim("GSE60450_Lactation-GenewiseCounts.txt.gz", 
                     row.names = "EntrezGeneID")

cat("Data downloaded successfully!\n")

# ============================================
# 2. CLEAN METADATA
# ============================================

cat("\nCleaning metadata...\n")

metadata <- data.frame(
  Sample = sample_info$geo_accession,
  SampleName = str_extract(sample_info$description.1, "MCL1-[A-Z]+"),
  CellType = ifelse(grepl("luminal", sample_info$`immunophenotype:ch1`), 
                    "luminal", "basal"),
  Status = case_when(
    sample_info$`developmental stage:ch1` == "virgin" ~ "virgin",
    sample_info$`developmental stage:ch1` == "18.5 day pregnancy" ~ "pregnant",
    sample_info$`developmental stage:ch1` == "2 day lactation" ~ "lactate"
  ),
  stringsAsFactors = FALSE
)

# Add combined group identifier
metadata$Group <- paste(metadata$CellType, metadata$Status, sep = ".")

# Set factor levels (virgin and basal as reference levels)
metadata$Status <- factor(metadata$Status, 
                          levels = c("virgin", "pregnant", "lactate"))
metadata$CellType <- factor(metadata$CellType, 
                            levels = c("basal", "luminal"))

# Set row names
rownames(metadata) <- metadata$SampleName

cat("Metadata cleaned!\n")

# ============================================
# 3. CONVERT ENTREZ IDs TO ENSEMBL IDs
# ============================================

cat("\nConverting Entrez IDs to Ensembl IDs...\n")

entrez_ids <- rownames(counts)
cat(paste("Total genes with Entrez IDs:", length(entrez_ids), "\n"))

# Map Entrez to Ensembl
ensembl_ids <- mapIds(org.Mm.eg.db,
                      keys = entrez_ids,
                      column = "ENSEMBL",
                      keytype = "ENTREZID",
                      multiVals = "first")

# Map Entrez to Gene Symbols
gene_symbols <- mapIds(org.Mm.eg.db,
                       keys = entrez_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

# Map Entrez to Gene Names (full names)
gene_names <- mapIds(org.Mm.eg.db,
                     keys = entrez_ids,
                     column = "GENENAME",
                     keytype = "ENTREZID",
                     multiVals = "first")

# Create annotation data frame
gene_annotation <- data.frame(
  EntrezID = entrez_ids,
  EnsemblID = ensembl_ids,
  GeneSymbol = gene_symbols,
  GeneName = gene_names,
  Length = counts$Length,
  stringsAsFactors = FALSE
)

# Remove genes without Ensembl IDs
genes_with_ensembl <- !is.na(gene_annotation$EnsemblID)
cat(paste("Genes successfully mapped to Ensembl:", sum(genes_with_ensembl), "\n"))
cat(paste("Genes without Ensembl mapping (removed):", 
          sum(!genes_with_ensembl), "\n"))

gene_annotation <- gene_annotation[genes_with_ensembl, ]

# Check for duplicate Ensembl IDs
n_duplicates <- sum(duplicated(gene_annotation$EnsemblID))
cat(paste("Duplicate Ensembl IDs found:", n_duplicates, "\n"))

# ============================================
# 4. HANDLE DUPLICATE ENSEMBL IDs
# ============================================

cat("\nHandling duplicate Ensembl IDs by summing counts...\n")

# Filter counts to only genes with Ensembl IDs
counts_filtered <- counts[genes_with_ensembl, ]

# Remove Length column
counts_only <- counts_filtered[, -1]

# Shorten column names to match metadata
sample_columns <- colnames(counts_only)
new_sample_names <- str_replace(
  str_extract(sample_columns, "MCL1\\.[A-Z]+"),
  "\\.", "-"
)
colnames(counts_only) <- new_sample_names

cat("Sample column names after cleaning:\n")
print(colnames(counts_only))

# Add EnsemblID column for aggregation
counts_for_agg <- counts_only
counts_for_agg$EnsemblID <- gene_annotation$EnsemblID

# Aggregate using base R (more reliable with column names)
counts_aggregated <- aggregate(
  . ~ EnsemblID, 
  data = counts_for_agg, 
  FUN = sum
)

cat(paste("\nGenes after aggregating duplicates:", nrow(counts_aggregated), "\n"))

# Set row names
rownames(counts_aggregated) <- counts_aggregated$EnsemblID

# Remove EnsemblID column to get final counts matrix
counts_final <- counts_aggregated[, -1]

# Verify column names match
cat("\nVerifying column names match:\n")
cat("Metadata sample names:\n")
print(metadata$SampleName)
cat("\nCounts column names:\n")
print(colnames(counts_final))

all_match <- all(metadata$SampleName %in% colnames(counts_final))
cat(paste("\nAll sample names present:", all_match, "\n"))

# Reorder columns to match metadata row order
counts_ordered <- counts_final[, metadata$SampleName]

cat("Counts matrix successfully reordered!\n")

# ============================================
# 5. CREATE GENE ANNOTATION (UNIQUE ENSEMBL IDs)
# ============================================

cat("\nCreating unique gene annotation...\n")

# Use base R approach to avoid first() function issues
gene_annotation_unique <- aggregate(
  cbind(GeneSymbol, GeneName, Length) ~ EnsemblID,
  data = gene_annotation,
  FUN = function(x) {
    if(is.numeric(x)) {
      return(round(mean(x, na.rm = TRUE)))
    } else {
      return(x[1])  # Take first value for character columns
    }
  }
)

# For EntrezID, concatenate multiple IDs
entrez_concat <- aggregate(
  EntrezID ~ EnsemblID,
  data = gene_annotation,
  FUN = function(x) paste(unique(x), collapse = ";")
)

# Merge back
gene_annotation_unique <- merge(gene_annotation_unique, entrez_concat, by = "EnsemblID")

# Reorder columns
gene_annotation_unique <- gene_annotation_unique[, c("EnsemblID", "EntrezID", "GeneSymbol", "GeneName", "Length")]

# Convert factors to characters if needed
gene_annotation_unique$GeneSymbol <- as.character(gene_annotation_unique$GeneSymbol)
gene_annotation_unique$GeneName <- as.character(gene_annotation_unique$GeneName)

cat(paste("Unique genes in annotation:", nrow(gene_annotation_unique), "\n"))

# ============================================
# 6. CREATE FINAL OUTPUT FILES
# ============================================

cat("\nCreating final output files...\n")

# Counts with explicit EnsemblID and GeneSymbol columns
counts_with_id <- data.frame(
  EnsemblID = rownames(counts_ordered),
  GeneSymbol = gene_annotation_unique$GeneSymbol[
    match(rownames(counts_ordered), gene_annotation_unique$EnsemblID)
  ],
  counts_ordered,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ============================================
# 7. SAVE CLEAN FILES
# ============================================

cat("\nSaving cleaned files...\n")

write.table(metadata, 
            file = "GSE60450_metadata_clean.txt",
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)
cat("✓ Saved: GSE60450_metadata_clean.txt\n")

write.table(counts_with_id,
            file = "GSE60450_counts_ensembl.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
cat("✓ Saved: GSE60450_counts_ensembl.txt\n")

write.table(gene_annotation_unique,
            file = "GSE60450_gene_annotation.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
cat("✓ Saved: GSE60450_gene_annotation.txt\n")
