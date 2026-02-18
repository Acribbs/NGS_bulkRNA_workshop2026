# ===============================================================
# Package Setup
# Checks for and installs all required packages automatically.
# Works on Windows, macOS, and Linux.
# ===============================================================

# --- Check R version ---
# Bioconductor 3.22 requires R >= 4.5. Older R versions will fail
# to find packages, producing confusing "not available" errors.
if (getRversion() < "4.4") {
  stop(
    "R version ", getRversion(), " detected. ",
    "This workshop requires R >= 4.4.\n",
    "Please update R from https://cran.r-project.org/ and try again."
  )
}

# --- Ensure a writable library path exists ---
# On some systems (e.g. Windows with R in Program Files), the default
# library is not writable. Create a user library if needed.
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}

# --- Install BiocManager if missing ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# --- Define required packages ---
bioc_packages <- c(
  "airway", "org.Hs.eg.db", "AnnotationDbi",
  "GEOquery", "org.Mm.eg.db"
)
cran_packages <- c("dplyr", "stringr")

# --- Install missing packages ---
is_installed <- function(pkg) requireNamespace(pkg, quietly = TRUE)

missing_bioc <- bioc_packages[!vapply(bioc_packages, is_installed, logical(1))]
if (length(missing_bioc) > 0) {
  cat("Installing Bioconductor packages:", paste(missing_bioc, collapse = ", "), "\n")
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

missing_cran <- cran_packages[!vapply(cran_packages, is_installed, logical(1))]
if (length(missing_cran) > 0) {
  cat("Installing CRAN packages:", paste(missing_cran, collapse = ", "), "\n")
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

# --- Verify everything installed ---
all_packages  <- c(bioc_packages, cran_packages)
still_missing <- all_packages[!vapply(all_packages, is_installed, logical(1))]
if (length(still_missing) > 0) {
  stop(
    "Could not install: ", paste(still_missing, collapse = ", "), "\n",
    "Possible causes:\n",
    "  - No internet connection\n",
    "  - Library path not writable (try running R as administrator,\n",
    "    or set R_LIBS_USER to a writable directory)\n",
    "  - R version incompatible with current Bioconductor release\n",
    "Please install the packages above manually and re-run this script."
  )
}

cat("All required packages are installed.\n")

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
# 0. MEMORY WARNING AND SETUP
# ============================================

cat("\n============================================\n")
cat("GSE60450 Dataset Download\n")
cat("============================================\n")
cat("WARNING: This dataset is large and requires significant memory.\n")
cat("If you encounter memory errors, try:\n")
cat("1. Increase memory limit before starting R:\n")
cat("   ulimit -v 50000000  # In terminal before starting R\n")
cat("2. Or set in R before downloading:\n")
cat("   Sys.setenv(R_MAX_VSIZE = \"50Gb\")\n")
cat("============================================\n\n")

# Try to increase memory limit if possible (Mac/Linux)
tryCatch({
  Sys.setenv(R_MAX_VSIZE = "50Gb")
  cat("Attempted to set R_MAX_VSIZE to 50Gb\n")
}, error = function(e) {
  cat("Note: Could not set R_MAX_VSIZE (this is normal on some systems)\n")
})

# ============================================
# 1. DOWNLOAD DATA FROM GEO
# ============================================

cat("Downloading data from GEO...\n")

# Download metadata with memory-efficient options
# AnnotGPL = FALSE: Don't download annotation platform data (saves memory)
cat("Note: This may take a few minutes and use significant memory...\n")
cat("Using memory-efficient options (AnnotGPL = FALSE)...\n")

# Initialize sample_info
sample_info <- NULL

# First attempt: Use AnnotGPL = FALSE to save memory
tryCatch({
  gse <- getGEO("GSE60450", 
                GSEMatrix = TRUE,
                AnnotGPL = FALSE)  # Don't download annotation platform (saves memory)
  gse <- gse[[1]]
  sample_info <- pData(gse)
  
  # Clear large objects from memory immediately after use
  rm(gse)
  gc(verbose = FALSE)
  
  cat("Metadata downloaded successfully!\n")
}, error = function(e) {
  cat("\nFirst download attempt failed with error:\n")
  cat(e$message, "\n\n")
  cat("Trying alternative method: downloading series matrix file directly...\n")
  
  # Alternative: Download series matrix file directly (more memory efficient)
  tryCatch({
    series_file <- "GSE60450_series_matrix.txt.gz"
    if (!file.exists(series_file)) {
      cat("Downloading series matrix file...\n")
      download.file(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60450/matrix/GSE60450_series_matrix.txt.gz",
        series_file,
        mode = "wb"
      )
    }
    
    # Parse the series matrix file
    cat("Parsing series matrix file...\n")
    gse <- getGEO(filename = series_file, AnnotGPL = FALSE)
    gse <- gse[[1]]
    sample_info <<- pData(gse)  # Use <<- to assign to parent scope
    
    # Clear memory
    rm(gse)
    gc(verbose = FALSE)
    
    cat("Metadata downloaded successfully using alternative method!\n")
  }, error = function(e2) {
    cat("\nAlternative download method also failed.\n")
    cat("Error:", e2$message, "\n\n")
    cat("Trying ultra memory-efficient method: parsing metadata only...\n")
    
    # Ultra memory-efficient: Parse only metadata lines from series matrix file
    tryCatch({
      series_file <- "GSE60450_series_matrix.txt.gz"
      if (!file.exists(series_file)) {
        cat("Series matrix file not found. Downloading...\n")
        download.file(
          "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60450/matrix/GSE60450_series_matrix.txt.gz",
          series_file,
          mode = "wb"
        )
      }
      
      cat("Parsing metadata from series matrix file (memory-efficient method)...\n")
      
      # Read file line by line, extract only metadata (lines starting with "!")
      con <- gzfile(series_file, "r")
      metadata_lines <- character()
      sample_names <- character()
      
      # Read until we hit the data table header (line starting with "ID_REF")
      while (TRUE) {
        line <- readLines(con, n = 1)
        if (length(line) == 0) break
        
        # Stop when we reach the data table
        # Look for ID_REF or any line that looks like a data header (starts with gene ID, has many tabs)
        if (grepl("^ID_REF", line) || (!grepl("^!", line) && length(strsplit(line, "\t")[[1]]) > 5)) {
          # Extract sample names from header
          header_parts <- strsplit(line, "\t")[[1]]
          # Remove the first column (ID_REF or gene ID)
          sample_names <- header_parts[-1]
          # Clean up sample names (remove quotes)
          sample_names <- gsub('^"|"$', '', sample_names)
          sample_names <- trimws(sample_names)
          cat("Extracted", length(sample_names), "sample names from data header\n")
          break
        }
        
        # Collect metadata lines
        if (grepl("^!", line)) {
          metadata_lines <- c(metadata_lines, line)
        }
      }
      close(con)
      
      # Parse metadata into a data frame
      cat("Extracting sample information from metadata...\n")
      
      # Parse each metadata line - values are tab-separated and may be quoted
      metadata_list <- list()
      n_samples <- 0
      characteristics_counter <- 0  # Track multiple characteristics lines
      
      # First pass: determine number of samples and extract key metadata
      for (line in metadata_lines) {
        if (grepl("^!Sample_", line)) {
          # Split by tab
          parts <- strsplit(line, "\t")[[1]]
          
          # Extract key from first part
          # Format can be: "!Sample_key = "value1"" OR "!Sample_key"	"value1"	"value2"
          first_part <- parts[1]
          
          # Check if key has " = " separator
          if (grepl(" = ", first_part)) {
            # Format: "!Sample_key = "value1""
            key_parts <- strsplit(first_part, " = ")[[1]]
            key <- sub("^!Sample_", "", key_parts[1])
            # Values start from second part (after " = ")
            values <- parts[-1]
            # First value might be in key_parts[2] if it's before the first tab
            if (length(key_parts) > 1) {
              first_value <- key_parts[2]
              # Remove quotes if present
              first_value <- gsub('^"|"$', '', first_value)
              values <- c(first_value, values)
            }
          } else {
            # Format: "!Sample_key"	"value1"	"value2"
            key <- sub("^!Sample_", "", first_part)
            # All values are in parts[-1]
            values <- parts[-1]
          }
          
          # Remove quotes and trim whitespace from all values
          values <- gsub('^"|"$', '', values)
          values <- trimws(values)
          
          # Store number of samples from first valid line
          if (n_samples == 0 && length(values) > 0) {
            n_samples <- length(values)
            cat(paste("Found", n_samples, "samples\n"))
          }
          
          # Handle multiple characteristics lines specially
          if (grepl("characteristics", key, ignore.case = TRUE)) {
            characteristics_counter <- characteristics_counter + 1
            # Check what type of characteristic this is
            if (any(grepl("immunophenotype", values, ignore.case = TRUE))) {
              # Extract cell type from values
              cell_type_values <- sub(".*immunophenotype:\\s*", "", values, ignore.case = TRUE)
              cell_type_values <- sub("\\s+cell.*", "", cell_type_values, ignore.case = TRUE)
              metadata_list[["immunophenotype"]] <- tolower(cell_type_values)
            } else if (any(grepl("developmental.stage|stage", values, ignore.case = TRUE))) {
              # Extract stage from values
              stage_values <- sub(".*developmental.stage:\\s*", "", values, ignore.case = TRUE)
              stage_values <- sub("\\s+day.*", "", stage_values, ignore.case = TRUE)
              # Normalize stage names
              stage_values <- case_when(
                grepl("virgin", stage_values, ignore.case = TRUE) ~ "virgin",
                grepl("pregnant|18", stage_values, ignore.case = TRUE) ~ "pregnant",
                grepl("lactat", stage_values, ignore.case = TRUE) ~ "lactate",
                TRUE ~ stage_values
              )
              metadata_list[["developmental_stage"]] <- stage_values
            }
          } else {
            # Clean up key name (keep underscores, remove other special chars)
            key_clean <- gsub("[^A-Za-z0-9_]", "", key)
            
            # Store if we have the right number of values
            if (length(values) == n_samples && n_samples > 0) {
              metadata_list[[key_clean]] <- values
            }
          }
        }
      }
      
      # If we couldn't determine from metadata, use sample_names length
      if (n_samples == 0 && length(sample_names) > 0) {
        n_samples <- length(sample_names)
        cat(paste("Using", n_samples, "samples from data header\n"))
      }
      
      # Create sample_info data frame
      if (length(metadata_list) > 0 && n_samples > 0) {
        # Convert to data frame
        sample_info <- as.data.frame(metadata_list, stringsAsFactors = FALSE)
        
        # Ensure we have the right number of rows
        if (nrow(sample_info) != n_samples) {
          cat("Warning: Row count mismatch. Adjusting...\n")
          if (nrow(sample_info) > n_samples) {
            sample_info <- sample_info[1:n_samples, ]
          }
        }
        
        # Set row names to geo_accession if available
        if ("geoaccession" %in% colnames(sample_info)) {
          rownames(sample_info) <- sample_info$geoaccession
        } else if (length(sample_names) == nrow(sample_info)) {
          rownames(sample_info) <- sample_names
        } else {
          rownames(sample_info) <- paste0("Sample", 1:nrow(sample_info))
        }
        
        cat("Successfully parsed metadata:\n")
        cat("  Columns:", paste(colnames(sample_info), collapse = ", "), "\n")
        cat("  Rows:", nrow(sample_info), "\n")
        
      } else {
        # Fallback: create minimal metadata from sample names
        cat("Warning: Could not parse metadata lines. Creating minimal metadata from sample names.\n")
        if (length(sample_names) > 0) {
          sample_info <- data.frame(
            geo_accession = sample_names,
            stringsAsFactors = FALSE
          )
          rownames(sample_info) <- sample_names
        } else {
          stop("Cannot create metadata: no sample names found in series matrix file.")
        }
      }
      
      sample_info <<- sample_info
      cat("Metadata extracted successfully using memory-efficient method!\n")
      cat(paste("Extracted", ncol(sample_info), "metadata columns for", nrow(sample_info), "samples\n"))
      
    }, error = function(e3) {
      cat("\nUltra memory-efficient method failed.\n")
      cat("Error:", e3$message, "\n")
      cat("Error details:\n")
      print(e3)
      cat("\n")
      
      # Check if file exists and try to diagnose
      if (file.exists(series_file)) {
        cat("Series matrix file exists. Checking format...\n")
        tryCatch({
          con <- gzfile(series_file, "r")
          first_lines <- readLines(con, n = 50)
          close(con)
          sample_lines <- grep("^!Sample_", first_lines, value = TRUE)
          cat("Found", length(sample_lines), "sample metadata lines\n")
          if (length(sample_lines) > 0) {
            cat("First sample line:", substr(sample_lines[1], 1, 100), "...\n")
          }
        }, error = function(e) {
          cat("Could not read file for diagnosis:", e$message, "\n")
        })
      } else {
        cat("Series matrix file does not exist.\n")
      }
      
      cat("\nSOLUTIONS:\n")
      cat("1. Increase R memory limit (Mac/Linux):\n")
      cat("   Before running R, increase memory limit in terminal:\n")
      cat("   ulimit -v 50000000  # Sets virtual memory limit to ~50GB\n")
      cat("   Then restart R/RStudio\n")
      cat("   OR in R, try: Sys.setenv(R_MAX_VSIZE = \"50Gb\")\n")
      cat("2. Download files manually:\n")
      cat("   - Visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450\n")
      cat("   - Download 'Series Matrix File(s)'\n")
      cat("   - Place GSE60450_series_matrix.txt.gz in working directory\n")
      cat("   - Re-run this script\n")
      cat("3. Use a machine with more RAM (32GB+ recommended)\n")
      cat("4. Contact workshop organizers for pre-processed metadata file\n\n")
      stop("GEO download failed. See solutions above.")
    })
  })
})

# Check if download was successful
if (is.null(sample_info)) {
  stop("Failed to download GEO metadata. See error messages above for solutions.")
}

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

# Check if sample_info exists and has data
if (is.null(sample_info) || nrow(sample_info) == 0) {
  stop("ERROR: sample_info is NULL or empty. Metadata download/extraction failed.\n",
       "Please check the error messages above and ensure the series matrix file was downloaded correctly.")
}

cat("Sample info columns:", paste(colnames(sample_info), collapse = ", "), "\n")
cat("Number of samples:", nrow(sample_info), "\n")

# Try to find the correct column names (they may vary)
geo_col <- if("geo_accession" %in% colnames(sample_info)) {
  sample_info$geo_accession
} else if("geoaccession" %in% colnames(sample_info)) {
  sample_info$geoaccession
} else {
  rownames(sample_info)
}

# Try to find description/source name column (may be named differently)
desc_col <- NULL
desc_col_name <- NULL
for (col in colnames(sample_info)) {
  if (grepl("description|sourcename|title", col, ignore.case = TRUE)) {
    desc_col <- sample_info[[col]]
    desc_col_name <- col
    break
  }
}

# Try to find immunophenotype/cell type column
celltype_col <- NULL
if ("immunophenotype" %in% colnames(sample_info)) {
  celltype_col <- sample_info$immunophenotype
} else {
  for (col in colnames(sample_info)) {
    if (grepl("immunophenotype", col, ignore.case = TRUE)) {
      celltype_col <- sample_info[[col]]
      break
    }
  }
}
# If not found, try to extract from source name
if (is.null(celltype_col) && !is.null(desc_col)) {
  # Source name format: "mammary gland, luminal cells, virgin" or "Luminal virgin #1"
  if (any(grepl("luminal|basal", desc_col, ignore.case = TRUE))) {
    # We'll extract from desc_col later
    celltype_col <- desc_col
  }
}

# Try to find developmental stage column
status_col <- NULL
if ("developmental_stage" %in% colnames(sample_info)) {
  status_col <- sample_info$developmental_stage
} else {
  for (col in colnames(sample_info)) {
    if (grepl("developmentalstage|stage", col, ignore.case = TRUE)) {
      status_col <- sample_info[[col]]
      break
    }
  }
}
# If not found, try to extract from source name or title
if (is.null(status_col) && !is.null(desc_col)) {
  # Source name or title contains stage info
  if (any(grepl("virgin|pregnant|lactat", desc_col, ignore.case = TRUE))) {
    status_col <- desc_col
  }
}

# Extract sample names from description/title
# The counts file uses descriptive names like "Luminal virgin #1", "Basal 18.5 dP #1", etc.
# So we should use the title directly to match the counts file format
if (!is.null(desc_col)) {
  # Check if desc_col contains titles like "Luminal virgin #1"
  # Format from series matrix: "Luminal virgin #1", "Basal 18.5 dP #1", "Luminal 2 dL #1"
  if (any(grepl("#", desc_col))) {
    # Looks like title format - use directly
    sample_names <- desc_col
  } else {
    # Might be source name format - try to construct title
    # Format: "mammary gland, luminal cells, virgin"
    # Construct: "Luminal virgin #1", "Basal 18.5 dP #1", etc.
    # But we need replicate numbers - use position for now
    sample_names <- paste(
      ifelse(grepl("luminal", desc_col, ignore.case = TRUE), "Luminal", "Basal"),
      case_when(
        grepl("virgin", desc_col, ignore.case = TRUE) ~ "virgin",
        grepl("pregnant|18", desc_col, ignore.case = TRUE) ~ "18.5 dP",
        grepl("lactat|2 d", desc_col, ignore.case = TRUE) ~ "2 dL",
        TRUE ~ "unknown"
      ),
      paste0("#", rep(1:2, each = 6)),  # 2 replicates per group, 6 groups
      sep = " "
    )
  }
} else {
  # Fallback: construct from cell type and status
  sample_names <- paste(
    ifelse(cell_types == "luminal", "Luminal", "Basal"),
    case_when(
      statuses == "virgin" ~ "virgin",
      statuses == "pregnant" ~ "18.5 dP",
      statuses == "lactate" ~ "2 dL",
      TRUE ~ "unknown"
    ),
    paste0("#", rep(1:2, each = 6)),
    sep = " "
  )
}

# Determine cell type
if (!is.null(celltype_col)) {
  # Check if celltype_col contains the actual cell type or needs parsing
  if (any(grepl("luminal cell|basal cell", celltype_col, ignore.case = TRUE))) {
    # Contains full description like "luminal cell population"
    cell_types <- ifelse(grepl("luminal", celltype_col, ignore.case = TRUE), 
                         "luminal", "basal")
  } else {
    # Might be in source name format, extract from there
    cell_types <- ifelse(grepl("luminal", celltype_col, ignore.case = TRUE), 
                         "luminal", "basal")
  }
} else if (!is.null(desc_col)) {
  # Extract from description/source name
  # Format: "mammary gland, luminal cells, virgin" or "Luminal virgin #1"
  cell_types <- ifelse(grepl("luminal", desc_col, ignore.case = TRUE), 
                       "luminal", "basal")
} else {
  # Default: assume basal (reference level)
  cell_types <- rep("basal", length(geo_col))
  cat("Warning: Could not determine cell types, defaulting to 'basal'\n")
}

# Determine status
if (!is.null(status_col)) {
  # Parse status from the status column
  # Format might be: "virgin", "18.5 day pregnancy", "2 day lactation", "developmental stage: 2 day lactation"
  statuses <- case_when(
    grepl("virgin", status_col, ignore.case = TRUE) ~ "virgin",
    grepl("pregnant|18\\.5|18.5", status_col, ignore.case = TRUE) ~ "pregnant",
    grepl("lactat|2 dL|2 day lact|2dL|2 d L", status_col, ignore.case = TRUE) ~ "lactate",
    TRUE ~ NA_character_
  )
} else if (!is.null(desc_col)) {
  # Extract from description/title
  # Format: "Luminal virgin #1", "Basal 18.5 dP #1", "Luminal 2 dL #1"
  # Note: "2 dL" means "2 day lactation"
  statuses <- case_when(
    grepl("virgin", desc_col, ignore.case = TRUE) ~ "virgin",
    grepl("pregnant|18\\.5|18.5|dP", desc_col, ignore.case = TRUE) ~ "pregnant",
    grepl("lactat|2 dL|2 day lact|2dL|dL", desc_col, ignore.case = TRUE) ~ "lactate",
    TRUE ~ NA_character_
  )
} else {
  stop("Cannot determine developmental stage. Please check metadata columns.\n",
       "Available columns: ", paste(colnames(sample_info), collapse = ", "))
}

# Check for NA statuses and try to fix them
if (any(is.na(statuses))) {
  cat("Warning: Some samples have NA status. Attempting to fix...\n")
  na_indices <- which(is.na(statuses))
  if (!is.null(desc_col)) {
    # Try again with more permissive patterns
    for (i in na_indices) {
      desc <- desc_col[i]
      if (grepl("virgin", desc, ignore.case = TRUE)) {
        statuses[i] <- "virgin"
      } else if (grepl("pregnant|18|dP", desc, ignore.case = TRUE)) {
        statuses[i] <- "pregnant"
      } else if (grepl("lactat|dL|2 d", desc, ignore.case = TRUE)) {
        statuses[i] <- "lactate"
      }
    }
  }
  # Report remaining NAs
  if (any(is.na(statuses))) {
    cat("Still have NA status for samples:", paste(which(is.na(statuses)), collapse = ", "), "\n")
    if (!is.null(desc_col)) {
      cat("Their descriptions:", paste(desc_col[is.na(statuses)], collapse = "; "), "\n")
    }
  }
}

# Create metadata data frame
# Ensure SampleName doesn't have problematic characters for file I/O
# Replace spaces with underscores or use a simpler format
sample_names_clean <- sample_names
# If sample names contain spaces or special chars, create MCL1 codes instead
if (any(grepl("[^A-Za-z0-9-]", sample_names))) {
  cat("Warning: Sample names contain special characters. Creating MCL1 codes...\n")
  # Create codes based on cell type and status: L/B + V/P/L + number
  # L = Luminal, B = Basal
  # V = Virgin, P = Pregnant, L = Lactate
  # But we need to match the actual codes in counts file (MCL1-DG, etc.)
  # For now, use a simpler approach: match by position or use geo accession
  # Actually, let's extract from title if it has MCL1 pattern, otherwise use geo
  if (!is.null(desc_col) && any(grepl("MCL1", desc_col))) {
    sample_names_clean <- str_extract(desc_col, "MCL1[-.][A-Z]+")
    sample_names_clean <- gsub("\\.", "-", sample_names_clean)
    # Fill NAs with geo accession
    sample_names_clean[is.na(sample_names_clean)] <- geo_col[is.na(sample_names_clean)]
  } else {
    # Use geo accession as sample name (will match by position or prefix)
    sample_names_clean <- geo_col
  }
}

metadata <- data.frame(
  Sample = geo_col,
  SampleName = sample_names_clean,
  CellType = cell_types,
  Status = statuses,
  stringsAsFactors = FALSE
)

# Check for NA statuses and warn
if (any(is.na(metadata$Status))) {
  cat("Warning: Some samples have NA status. Checking...\n")
  na_status <- which(is.na(metadata$Status))
  cat("Samples with NA status:", paste(metadata$SampleName[na_status], collapse = ", "), "\n")
  if (!is.null(desc_col)) {
    cat("Their descriptions:", paste(desc_col[na_status], collapse = "; "), "\n")
  }
}

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
# The counts file has columns like "MCL1-DG_BC2CTUACXX_ACTTGA_L002_R1"
# We need to convert these to match the metadata SampleName format
# Metadata uses titles like "Luminal virgin #1" from the title column
sample_columns <- colnames(counts_only)
cat("Original sample column names:\n")
print(sample_columns[1:min(3, length(sample_columns))])

# First, try to match by position if metadata is already created
# Otherwise, we'll need to construct names from the counts file codes
# The counts file codes need to be mapped to descriptive names
# For now, we'll extract MCL1 codes and match them later to metadata
# Extract MCL1-XX or MCL1.XX pattern (handles both formats)
new_sample_names <- str_extract(sample_columns, "MCL1[-.][A-Z]+")
# Replace dot with dash if present to standardize
new_sample_names <- str_replace(new_sample_names, "\\.", "-")

# Check if extraction worked
if (any(is.na(new_sample_names))) {
  cat("Warning: Pattern extraction failed for some samples. Trying alternative method...\n")
  # Alternative: use first part before underscore, then extract pattern
  temp_names <- sub("_.*", "", sample_columns)
  # Try to extract MCL1-XX or MCL1.XX pattern
  new_sample_names <- str_extract(temp_names, "MCL1[-.][A-Z]+")
  new_sample_names <- str_replace(new_sample_names, "\\.", "-")
  
  # If still failing, just use the prefix before underscore
  na_indices <- is.na(new_sample_names)
  if (any(na_indices)) {
    new_sample_names[na_indices] <- temp_names[na_indices]
  }
}

# Final check
if (any(is.na(new_sample_names))) {
  stop("ERROR: Could not extract sample names from counts file columns.\n",
       "Original columns: ", paste(sample_columns, collapse = ", "), "\n",
       "Please check the counts file format.")
}

# Store the MCL1 codes temporarily - we'll match them to metadata titles later
colnames(counts_only) <- new_sample_names

cat("\nSample column names after cleaning (MCL1 codes):\n")
print(colnames(counts_only))
cat("Note: These will be matched to metadata SampleName after metadata is created\n")

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

# Try to match sample names - they might have different formats
# Counts: "MCL1-DG_BC2CTUACXX_ACTTGA_L002_R1"
# Metadata: "MCL1-DG" or similar
# Extract the prefix (MCL1-XX) from counts column names
counts_prefixes <- sub("_.*", "", colnames(counts_final))
cat("\nCounts prefixes (first part before underscore):\n")
print(counts_prefixes)

# Try matching by prefix
matched_indices <- match(metadata$SampleName, counts_prefixes)
if (any(is.na(matched_indices))) {
  # Try the other way - extract prefix from metadata sample names
  metadata_prefixes <- sub("_.*", "", metadata$SampleName)
  matched_indices <- match(metadata_prefixes, counts_prefixes)
  
  if (any(is.na(matched_indices))) {
    # Try exact match first
    matched_indices <- match(metadata$SampleName, colnames(counts_final))
    
    if (any(is.na(matched_indices))) {
      cat("\nWARNING: Sample name matching issues detected.\n")
      cat("Trying to match by position (assuming same order)...\n")
      
      if (length(metadata$SampleName) == length(colnames(counts_final))) {
        # Match by position
        matched_indices <- 1:length(metadata$SampleName)
        cat("Matched", length(matched_indices), "samples by position\n")
      } else {
        stop("Cannot match samples: different numbers of samples in metadata (", 
             length(metadata$SampleName), ") and counts (", 
             length(colnames(counts_final)), ")\n",
             "Metadata samples: ", paste(metadata$SampleName, collapse = ", "), "\n",
             "Counts columns: ", paste(colnames(counts_final), collapse = ", "), "\n")
      }
    }
  }
}

# Reorder columns to match metadata row order
# Since metadata SampleName uses descriptive titles and counts use MCL1 codes,
# we'll match by position (order is consistent) and update column names
if (all(!is.na(matched_indices))) {
  counts_ordered <- counts_final[, matched_indices, drop = FALSE]
  # Update column names to match metadata SampleName (descriptive titles)
  colnames(counts_ordered) <- metadata$SampleName
  cat("Successfully matched and reordered", ncol(counts_ordered), "samples\n")
  cat("Updated counts column names to match metadata SampleName\n")
} else if (length(metadata$SampleName) == length(colnames(counts_final))) {
  # Match by position if same number of samples
  cat("Matching samples by position (same order assumed)...\n")
  counts_ordered <- counts_final
  # Update column names to match metadata SampleName
  colnames(counts_ordered) <- metadata$SampleName
  cat("Successfully matched", ncol(counts_ordered), "samples by position\n")
} else {
  stop("Failed to match samples. Please check sample names.\n",
       "Metadata samples (", length(metadata$SampleName), "): ", 
       paste(metadata$SampleName, collapse = ", "), "\n",
       "Counts columns (", length(colnames(counts_final)), "): ", 
       paste(colnames(counts_final), collapse = ", "), "\n")
}

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

# Verify counts_ordered has sample columns
if (ncol(counts_ordered) == 0) {
  stop("ERROR: counts_ordered has no sample columns. This usually means sample names don't match between metadata and counts.\n",
       "Metadata sample names: ", paste(metadata$SampleName, collapse = ", "), "\n",
       "Counts column names: ", paste(colnames(counts_final), collapse = ", "), "\n",
       "Please check that sample names match between files.")
}

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

cat("Counts file will have", ncol(counts_with_id), "columns (EnsemblID, GeneSymbol, and", 
    ncol(counts_ordered), "sample columns)\n")

# ============================================
# 7. SAVE CLEAN FILES
# ============================================

cat("\nSaving cleaned files...\n")

# Write metadata file using write.table with proper formatting
# Ensure no NA values in Status column before writing
if (any(is.na(metadata$Status))) {
  cat("Warning: Some Status values are NA. These will be written as 'NA' string.\n")
  cat("Samples with NA Status:", paste(metadata$SampleName[is.na(metadata$Status)], collapse = ", "), "\n")
}

write.table(metadata, 
            file = "GSE60450_metadata_clean.txt",
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE,
            na = "")  # Write empty string instead of NA
cat("✓ Saved: GSE60450_metadata_clean.txt\n")
cat("  Samples:", nrow(metadata), "\n")
cat("  Columns:", paste(colnames(metadata), collapse = ", "), "\n")

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
