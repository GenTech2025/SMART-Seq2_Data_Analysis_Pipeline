# =============================================================
# 📚 Library Installation Script
# This script installs all packages listed in the provided library calls.
# NOTE: Many of these packages (scater, SingleCellExperiment, etc.) 
# belong to Bioconductor and require specific installation steps.
# =============================================================

# Function to check if a package is installed and install it if not found
install_if_needed <- function(package_name, source) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    cat("\n[Attempting to install:", package_name, "]\n")
    
    # Check which repository to use (Bioconductor or CRAN)
    if (source == "BiocManager") {
      tryCatch({
        BiocManager::install(package_name)
      }, error = function(e) {
        cat("--- ERROR: Could not install", package_name, ". Check internet connection or dependencies.\n")
      })
    } else if (source == "CRAN") {
      tryCatch({
        install.packages(package_name)
      }, error = function(e) {
        cat("--- ERROR: Could not install", package_name, ". Check internet connection or dependencies.\n")
      })
    } else {
      cat("--- WARNING: Unknown source for package:", package_name, "\n")
    }
  } else {
    cat("[INFO] Package", package_name, "is already installed. Skipping.\n")
  }
}


# -------------------------------------------
# STEP 1: Install the Bioconductor Manager (Essential prerequisite)
# -------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cat("\n=============================================")
cat("\n[STARTING PACKAGE INSTALLATION PROCESS]")
cat("\n=============================================\n")


# -------------------------------------------
# STEP 2: Install Bioconductor Packages
# (These are single-cell/genomics packages)
# -------------------------------------------

bioc_packages <- c(
  "SingleCellExperiment", "scater", "scran", "biomaRt", 
  "AnnotationDbi", "Seurat", "fgsea", "GSEABase", 
  "msigdb" # msigdb is often a data dependency, but we include it here
)

cat("\n\n[--- Installing BIOCONDUCTOR Packages ---\n]")
for (pkg in bioc_packages) {
  install_if_needed(pkg, "BiocManager")
}


# -------------------------------------------
# STEP 3: Install CRAN Packages
# (These are general purpose packages like ggplot2 and tidyverse)
# -------------------------------------------

cran_packages <- c(
  "tidyverse", "ggplot2", "jsonlite", "pheatmap", "fgsea", # fgsea is sometimes bioconductor, but often installed via cran dependencies too.
  "scMerge", # Note: scMerge might require specific setup, using standard install for attempt
  "org.Hs.eg.db" # This package is technically a Bioc version of an ID mapping tool, but included here for completeness if the user's environment treats it as CRAN-friendly.
)

cat("\n\n[--- Installing CRAN Packages ---\n]")
for (pkg in cran_packages) {
  install_if_needed(pkg, "CRAN")
}


# -------------------------------------------
# STEP 4: Final Cleanup and Confirmation
# -------------------------------------------

cat("\n=============================================")
cat("\n✅ Installation script finished!")
cat("\nPlease restart your R session to ensure all libraries are properly loaded.")
cat("=============================================\n")
# 