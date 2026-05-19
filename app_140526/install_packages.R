# Installs missing packages required by the RNA-seq dashboard.
# Run once, or let run_RNAseq_dashboard.bat call it automatically.

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_packages <- c(
  "shiny", "shinyFiles", "shinythemes", "DT", "ggplot2", "dplyr", "readr", "stringr", "tibble",
  "tidyr", "ggrepel", "RColorBrewer", "pheatmap", "ashr", "msigdbr", "VennDiagram"
)

bioc_packages <- c(
  "DESeq2", "tximport", "topGO", "org.At.tair.db", "GO.db", "AnnotationDbi",
  "SummarizedExperiment", "rrvgo"
)

install_missing_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing CRAN packages: ", paste(missing, collapse = ", "))
    install.packages(missing, dependencies = TRUE)
  } else {
    message("All CRAN packages are already installed.")
  }
}

install_missing_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing Bioconductor packages: ", paste(missing, collapse = ", "))
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  } else {
    message("All Bioconductor packages are already installed.")
  }
}

install_missing_cran(cran_packages)
install_missing_bioc(bioc_packages)

message("Package check finished.")
