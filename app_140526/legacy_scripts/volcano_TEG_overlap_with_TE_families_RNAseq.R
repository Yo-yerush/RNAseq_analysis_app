suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# If you run this script from the app root, load shared helpers (standardize_de_table, clean_tair_id).
if (file.exists(file.path("R", "helpers.R"))) {
  source(file.path("R", "helpers.R"), local = TRUE)
}

normalize_gene_id <- function(x) {
  if (exists("clean_tair_id", mode = "function")) return(clean_tair_id(x))
  x <- trimws(as.character(x))
  x <- sub("\\..*$", "", x)
  toupper(x)
}

load_te_gene_map <- function(description_file) {
  if (!file.exists(description_file)) {
    stop("Description file not found: ", description_file)
  }
  desc <- read.csv(description_file, check.names = FALSE)
  needed <- c("gene_id", "gene_model_type", "Derives_from")
  missing <- setdiff(needed, names(desc))
  if (length(missing)) {
    stop("Description file is missing columns: ", paste(missing, collapse = ", "))
  }
  desc <- desc %>%
    dplyr::filter(gene_model_type == "transposable_element_gene") %>%
    dplyr::select(gene_id, Derives_from) %>%
    dplyr::filter(!is.na(Derives_from) & Derives_from != "")
  desc$gene_id <- normalize_gene_id(desc$gene_id)
  desc
}

load_te_file <- function(te_file) {
  if (!file.exists(te_file)) {
    stop("TAIR10 TE file not found: ", te_file)
  }
  te <- read.csv(te_file, sep = "\t", check.names = FALSE)
  needed <- c("Transposon_Name", "Transposon_Family", "Transposon_Super_Family")
  missing <- setdiff(needed, names(te))
  if (length(missing)) {
    stop("TAIR10 TE file is missing columns: ", paste(missing, collapse = ", "))
  }
  te %>%
    dplyr::rename(Derives_from = Transposon_Name) %>%
    dplyr::select(Derives_from, Transposon_Family, Transposon_Super_Family)
}

build_te_overlap_table <- function(de_df,
                                   description_file = file.path("description_files", "Methylome.At_description_file.csv.gz"),
                                   te_file = file.path("description_files", "TAIR10_Transposable_Elements.txt")) {
  if (is.null(de_df) || nrow(de_df) == 0) stop("DE table is empty.")
  if (exists("standardize_de_table", mode = "function")) {
    de_df <- standardize_de_table(de_df)
  }
  needed <- c("gene_id", "log2FoldChange", "padj")
  missing <- setdiff(needed, names(de_df))
  if (length(missing)) stop("DE table is missing columns: ", paste(missing, collapse = ", "))
  de_df$gene_id <- normalize_gene_id(de_df$gene_id)
  de_df <- de_df[!is.na(de_df$padj), c("gene_id", "log2FoldChange", "padj"), drop = FALSE]

  gene_map <- load_te_gene_map(description_file)
  te_meta <- load_te_file(te_file)
  te_map <- merge(gene_map, te_meta, by = "Derives_from")
  merge(te_map, de_df, by = "gene_id")
}

make_retro_te_volcano <- function(de_df,
                                  description_file = file.path("description_files", "Methylome.At_description_file.csv.gz"),
                                  te_file = file.path("description_files", "TAIR10_Transposable_Elements.txt"),
                                  families = NULL,
                                  super_families = c("LTR/Copia", "LTR/Gypsy", "LINE/L1"),
                                  padj_cutoff = 0.05,
                                  lfc_cutoff = 1,
                                  point_size = 0.75,
                                  point_alpha = 0.4) {
  te_df <- build_te_overlap_table(de_df, description_file, te_file)
  if (!is.null(families) && length(families) > 0) {
    retro_df <- te_df %>% dplyr::filter(Transposon_Family %in% families)
    color_col <- "Transposon_Family"
    label_name <- "Transposon\nFamily"
  } else if (!is.null(super_families) && length(super_families) > 0) {
    retro_df <- te_df %>% dplyr::filter(Transposon_Super_Family %in% super_families)
    color_col <- "Transposon_Super_Family"
    label_name <- "Transposon\nSuper-Family"
  } else {
    stop("Select at least one TE family or super-family.")
  }
  if (nrow(retro_df) == 0) stop("No TE genes found for the requested families.")

  retro_df$geneCat <- ifelse(
    retro_df$padj < padj_cutoff & retro_df$log2FoldChange > lfc_cutoff, "Upregulated",
    ifelse(retro_df$padj < padj_cutoff & retro_df$log2FoldChange < -lfc_cutoff, "Downregulated", "nonDE")
  )
  retro_df$geneCat <- factor(retro_df$geneCat, levels = c("Upregulated", "Downregulated", "nonDE"))

  p <- ggplot(retro_df, aes_string(x = "log2FoldChange", y = "-log10(padj)", color = color_col)) +
    geom_point(alpha = point_alpha, size = point_size) +
    xlab("log2(Fold-Change)") +
    ylab("-log10(padj)") +
    theme_bw(base_family = "serif") +
    theme(legend.position = "right") + 
    labs(color = label_name) +
    scale_y_continuous(breaks = c(0, 5, 10, 15), labels = c("0", "5", "10", "15")) +
    geom_segment(aes(x = lfc_cutoff, y = -log10(padj_cutoff), xend = lfc_cutoff, yend = Inf),
                 col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed") +
    geom_segment(aes(x = -lfc_cutoff, y = -log10(padj_cutoff), xend = -lfc_cutoff, yend = Inf),
                 col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed") +
    geom_segment(aes(x = lfc_cutoff, y = -log10(padj_cutoff), xend = Inf, yend = -log10(padj_cutoff)),
                 col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed") +
    geom_segment(aes(x = -Inf, y = -log10(padj_cutoff), xend = -lfc_cutoff, yend = -log10(padj_cutoff)),
                 col = "gray20", alpha = 0.6, size = 0.4, linetype = "dashed")

  list(data = retro_df, plot = p)
}

# Example usage (from app root):
# df <- read.csv("example_data/all_genes_results_mto1_vs_wt.csv", check.names = FALSE)
# out <- make_retro_te_volcano(df)
# out$plot

run_te_enrichment <- function(de_df,
                              description_file = file.path("description_files", "Methylome.At_description_file.csv.gz"),
                              te_file = file.path("description_files", "TAIR10_Transposable_Elements.txt"),
                              pvalue_cutoff = 0.05) {
  
  te_df <- build_te_overlap_table(de_df, description_file, te_file)
  if (nrow(te_df) == 0) return(NULL)
  
  if (!"pValue" %in% names(de_df)) de_df$pValue <- de_df$padj
  # Use the whole transcriptome for background p-values
  all_genes <- de_df[!is.na(de_df$pValue) & !is.na(de_df$gene_id), ]
  geneList <- all_genes$pValue
  names(geneList) <- all_genes$gene_id
  
  lfcList <- all_genes$log2FoldChange
  names(lfcList) <- all_genes$gene_id
  
  superfamilies <- unique(te_df$Transposon_Super_Family)
  
  res_list <- lapply(superfamilies, function(sf) {
    sf_genes <- te_df$gene_id[te_df$Transposon_Super_Family == sf]
    list.genes.in.sf <- intersect(names(geneList), sf_genes)
    if (length(list.genes.in.sf) < 3) return(NULL)
    
    list.genes.not.in.sf <- setdiff(names(geneList), list.genes.in.sf)
    
    scores.in.sf <- geneList[list.genes.in.sf]
    scores.not.in.sf <- geneList[list.genes.not.in.sf]
    
    p.value <- suppressWarnings(
      wilcox.test(scores.in.sf, scores.not.in.sf, alternative = "less")$p.value
    )
    
    sig_genes <- list.genes.in.sf[scores.in.sf < pvalue_cutoff]
    up_genes <- sum(lfcList[sig_genes] > 0, na.rm = TRUE)
    down_genes <- sum(lfcList[sig_genes] < 0, na.rm = TRUE)
    
    data.frame(
      superfamily = sf,
      p.value = p.value,
      Significant = length(sig_genes),
      Upregulated = up_genes,
      Downregulated = down_genes,
      Annotated = length(list.genes.in.sf),
      stringsAsFactors = FALSE
    )
  })
  
  outdat <- do.call(rbind, res_list)
  if (is.null(outdat) || nrow(outdat) == 0) return(NULL)
  
  outdat <- outdat[order(outdat$p.value), ]
  rownames(outdat) <- NULL
  return(outdat)
}

plot_te_enrichment <- function(te_res_df, p_value_threshold = 0.05) {
  if (is.null(te_res_df) || nrow(te_res_df) == 0) return(NULL)
  
  plot_df <- te_res_df[te_res_df$p.value <= p_value_threshold, ]
  if (nrow(plot_df) == 0) return(NULL)
  
  up_df <- plot_df[plot_df$Upregulated > 0, ]
  if (nrow(up_df) > 20) up_df <- head(up_df, 20)
  
  down_df <- plot_df[plot_df$Downregulated > 0, ]
  if (nrow(down_df) > 20) down_df <- head(down_df, 20)
  
  if (nrow(up_df) == 0 && nrow(down_df) == 0) return(NULL)
  
  up_plot <- NULL
  if (nrow(up_df) > 0) {
    up_plot <- ggplot(up_df, aes(x = Upregulated, y = reorder(superfamily, Upregulated), color = p.value, size = Annotated)) +
      geom_point() +
      scale_color_gradient("p.value", low = "#cf534c", high = "black") +
      theme_bw() +
      theme(
        legend.key.size = unit(0.25, 'cm'),
        legend.title = element_text(size = 9.5),
        text = element_text(family = "serif")
      ) +
      labs(title = "Upregulated", x = "Significant Genes", y = "") +
      guides(color = guide_colorbar(order = 1, barheight = 4))
  }
  
  down_plot <- NULL
  if (nrow(down_df) > 0) {
    down_plot <- ggplot(down_df, aes(x = Downregulated, y = reorder(superfamily, Downregulated), color = p.value, size = Annotated)) +
      geom_point() +
      scale_color_gradient("p.value", low = "#6397eb", high = "black") +
      theme_bw() +
      theme(
        legend.key.size = unit(0.25, 'cm'),
        legend.title = element_text(size = 9.5),
        text = element_text(family = "serif")
      ) +
      labs(title = "Downregulated", x = "Significant Genes", y = "") +
      guides(color = guide_colorbar(order = 1, barheight = 4))
  }
  
  if (!is.null(up_plot) && !is.null(down_plot)) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(up_plot + down_plot)
    } else if (requireNamespace("gridExtra", quietly = TRUE)) {
      return(gridExtra::arrangeGrob(up_plot, down_plot, ncol = 2))
    } else {
      warning("Please install 'patchwork' or 'gridExtra' to see side-by-side plots.")
      return(up_plot)
    }
  } else if (!is.null(up_plot)) {
    return(up_plot)
  } else {
    return(down_plot)
  }
}
