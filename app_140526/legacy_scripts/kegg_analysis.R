# KEGG Analysis Helper Functions for RNA-seq Dashboard
# Adapted from user's run_KEGG_010425.r and pathview scripts.

suppressPackageStartupMessages({
  library(KEGGREST)
  library(dplyr)
  library(ggplot2)
})

########################################################
# 1. Fetch & Cache KEGG Genes for Arabidopsis
########################################################
get_kegg_genes_cached <- function(cache_file = "description_files/kegg_ath_genes_cache.rds") {
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  
  message("Downloading KEGG pathways for Arabidopsis (ath). This may take a few minutes...")
  pathways.list <- KEGGREST::keggList("pathway", "ath")
  pathway.codes <- sub("path:", "", names(pathways.list))
  
  genes.by.pathway <- sapply(pathway.codes, function(pwid) {
    pw <- KEGGREST::keggGet(pwid)
    if (is.null(pw[[1]]$GENE)) return(NA)
    pw2 <- pw[[1]]$GENE[c(TRUE, FALSE)] # Odd indices contain gene IDs
    pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = TRUE), function(x) x[1]))
    return(pw2)
  }, simplify = FALSE)
  
  # Clean up pathway names
  pathway.names <- as.character(pathways.list)
  pathway.names <- gsub(" - Arabidopsis thaliana \\(thale cress\\)", "", pathway.names)
  names(pathway.names) <- pathway.codes
  
  res <- list(
    genes_by_pathway = genes.by.pathway,
    pathway_names = pathway.names
  )
  
  # Try to save to cache
  tryCatch(saveRDS(res, cache_file), error = function(e) warning("Could not save KEGG cache: ", e$message))
  return(res)
}


########################################################
# 2. Run Wilcoxon Enrichment Test
########################################################
run_kegg_enrichment <- function(de_df, pvalue_cutoff = 0.05) {
  # Requires gene_id, padj/pValue, log2FoldChange
  if (!"pValue" %in% names(de_df)) de_df$pValue <- de_df$padj
  
  # Remove NAs
  de_df <- de_df[!is.na(de_df$gene_id) & !is.na(de_df$pValue), ]
  
  # Ensure pValues > 0 for log math if needed, and cap at 0.999 for non-sig
  geneList <- de_df$pValue
  names(geneList) <- de_df$gene_id
  
  # Also store log2FC for up/down classification
  lfcList <- de_df$log2FoldChange
  names(lfcList) <- de_df$gene_id
  
  kegg_data <- get_kegg_genes_cached()
  genes.by.pathway <- kegg_data$genes_by_pathway
  pathways.list <- kegg_data$pathway_names
  
  # Wilcoxon test for each pathway
  res_list <- lapply(names(genes.by.pathway), function(pathway) {
    pathway.genes <- genes.by.pathway[[pathway]]
    if (any(is.na(pathway.genes))) return(NULL)
    
    list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
    if (length(list.genes.in.pathway) < 3) return(NULL) # Skip very small pathways
    
    list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
    
    scores.in.pathway <- geneList[list.genes.in.pathway]
    scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
    
    # Wilcoxon test: are p-values in pathway significantly smaller?
    p.value <- suppressWarnings(
      wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
    )
    
    sig_genes <- list.genes.in.pathway[scores.in.pathway < pvalue_cutoff]
    up_genes <- sum(lfcList[sig_genes] > 0, na.rm = TRUE)
    down_genes <- sum(lfcList[sig_genes] < 0, na.rm = TRUE)
    
    data.frame(
      pathway.code = pathway,
      pathway.name = pathways.list[pathway],
      p.value = p.value,
      Significant = length(sig_genes),
      Upregulated = up_genes,
      Downregulated = down_genes,
      Annotated = length(list.genes.in.pathway),
      stringsAsFactors = FALSE
    )
  })
  
  outdat <- do.call(rbind, res_list)
  if (is.null(outdat) || nrow(outdat) == 0) return(NULL)
  
  outdat <- outdat[order(outdat$p.value), ]
  rownames(outdat) <- NULL
  return(outdat)
}

########################################################
# 3. Bubble Plot for Enriched Pathways
########################################################
plot_kegg_bubble <- function(kegg_res_df, p_value_threshold = 0.05) {
  if (is.null(kegg_res_df) || nrow(kegg_res_df) == 0) return(NULL)
  
  plot_df <- kegg_res_df[kegg_res_df$p.value <= p_value_threshold, ]
  if (nrow(plot_df) == 0) return(NULL)
  
  up_df <- plot_df[plot_df$Upregulated > 0, ]
  if (nrow(up_df) > 20) up_df <- head(up_df, 20)
  
  down_df <- plot_df[plot_df$Downregulated > 0, ]
  if (nrow(down_df) > 20) down_df <- head(down_df, 20)
  
  if (nrow(up_df) == 0 && nrow(down_df) == 0) return(NULL)
  
  up_plot <- NULL
  if (nrow(up_df) > 0) {
    up_plot <- ggplot(up_df, aes(x = Upregulated, y = reorder(pathway.name, Upregulated), color = p.value, size = Annotated)) +
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
    down_plot <- ggplot(down_df, aes(x = Downregulated, y = reorder(pathway.name, Downregulated), color = p.value, size = Annotated)) +
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
  
  # Combine plots using gridExtra or just return one if the other is empty
  if (!is.null(up_plot) && !is.null(down_plot)) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(up_plot + down_plot)
    } else if (requireNamespace("gridExtra", quietly = TRUE)) {
      return(gridExtra::arrangeGrob(up_plot, down_plot, ncol = 2))
    } else {
      # Fallback: just return up plot with a warning
      warning("Please install 'patchwork' or 'gridExtra' to see side-by-side plots.")
      return(up_plot)
    }
  } else if (!is.null(up_plot)) {
    return(up_plot)
  } else {
    return(down_plot)
  }
}
