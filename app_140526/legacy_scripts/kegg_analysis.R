# KEGG Analysis Helper Functions for RNA-seq Dashboard
# Adapted from user's run_KEGG_010425.r and pathview scripts.

suppressPackageStartupMessages({
  library(KEGGREST)
  library(dplyr)
  library(ggplot2)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

if (!exists("plot_theme_choice", mode = "function")) {
  plot_theme_choice <- function(plot_theme = "bw", base_size = 12, font_family = "serif") {
    if (is.null(plot_theme) || !nzchar(plot_theme)) plot_theme <- "bw"
    if (is.null(font_family) || !nzchar(font_family)) font_family <- "serif"
    switch(tolower(plot_theme),
      classic = ggplot2::theme_classic(base_size = base_size, base_family = font_family),
      minimal = ggplot2::theme_minimal(base_size = base_size, base_family = font_family),
      linedraw = ggplot2::theme_linedraw(base_size = base_size, base_family = font_family),
      light = ggplot2::theme_light(base_size = base_size, base_family = font_family),
      gray = ggplot2::theme_gray(base_size = base_size, base_family = font_family),
      dark = ggplot2::theme_dark(base_size = base_size, base_family = font_family),
      void = ggplot2::theme_void(base_size = base_size, base_family = font_family),
      bw = ggplot2::theme_bw(base_size = base_size, base_family = font_family),
      ggplot2::theme_bw(base_size = base_size, base_family = font_family)
    )
  }
}

find_kegg_species_code <- function(query) {
  query <- trimws(as.character(query))
  if (!nzchar(query)) return(NA_character_)
  hits <- tryCatch(KEGGREST::keggFind("organism", query), error = function(e) NULL)
  if (is.null(hits) || length(hits) == 0) return(NA_character_)
  first <- strsplit(names(hits)[1], ":", fixed = TRUE)[[1]]
  if (length(first) >= 2) first[2] else names(hits)[1]
}

########################################################
# 1. Fetch & Cache KEGG Genes
########################################################
get_kegg_genes_cached <- function(kegg_species = "ath", cache_file = NULL) {
  if (is.null(cache_file)) {
    cache_file <- file.path("kegg_genes_cache", paste0("kegg_", kegg_species, "_genes_cache.rds"))
  }
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  
  message("Downloading KEGG pathways for ", kegg_species, ". This may take a few minutes...")
  pathways.list <- KEGGREST::keggList("pathway", kegg_species)
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
  pathway.names <- sub(" - .*", "", pathway.names)
  names(pathway.names) <- pathway.codes
  
  res <- list(
    genes_by_pathway = genes.by.pathway,
    pathway_names = pathway.names
  )
  
  # Try to save to cache
  tryCatch(saveRDS(res, cache_file), error = function(e) warning("Could not save KEGG cache: ", e$message))
  return(res)
}

map_de_ids_for_kegg <- function(de_df, gene_id_type = NULL, orgdb = NULL) {
  keytype <- toupper(gene_id_type %||% "")
  if (keytype %in% c("", "ENTREZID", "TAIR")) return(de_df)
  if (is.null(orgdb) || !nzchar(orgdb)) return(de_df)
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace(orgdb, quietly = TRUE)) return(de_df)

  orgdb_obj <- get(orgdb, envir = asNamespace(orgdb))
  if (!keytype %in% AnnotationDbi::keytypes(orgdb_obj) || !"ENTREZID" %in% AnnotationDbi::columns(orgdb_obj)) {
    return(de_df)
  }

  lookup <- data.frame(
    original_gene_id = as.character(de_df$gene_id),
    lookup_id = trimws(as.character(de_df$gene_id)),
    stringsAsFactors = FALSE
  )
  if (keytype == "ENSEMBL") lookup$lookup_id <- sub("\\.[0-9]+$", "", lookup$lookup_id)
  lookup <- lookup[!is.na(lookup$lookup_id) & lookup$lookup_id != "", , drop = FALSE]
  if (nrow(lookup) == 0) return(de_df)

  mapped <- tryCatch(
    suppressMessages(AnnotationDbi::select(
      orgdb_obj,
      keys = unique(lookup$lookup_id),
      keytype = keytype,
      columns = "ENTREZID"
    )),
    error = function(e) NULL
  )
  if (is.null(mapped) || nrow(mapped) == 0 || !"ENTREZID" %in% names(mapped)) return(de_df)

  mapped <- mapped[!is.na(mapped$ENTREZID) & mapped$ENTREZID != "", , drop = FALSE]
  if (nrow(mapped) == 0) return(de_df)

  lookup_mapped <- dplyr::inner_join(
    lookup,
    mapped,
    by = stats::setNames(keytype, "lookup_id")
  )
  lookup_mapped <- lookup_mapped[!duplicated(lookup_mapped$original_gene_id), , drop = FALSE]
  if (nrow(lookup_mapped) == 0) return(de_df)

  out <- dplyr::inner_join(
    de_df,
    lookup_mapped[, c("original_gene_id", "ENTREZID"), drop = FALSE],
    by = c("gene_id" = "original_gene_id")
  )
  if (nrow(out) == 0) return(de_df)
  out$original_gene_id <- out$gene_id
  out$gene_id <- as.character(out$ENTREZID)
  out$ENTREZID <- NULL
  out
}


########################################################
# 2. Run Wilcoxon Enrichment Test
########################################################
run_kegg_enrichment <- function(de_df, pvalue_cutoff = 0.05, kegg_species = "ath",
                                gene_id_type = NULL, orgdb = NULL) {
  # Requires gene_id, padj/pValue, log2FoldChange
  if (!"pValue" %in% names(de_df)) de_df$pValue <- de_df$padj
  original_n <- nrow(de_df)
  original_ids <- unique(as.character(de_df$gene_id))
  de_df <- map_de_ids_for_kegg(de_df, gene_id_type = gene_id_type, orgdb = orgdb)
  mapped_n <- length(unique(as.character(de_df$gene_id)))
  
  # Remove NAs
  de_df <- de_df[!is.na(de_df$gene_id) & !is.na(de_df$pValue), ]
  
  # Ensure pValues > 0 for log math if needed, and cap at 0.999 for non-sig
  geneList <- de_df$pValue
  names(geneList) <- de_df$gene_id
  
  # Also store log2FC for up/down classification
  lfcList <- de_df$log2FoldChange
  names(lfcList) <- de_df$gene_id
  
  kegg_data <- get_kegg_genes_cached(kegg_species = kegg_species)
  genes.by.pathway <- kegg_data$genes_by_pathway
  pathways.list <- kegg_data$pathway_names
  kegg_gene_universe <- unique(unlist(genes.by.pathway, use.names = FALSE))
  kegg_gene_universe <- kegg_gene_universe[!is.na(kegg_gene_universe)]
  overlap_n <- length(intersect(unique(names(geneList)), kegg_gene_universe))
  if (overlap_n == 0) {
    stop(
      "No overlap between the loaded gene IDs and KEGG ", kegg_species, " pathway genes. ",
      "For human KEGG this usually requires Entrez IDs; selected Gene ID type was '", gene_id_type %||% "not set", "'. ",
      "Mapped IDs tested: ", mapped_n, " from ", length(original_ids), " input IDs."
    )
  }
  
  # Wilcoxon test for each pathway
  res_list <- lapply(names(genes.by.pathway), function(pathway) {
    pathway.genes <- genes.by.pathway[[pathway]]
    if (any(is.na(pathway.genes))) return(NULL)
    
    list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
    if (length(list.genes.in.pathway) < 3) return(NULL) # Skip very small pathways
    
    list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
    if (length(list.genes.not.in.pathway) < 3) return(NULL)
    
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
  outdat$Gene_ID_type <- gene_id_type %||% ""
  outdat$Mapped_gene_ids <- mapped_n
  outdat$KEGG_overlap_gene_ids <- overlap_n
  outdat$Input_rows <- original_n
  rownames(outdat) <- NULL
  return(outdat)
}

########################################################
# 3. Bubble Plot for Enriched Pathways
########################################################
plot_kegg_bubble <- function(kegg_res_df, p_value_threshold = 0.05,
                             color_up = "#B2182B", color_down = "#2166AC",
                             plot_theme = "classic", font_family = "serif") {
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
      scale_color_gradient("p.value", low = color_up, high = "#979797") +
      plot_theme_choice(plot_theme, base_size = 12, font_family = font_family) +
      theme(
        legend.key.size = unit(0.25, 'cm'),
        legend.title = element_text(size = 9.5)
      ) +
      labs(title = "Upregulated", x = "Significant Genes", y = "") +
      guides(color = guide_colorbar(order = 1, barheight = 4))
  }
  
  down_plot <- NULL
  if (nrow(down_df) > 0) {
    down_plot <- ggplot(down_df, aes(x = Downregulated, y = reorder(pathway.name, Downregulated), color = p.value, size = Annotated)) +
      geom_point() +
      scale_color_gradient("p.value", low = color_down, high = "#979797") +
      plot_theme_choice(plot_theme, base_size = 12, font_family = font_family) +
      theme(
        legend.key.size = unit(0.25, 'cm'),
        legend.title = element_text(size = 9.5)
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
