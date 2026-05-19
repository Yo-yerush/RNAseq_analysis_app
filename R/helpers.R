# Helper functions for the local RNA-seq dashboard app
# Adapted from the user's existing DESeq2/GO/REVIGO scripts into app-safe functions.

clean_tair_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\..*$", "", x)
  toupper(x)
}

as_numeric_p <- function(x) {
  if (is.numeric(x)) return(x)
  x <- trimws(as.character(x))
  x <- gsub(",", "", x)
  x <- gsub("^<\\s*", "", x)
  x <- gsub("^>\\s*", "", x)
  suppressWarnings(as.numeric(x))
}

first_existing_col <- function(df, candidates) {
  nms <- names(df)
  nms_low <- tolower(nms)
  cand_low <- tolower(candidates)
  hit <- match(cand_low, nms_low)
  hit <- hit[!is.na(hit)]
  if (length(hit) == 0) return(NULL)
  nms[hit[1]]
}

read_any_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv", "txt", "tsv")) {
    if (requireNamespace("readr", quietly = TRUE)) {
      if (ext == "csv") {
        return(readr::read_csv(path, show_col_types = FALSE, guess_max = 100000))
      } else {
        return(readr::read_tsv(path, show_col_types = FALSE, guess_max = 100000))
      }
    }
  }
  if (ext == "csv") {
    utils::read.csv(path, check.names = FALSE)
  } else {
    utils::read.table(path, header = TRUE, sep = "\t", check.names = FALSE, quote = "", comment.char = "")
  }
}

# Cached description table – loaded once per session
.desc_cache <- NULL
load_description_file <- function(
    path = file.path("description_files", "Methylome.At_description_file.csv.gz")
) {
  if (!is.null(.desc_cache)) return(.desc_cache)
  if (!file.exists(path)) {
    message("Description file not found: ", path, " — skipping annotation merge.")
    return(NULL)
  }
  tryCatch({
    desc     <- read.csv(path, check.names = FALSE)
    gene_col <- if ("gene_id" %in% names(desc)) "gene_id" else names(desc)[1]
    names(desc)[names(desc) == gene_col] <- "gene_id"
    desc$gene_id <- clean_tair_id(desc$gene_id)
    desc <- desc[!duplicated(desc$gene_id), , drop = FALSE]
    .desc_cache <<- desc
    desc
  }, error = function(e) {
    message("Could not load description file: ", e$message)
    NULL
  })
}

# Left-join DE table with description file.
# Adds annotation columns (GO terms, Symbol, KEGG, etc.) that are NOT already present.
merge_with_description <- function(de_df) {
  desc <- load_description_file()
  if (is.null(desc)) return(de_df)
  new_cols <- setdiff(names(desc), c(names(de_df)))
  if (length(new_cols) == 0) return(de_df)
  desc_sub <- desc[, c("gene_id", new_cols), drop = FALSE]
  merged   <- merge(de_df, desc_sub, by = "gene_id", all.x = TRUE, sort = FALSE)
  # Restore original row order
  idx <- match(de_df$gene_id, merged$gene_id)
  merged[idx[!is.na(idx)], , drop = FALSE]
}

standardize_de_table <- function(df) {
  df <- as.data.frame(df, check.names = FALSE)


  gene_col <- first_existing_col(df, c("gene_id", "gene", "Gene", "GeneID", "TAIR", "locus_tag", "row.names", "rownames", "id"))
  lfc_col  <- first_existing_col(df, c("log2FoldChange", "log2FC", "logFC", "LFC", "fold_change_log2"))
  padj_col <- first_existing_col(df, c("padj", "FDR", "qvalue", "q_value", "adj.P.Val", "p.adjust", "p_adj"))
  p_col    <- first_existing_col(df, c("pValue", "pvalue", "p.value", "P.Value", "pval"))
  base_col <- first_existing_col(df, c("baseMean", "mean", "meanExpression", "MeanExpression"))

  if (is.null(gene_col) || is.null(lfc_col) || is.null(padj_col)) {
    if (ncol(df) >= 3) {
      names(df)[1:3] <- c("gene_id", "log2FoldChange", "padj")
      gene_col <- "gene_id"
      lfc_col  <- "log2FoldChange"
      padj_col <- "padj"
    } else {
      stop("Could not find required columns and the input has fewer than 3 columns.")
    }
  }

  if (gene_col != "gene_id") names(df)[names(df) == gene_col] <- "gene_id"
  if (lfc_col != "log2FoldChange") names(df)[names(df) == lfc_col] <- "log2FoldChange"
  if (padj_col != "padj") names(df)[names(df) == padj_col] <- "padj"
  if (!is.null(p_col) && p_col != "pValue") names(df)[names(df) == p_col] <- "pValue"
  if (!is.null(base_col) && base_col != "baseMean") names(df)[names(df) == base_col] <- "baseMean"

  # Drop heavy extra columns from uploaded CSVs to save memory
  keep_cols <- c("gene_id", "log2FoldChange", "padj", "pValue", "baseMean")
  df <- df[, intersect(keep_cols, names(df)), drop = FALSE]

  df$gene_id <- clean_tair_id(df$gene_id)
  df$log2FoldChange <- suppressWarnings(as.numeric(df$log2FoldChange))
  df$padj <- as_numeric_p(df$padj)
  if (!"pValue" %in% names(df)) df$pValue <- df$padj
  df$pValue <- as_numeric_p(df$pValue)
  if ("baseMean" %in% names(df)) df$baseMean <- suppressWarnings(as.numeric(df$baseMean))
  df <- df[!is.na(df$gene_id) & df$gene_id != "", , drop = FALSE]
  # Merge with description file to add GO terms, Symbol, KEGG, etc.
  df <- merge_with_description(df)
  df
}

classify_de <- function(df, alpha = 0.05, lfc_cutoff = 1) {
  df$DE_class <- "not_significant"
  df$DE_class[!is.na(df$padj) & df$padj < alpha & df$log2FoldChange >= lfc_cutoff] <- "up"
  df$DE_class[!is.na(df$padj) & df$padj < alpha & df$log2FoldChange <= -lfc_cutoff] <- "down"
  df
}

make_volcano_plot <- function(df, alpha = 0.05, lfc_cutoff = 1, title = "Volcano plot",
                              point_size = 1, point_alpha = 0.65,
                              color_up = "#B2182B", color_down = "#2166AC", color_ns = "grey70") {
  df <- classify_de(df, alpha, lfc_cutoff)
  df$neg_log10_padj <- -log10(pmax(df$padj, .Machine$double.xmin))
  ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange, y = neg_log10_padj, color = DE_class)) +
    ggplot2::geom_point(alpha = point_alpha, size = point_size) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed", linewidth = 0.3) +
    ggplot2::scale_color_manual(values = c(up = color_up, down = color_down, not_significant = color_ns)) +
    ggplot2::theme_classic(base_size = 12, base_family = "serif") +
    ggplot2::labs(x = "log2 fold change", y = "-log10 adjusted p-value", color = "Class", title = title)
}

make_ma_plot <- function(df, alpha = 0.05, lfc_cutoff = 1, title = "MA plot",
                         point_size = 1, point_alpha = 0.65,
                         color_up = "#B2182B", color_down = "#2166AC", color_ns = "grey70") {
  if (!"baseMean" %in% names(df)) {
    stop("A true MA plot requires a baseMean/mean-expression column. It is available after running DESeq2 from RSEM files, or in uploaded CSVs that include baseMean.")
  }
  df <- df[!is.na(df$baseMean) & df$baseMean > 0, , drop = FALSE]
  df <- classify_de(df, alpha, lfc_cutoff)
  ns_df <- df[df$DE_class == "not_significant", , drop = FALSE]
  sig_df <- df[df$DE_class %in% c("up", "down"), , drop = FALSE]
  ggplot2::ggplot(df, ggplot2::aes(x = baseMean, y = log2FoldChange, color = DE_class)) +
    ggplot2::geom_point(data = ns_df, alpha = point_alpha, size = point_size) +
    ggplot2::geom_point(data = sig_df, alpha = point_alpha, size = point_size) +
    ggplot2::geom_hline(yintercept = c(-lfc_cutoff, 0, lfc_cutoff), linetype = c("dashed", "solid", "dashed"), linewidth = 0.3) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = c(up = color_up, down = color_down, not_significant = color_ns)) +
    ggplot2::theme_classic(base_size = 12, base_family = "serif") +
    ggplot2::labs(x = "Mean of normalized counts", y = "log2 fold change", color = "Class", title = title)
}

make_pca_plot <- function(pca_df, title = "PCA", point_size = 3.2, point_alpha = 0.9,
                          show_labels = TRUE, conditions = NULL) {
  req_cols <- c("PC1", "PC2", "condition", "name")
  if (!all(req_cols %in% names(pca_df))) stop("PCA table is missing required columns.")
  # Filter to selected conditions if provided
  if (!is.null(conditions) && length(conditions) > 0) {
    pca_df <- pca_df[pca_df$condition %in% conditions, , drop = FALSE]
  }
  if (nrow(pca_df) == 0) stop("No samples match the selected conditions.")
  # Use sample_label for display if available, otherwise fall back to name
  label_col <- if ("sample_label" %in% names(pca_df) && !all(is.na(pca_df$sample_label))) "sample_label" else "name"
  pca_df$display_label <- pca_df[[label_col]]
  pc1lab <- if ("percentVar1" %in% names(pca_df)) paste0("PC1: ", unique(pca_df$percentVar1)[1], "% variance") else "PC1"
  pc2lab <- if ("percentVar2" %in% names(pca_df)) paste0("PC2: ", unique(pca_df$percentVar2)[1], "% variance") else "PC2"
  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = condition, label = display_label)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey70", linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, color = "grey70", linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha)
  if (isTRUE(show_labels)) {
    p <- p + ggrepel::geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = 30)
  }
  p +
    ggplot2::theme_classic(base_size = 12, base_family = "serif") +
    ggplot2::labs(x = pc1lab, y = pc2lab, color = "Condition", title = title)
}

normalize_coldata <- function(coldata, sample_col = NULL, condition_col = NULL, label_col = NULL) {
  coldata <- as.data.frame(coldata, check.names = FALSE)
  if (is.null(sample_col)) {
    sample_col <- first_existing_col(coldata, c("sample_id", "x", "file", "SampleID", "sample_name", "sample"))
  }
  if (is.null(condition_col)) {
    condition_col <- first_existing_col(coldata, c("condition", "exp", "group", "genotype", "treatment"))
  }
  if (is.null(label_col)) {
    label_col <- first_existing_col(coldata, c("sample_label", "label", "sample", "Sample"))
  }
  if (is.null(sample_col)) stop("Could not find a sample/file column in colData.")
  if (is.null(condition_col)) stop("Could not find a condition/group column in colData.")
  out <- data.frame(
    sample_id = as.character(coldata[[sample_col]]),
    condition = as.character(coldata[[condition_col]]),
    stringsAsFactors = FALSE
  )
  if (!is.null(label_col) && label_col %in% names(coldata)) {
    out$sample_label <- as.character(coldata[[label_col]])
  } else {
    out$sample_label <- out$sample_id
  }
  out$sample_id <- sub("\\.genes\\.results$", "", out$sample_id)
  out$sample_id <- basename(out$sample_id)
  out <- out[!is.na(out$sample_id) & out$sample_id != "", , drop = FALSE]
  out
}

scan_rsem_files <- function(folder) {
  if (is.null(folder) || !dir.exists(folder)) return(data.frame())
  files <- list.files(folder, pattern = "\\.genes\\.results$", full.names = TRUE)
  data.frame(
    sample_id = sub("\\.genes\\.results$", "", basename(files)),
    condition = "condition_1",
    sample_label = sub("\\.genes\\.results$", "", basename(files)),
    file = files,
    stringsAsFactors = FALSE
  )
}

run_deseq2_from_rsem <- function(folder, coldata, treatment, control, lfc_shrink = FALSE, min_count = 10) {
  required <- c("DESeq2", "tximport", "tibble")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "))

  rsem_table <- scan_rsem_files(folder)
  if (nrow(rsem_table) == 0) stop("No .genes.results files were found in the selected folder.")

  coldata <- normalize_coldata(coldata, sample_col = "sample_id", condition_col = "condition", label_col = "sample_label")
  coldata$sample_id <- sub("\\.genes\\.results$", "", basename(coldata$sample_id))
  rownames(coldata) <- coldata$sample_id
  rsem_table <- rsem_table[match(coldata$sample_id, rsem_table$sample_id), , drop = FALSE]
  if (any(is.na(rsem_table$file))) {
    missing_samples <- coldata$sample_id[is.na(rsem_table$file)]
    stop("These colData sample_id values do not match .genes.results files: ", paste(missing_samples, collapse = ", "))
  }

  if (!treatment %in% coldata$condition) stop("Treatment is not present in colData condition column.")
  if (!control %in% coldata$condition) stop("Control is not present in colData condition column.")

  coldata$condition <- factor(coldata$condition)
  coldata$condition <- stats::relevel(coldata$condition, ref = control)
  files <- rsem_table$file
  names(files) <- coldata$sample_id

  txi <- tximport::tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
  txi$length[txi$length == 0] <- 1
  dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
  keep <- rowSums(DESeq2::counts(dds)) >= min_count
  dds <- dds[keep, ]
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("condition", treatment, control), alpha = 0.05)

  if (isTRUE(lfc_shrink)) {
    if (requireNamespace("ashr", quietly = TRUE)) {
      res <- DESeq2::lfcShrink(dds, contrast = c("condition", treatment, control), res = res, type = "ashr")
    } else {
      warning("Package 'ashr' is not installed, so lfcShrink was skipped.")
    }
  }

  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- tibble::as_tibble(res_df)
  if ("pvalue" %in% names(res_df)) names(res_df)[names(res_df) == "pvalue"] <- "pValue"
  res_df <- res_df[, c("gene_id", setdiff(names(res_df), "gene_id"))]
  res_df <- standardize_de_table(res_df)
  res_df <- res_df[order(res_df$padj), , drop = FALSE]

  norm_counts <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
  norm_counts$gene_id <- clean_tair_id(rownames(norm_counts))
  norm_counts <- norm_counts[, c("gene_id", setdiff(names(norm_counts), "gene_id"))]

  vst_obj <- DESeq2::vst(dds, blind = FALSE)
  pca_data <- DESeq2::plotPCA(vst_obj, intgroup = "condition", returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"), 1)
  pca_data$name <- rownames(pca_data)
  pca_data$sample_label <- coldata$sample_label[match(pca_data$name, coldata$sample_id)]
  pca_data$percentVar1 <- percentVar[1]
  pca_data$percentVar2 <- percentVar[2]

  list(
    de_table = res_df,
    norm_counts = norm_counts,
    pca_table = pca_data,
    summary = paste(capture.output(summary(res)), collapse = "\n"),
    coldata = coldata
  )
}

detect_go_bp_col <- function(df) {
  nms <- names(df)
  low <- tolower(nms)
  exact <- which(low %in% c("go.biological.process", "gene.ontology..biological.process.", "gene_ontology_biological_process", "go_bp", "go.biological_process"))
  if (length(exact)) return(nms[exact[1]])
  hit <- grep("go.*biological.*process|biological.*process|go_bp", low)
  if (length(hit)) return(nms[hit[1]])
  NULL
}

add_go_bp_column <- function(df) {
  bp_col <- detect_go_bp_col(df)
  if (!is.null(bp_col)) {
    df$GO_BP_terms <- as.character(df[[bp_col]])
    return(df)
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace("org.At.tair.db", quietly = TRUE)) {
    df$GO_BP_terms <- NA_character_
    return(df)
  }
  genes <- unique(clean_tair_id(df$gene_id))
  annot <- AnnotationDbi::select(org.At.tair.db::org.At.tair.db, keys = genes, keytype = "TAIR", columns = c("GO", "ONTOLOGY"))
  annot <- annot[!is.na(annot$GO) & annot$ONTOLOGY == "BP", , drop = FALSE]
  collapsed <- stats::aggregate(GO ~ TAIR, data = annot, FUN = function(x) paste(unique(x), collapse = "; "))
  names(collapsed) <- c("gene_id", "GO_BP_terms")
  df <- merge(df, collapsed, by = "gene_id", all.x = TRUE, sort = FALSE)
  df
}

run_topgo_enrichment <- function(de_df, direction = c("up", "down", "all"), ontology = "BP", alpha = 0.05, lfc_cutoff = 1,
                                 p_cutoff = 0.01, algorithm = "weight01", statistic = "fisher") {
  direction <- match.arg(direction)
  required <- c("topGO", "GO.db", "org.At.tair.db")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages(library(topGO))
  suppressPackageStartupMessages(library(GO.db))
  de_df <- classify_de(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  bg <- unique(clean_tair_id(de_df$gene_id))
  if (direction == "up") interesting <- unique(de_df$gene_id[de_df$DE_class == "up"])
  if (direction == "down") interesting <- unique(de_df$gene_id[de_df$DE_class == "down"])
  if (direction == "all") interesting <- unique(de_df$gene_id[de_df$DE_class %in% c("up", "down")])
  interesting <- intersect(clean_tair_id(interesting), bg)
  if (length(interesting) < 2) stop("Too few significant genes for GO enrichment in this direction.")

  geneList <- factor(as.integer(bg %in% interesting))
  names(geneList) <- bg
  GOdata <- methods::new("topGOdata", ontology = ontology, allGenes = geneList,
                         geneSelectionFun = function(x) x == 1,
                         annot = get("annFUN.org", envir = asNamespace("topGO")), mapping = "org.At.tair.db", ID = "entrez")
  result <- topGO::runTest(GOdata, algorithm = algorithm, statistic = statistic)
  all_res <- topGO::GenTable(GOdata, pValue = result, topNodes = length(result@score))
  all_res$pValue_num <- as_numeric_p(all_res$pValue)
  all_res$Significant <- suppressWarnings(as.numeric(all_res$Significant))
  all_res$Expected <- suppressWarnings(as.numeric(all_res$Expected))
  all_res$FoldEnrichment <- all_res$Significant / pmax(all_res$Expected, 1e-12)
  all_res$Direction <- direction
  all_res$Ontology <- ontology
  all_res <- all_res[order(all_res$pValue_num), , drop = FALSE]
  sig <- all_res[!is.na(all_res$pValue_num) & all_res$pValue_num <= p_cutoff, , drop = FALSE]
  sig
}

make_go_bubble_plot <- function(go_df, title = "GO enrichment", top_n = 20, direction = NULL, point_alpha = 0.82) {
  if (is.null(go_df) || nrow(go_df) == 0) stop("No GO terms to plot.")
  go_df <- go_df[order(go_df$pValue_num), , drop = FALSE]
  go_df <- head(go_df, top_n)
  go_df$Term_short <- stringr::str_trunc(go_df$Term, 75)
  go_df$neg_log10_p <- -log10(pmax(go_df$pValue_num, .Machine$double.xmin))
  high_col <- "#b2182b"
  if (!is.null(direction)) {
    dir_low <- tolower(direction)
    if (dir_low == "down") high_col <- "#2166ac"
    if (dir_low == "all") high_col <- "#000000"
  }
  ggplot2::ggplot(go_df, ggplot2::aes(x = FoldEnrichment, y = stats::reorder(Term_short, FoldEnrichment))) +
    ggplot2::geom_point(ggplot2::aes(size = Significant, color = neg_log10_p), alpha = point_alpha) +
    ggplot2::theme_classic(base_size = 12, base_family = "serif") +
    ggplot2::scale_color_gradient(low = "#d9d9d9", high = high_col) +
    ggplot2::labs(x = "Fold enrichment", y = NULL, size = "Genes", color = "-log10(p)", title = title)
}

go_offspring_terms <- function(go_id) {
  if (!requireNamespace("GO.db", quietly = TRUE)) stop("Package GO.db is required.")
  terms <- as.character(GO.db::GOBPOFFSPRING[[go_id]])
  unique(c(go_id, terms[!is.na(terms)]))
}

go_term_title <- function(go_id) {
  if (!requireNamespace("GO.db", quietly = TRUE)) return(go_id)
  out <- tryCatch(GO.db::GOTERM[[go_id]]@Term, error = function(e) NA_character_)
  ifelse(is.na(out), go_id, out)
}

genes_matching_go_terms <- function(df, go_terms, go_col = "GO_BP_terms") {
  if (!go_col %in% names(df)) return(character())
  txt <- as.character(df[[go_col]])
  txt[is.na(txt)] <- ""
  hits <- grepl(paste(go_terms, collapse = "|"), txt)
  unique(df$gene_id[hits])
}

make_go_offspring_summary <- function(de_df, parent_go_ids, alpha = 0.05, lfc_cutoff = 1) {
  de_df <- add_go_bp_column(de_df)
  de_df <- classify_de(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  parent_go_ids <- trimws(unlist(strsplit(parent_go_ids, ",|;|\\s+")))
  parent_go_ids <- parent_go_ids[parent_go_ids != ""]
  rows <- lapply(parent_go_ids, function(go_id) {
    terms <- go_offspring_terms(go_id)
    genes <- genes_matching_go_terms(de_df, terms)
    sub <- de_df[de_df$gene_id %in% genes, , drop = FALSE]
    sig <- sub[sub$DE_class %in% c("up", "down"), , drop = FALSE]
    data.frame(
      Parent_GO_ID = go_id,
      Category = go_term_title(go_id),
      Offspring_terms = length(terms),
      Total = length(unique(sub$gene_id)),
      Upregulated = length(unique(sig$gene_id[sig$DE_class == "up"])),
      Downregulated = length(unique(sig$gene_id[sig$DE_class == "down"])),
      Significant = length(unique(sig$gene_id)),
      Percentage = ifelse(length(unique(sub$gene_id)) > 0, paste0(round(100 * length(unique(sig$gene_id)) / length(unique(sub$gene_id)), 2), "%"), "NA"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

make_abiotic_stress_table <- function(de_df, dataset = c("all", "up", "down"), alpha = 0.05, lfc_cutoff = 1) {
  dataset <- match.arg(dataset)
  de_df <- add_go_bp_column(de_df)
  de_df <- classify_de(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  bg_total <- unique(de_df$gene_id)
  if (dataset == "all") bg_sig <- unique(de_df$gene_id[de_df$DE_class %in% c("up", "down")])
  if (dataset == "up") bg_sig <- unique(de_df$gene_id[de_df$DE_class == "up"])
  if (dataset == "down") bg_sig <- unique(de_df$gene_id[de_df$DE_class == "down"])

  stress_parents <- c(
    Cold = "GO:0009409",
    Osmotic = "GO:0006970",
    Salt = "GO:1902074",
    Water_deprivation = "GO:0009414",
    DNA_damage = "GO:0006974",
    Oxidative = "GO:0006979",
    UVB = "GO:0010224",
    Wounding = "GO:0009611",
    Heat = "GO:0009408"
  )
  stress_names <- c(
    Cold = "Cold stress", Osmotic = "Osmotic stress", Salt = "Salt stress",
    Water_deprivation = "Water deprivation", DNA_damage = "DNA damage",
    Oxidative = "Oxidative stress", UVB = "UV-B stress", Wounding = "Wounding", Heat = "Heat stress"
  )

  rows <- lapply(names(stress_parents), function(key) {
    terms <- go_offspring_terms(stress_parents[[key]])
    term_genes <- genes_matching_go_terms(de_df, terms)
    a <- length(intersect(term_genes, bg_sig))
    b <- length(setdiff(term_genes, bg_sig))
    c <- length(setdiff(bg_sig, term_genes))
    d <- max(length(bg_total) - a - b - c, 0)
    mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                  dimnames = list(c("In_term", "Not_in_term"), c("Significant", "Not_significant")))
    ft <- tryCatch(stats::fisher.test(mat, alternative = "greater"), error = function(e) NULL)
    bg_rate <- length(bg_sig) / max(length(bg_total), 1)
    term_rate <- a / max(length(term_genes), 1)
    data.frame(
      Stress_type = key,
      Test_name = stress_names[[key]],
      Parent_GO_ID = stress_parents[[key]],
      Sig_in_term = a,
      Total_in_term = length(term_genes),
      Sig_in_background = length(bg_sig),
      Total_background = length(bg_total),
      Fold_enrichment = ifelse(bg_rate > 0, term_rate / bg_rate, NA_real_),
      P_value = if (is.null(ft)) NA_real_ else ft$p.value,
      Odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate),
      Significant = if (is.null(ft)) FALSE else ft$p.value < 0.05,
      Dataset = dataset,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$P_value), , drop = FALSE]
}

make_abiotic_stress_plot <- function(stress_df, title = "Abiotic stress enrichment") {
  if (is.null(stress_df) || nrow(stress_df) == 0) stop("No stress enrichment table to plot.")
  stress_df$Label <- paste0("p=", signif(stress_df$P_value, 2), "\n", stress_df$Sig_in_term, "/", stress_df$Total_in_term)
  ggplot2::ggplot(stress_df, ggplot2::aes(x = stats::reorder(Test_name, Fold_enrichment), y = Fold_enrichment, fill = Significant)) +
    ggplot2::geom_col(alpha = 0.85) +
    ggplot2::geom_text(ggplot2::aes(label = Label), hjust = -0.05, size = 3) +
    ggplot2::coord_flip() +
    ggplot2::theme_classic(base_size = 12, base_family = "serif") +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey70")) +
    ggplot2::labs(x = NULL, y = "Fold enrichment", fill = "p < 0.05", title = title) +
    ggplot2::expand_limits(y = max(stress_df$Fold_enrichment, na.rm = TRUE) * 1.25)
}

run_rrvgo_reduce <- function(go_df, ontology = "BP", top_n = 80, threshold = 0.7, title = "REVIGO-like semantic reduction") {
  required <- c("rrvgo", "org.At.tair.db", "ggplot2")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "))
  if (is.null(go_df) || nrow(go_df) < 2) stop("At least two GO terms are required for semantic reduction.")
  go_df <- go_df[order(go_df$pValue_num), , drop = FALSE]
  go_df <- head(go_df, top_n)
  scores <- stats::setNames(-log10(pmax(go_df$pValue_num, .Machine$double.xmin)), go_df$GO.ID)
  simMatrix <- rrvgo::calculateSimMatrix(names(scores), orgdb = "org.At.tair.db", ont = ontology, method = "Rel")
  reducedTerms <- rrvgo::reduceSimMatrix(simMatrix, scores, threshold = threshold, orgdb = "org.At.tair.db")
  plot_obj <- tryCatch({
    rrvgo::scatterPlot(simMatrix, reducedTerms, algorithm = "umap") + ggplot2::ggtitle(title)
  }, error = function(e) {
    rrvgo::scatterPlot(simMatrix, reducedTerms, algorithm = "pca") + ggplot2::ggtitle(paste0(title, " (PCA layout)"))
  })
  list(plot = plot_obj, table = as.data.frame(reducedTerms))
}

save_plot_png <- function(plot_obj, file, width = 9, height = 6, dpi = 300) {
  ggplot2::ggsave(filename = file, plot = plot_obj, width = width, height = height, dpi = dpi, bg = "white")
}

render_svg_plot <- function(plot_obj, w, h) {
  if (is.null(plot_obj)) return(NULL)
  tf <- tempfile(fileext = ".svg")
  on.exit(unlink(tf), add = TRUE)
  tryCatch({
    ggplot2::ggsave(tf, plot = plot_obj, width = max(w, 100)/100, height = max(h, 100)/100, device = "svg", bg = "white")
    svg_code <- paste(readLines(tf, warn = FALSE), collapse = "\n")
    shiny::HTML(svg_code)
  }, error = function(e) {
    shiny::tags$div(style = "color: red; padding: 20px;", paste("Error rendering SVG:", e$message))
  })
}

download_plot_ui <- function(id, default_btn_label = "Download") {
  shiny::div(class = "download-row", style = "display: flex; align-items: center; gap: 10px; margin-top: 6px;",
    shiny::selectInput(paste0("format_", id), NULL, choices = c("SVG" = "svg", "PNG" = "png", "PDF" = "pdf"), width = "80px", selectize = FALSE),
    shiny::downloadButton(paste0("download_", id), default_btn_label)
  )
}

download_plot_server <- function(plot_reactive, format_reactive, filename_prefix, w_input, h_input) {
  shiny::downloadHandler(
    filename = function() {
      prefix <- if (is.function(filename_prefix)) filename_prefix() else filename_prefix
      paste0(prefix, "_", Sys.Date(), ".", format_reactive())
    },
    content = function(file) {
      shiny::req(plot_reactive())
      w_in <- max(w_input(), 100) / 96
      h_in <- max(h_input(), 100) / 96
      fmt  <- format_reactive()
      if (fmt == "png") {
        ggplot2::ggsave(file, plot = plot_reactive(),
                        width = w_in, height = h_in, units = "in", dpi = 96, bg = "white")
      } else if (fmt == "svg") {
        grDevices::svg(filename = file, width = w_in, height = h_in, bg = "white")
        print(plot_reactive())
        grDevices::dev.off()
      } else if (fmt == "pdf") {
        grDevices::pdf(file = file, width = w_in, height = h_in, bg = "white")
        print(plot_reactive())
        grDevices::dev.off()
      }
    }
  )
}
