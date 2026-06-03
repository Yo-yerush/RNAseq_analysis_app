# Helper functions for the local RNA-seq dashboard app
# Adapted from the user's existing DESeq2/GO/REVIGO scripts into app-safe functions.

`%||%` <- function(a, b) if (!is.null(a)) a else b

clean_tair_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\..*$", "", x)
  toupper(x)
}

clean_gene_id <- function(x) {
  x <- trimws(as.character(x))
  sub("\\.[0-9]+$", "", x)
}

gene_join_key <- function(x) {
  x <- clean_gene_id(x)
  x <- sub("^.*:", "", x)
  toupper(x)
}

normalize_orgdb_gene_keys <- function(x, orgdb = NULL, keytype = NULL) {
  keytype <- toupper(keytype %||% "")
  out <- clean_gene_id(x)

  # org.EcK12.eg.db stores most E. coli K-12 Blattner b-numbers as ECK aliases.
  # RefSeq accessions such as NP_414542 are handled by REFSEQ/ACCNUM directly.
  if (identical(orgdb, "org.EcK12.eg.db") && identical(keytype, "ALIAS")) {
    hit <- grepl("^b[0-9]{4}$", out, ignore.case = TRUE)
    out[hit] <- paste0("ECK", sub("^[bB]", "", out[hit]))
  }

  out
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

count_fixed_matches <- function(x, pattern) {
  matches <- gregexpr(pattern, x, fixed = TRUE)[[1]]
  if (length(matches) == 1 && matches[1] == -1L) 0L else length(matches)
}

is_remote_path <- function(path) {
  grepl("^(https?|ftp)://", path, ignore.case = TRUE)
}

DEFAULT_ARABIDOPSIS_DESCRIPTION_URL <- "https://github.com/Yo-yerush/RA_lab_db/raw/refs/heads/main/description_files/At_custom_description_file.csv.gz"
DEFAULT_HUMAN_DESCRIPTION_URL <- "https://github.com/Yo-yerush/RA_lab_db/raw/refs/heads/main/description_files/Hs_description_file.csv.gz"
DEFAULT_MG1655_DESCRIPTION_URL <- "https://raw.githubusercontent.com/Yo-yerush/RA_lab_db/refs/heads/main/description_files/MG1655_description_file.csv"
DEFAULT_TAIR_TE_URL <- "https://raw.githubusercontent.com/Yo-yerush/RA_lab_db/refs/heads/main/description_files/TAIR10_Transposable_Elements.txt"
DEFAULT_AT_GENE_FAMILIES_URL <- "https://raw.githubusercontent.com/Yo-yerush/RA_lab_db/refs/heads/main/description_files/gene_families_sep_29_09_update.txt"
DEFAULT_HGNC_FAMILY_URL <- "https://storage.googleapis.com/public-download-files/hgnc/csv/csv/genefamily_db_tables/family.csv"
DEFAULT_HGNC_GENE_HAS_FAMILY_URL <- "https://storage.googleapis.com/public-download-files/hgnc/csv/csv/genefamily_db_tables/gene_has_family.csv"
DEFAULT_HGNC_COMPLETE_SET_URL <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
PMN_PATHWAYS_BASE_URL <- "https://plantcyc-ftp.storage.googleapis.com/pmn/Pathways/Data_dumps/PMN15.5_January2023/pathways"
PMN_PATHWAYS_DATE <- "20230103"

default_description_file_url <- function(tax_id = 3702) {
  tax_id <- suppressWarnings(as.integer(tax_id))
  if (identical(tax_id, 3702L)) return(DEFAULT_ARABIDOPSIS_DESCRIPTION_URL)
  if (identical(tax_id, 9606L)) return(DEFAULT_HUMAN_DESCRIPTION_URL)
  if (identical(tax_id, 511145L)) return(DEFAULT_MG1655_DESCRIPTION_URL)
  NULL
}

infer_table_delimiter <- function(path) {
  path_low <- tolower(sub("[?#].*$", "", path))
  ext <- tolower(tools::file_ext(path_low))
  if (grepl("\\.csv\\.gz$", path_low) || ext == "csv") return(",")
  if (grepl("\\.(tsv|txt)\\.gz$", path_low) || ext %in% c("tsv", "txt")) return("\t")
  "\t"
}

detect_table_delimiter <- function(path) {
  con <- if (grepl("\\.gz$", tolower(path))) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, n = 20, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) return("\t")

  header <- lines[1]
  tab_count <- count_fixed_matches(header, "\t")
  comma_count <- count_fixed_matches(header, ",")
  if (tab_count == 0 && comma_count == 0) return(NULL)
  if (comma_count > tab_count) "," else "\t"
}

read_delimited_table <- function(path, delim) {
  if (requireNamespace("readr", quietly = TRUE)) {
    return(readr::read_delim(path, delim = delim, show_col_types = FALSE, guess_max = 100000))
  }
  if (is_remote_path(path)) {
    con <- if (grepl("\\.gz([?#].*)?$", tolower(path))) gzcon(url(path, open = "rb")) else url(path, open = "rt")
    on.exit(close(con), add = TRUE)
  } else if (grepl("\\.gz$", tolower(path))) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
  } else {
    con <- path
  }
  utils::read.table(con, header = TRUE, sep = delim, check.names = FALSE, quote = "\"", comment.char = "")
}

read_excel_table <- function(path) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Missing required package: readxl. Run install_packages.R, then restart the app.")
  }
  readxl::read_excel(path)
}

read_any_table <- function(path, delim = NULL, source_name = path) {
  ext <- tolower(tools::file_ext(source_name))
  path_low <- tolower(source_name)
  if (ext %in% c("xlsx", "xls")) {
    return(read_excel_table(path))
  }
  if (is.null(delim)) {
    delim <- if (is_remote_path(path)) infer_table_delimiter(path) else detect_table_delimiter(path)
  }
  if (is.null(delim)) {
    delim <- infer_table_delimiter(path)
  }
  if (grepl("\\.csv\\.gz$", path_low)) {
    return(read_delimited_table(path, delim))
  }
  if (grepl("\\.(tsv|txt)\\.gz$", path_low)) {
    return(read_delimited_table(path, delim))
  }
  if (ext %in% c("csv", "txt", "tsv")) {
    return(read_delimited_table(path, delim))
  }
  read_delimited_table(path, delim)
}

# Cached description table – loaded once per session
.desc_cache <- NULL
.desc_cache_path <- NULL

normalize_annotation_table <- function(desc) {
  desc <- as.data.frame(desc, check.names = FALSE)
  gene_col <- first_existing_col(desc, c(
    "gene_id", "gene", "Gene", "GeneID", "Gene.ID", "TAIR", "locus_tag",
    "Locus", "id", "ID", "query", "input_id"
  ))
  if (is.null(gene_col)) {
    gene_col <- names(desc)[1]
  }
  names(desc)[names(desc) == gene_col] <- "gene_id"
  desc$gene_id <- clean_gene_id(desc$gene_id)
  desc$annotation_key <- gene_join_key(desc$gene_id)
  desc <- desc[!is.na(desc$annotation_key) & desc$annotation_key != "", , drop = FALSE]
  desc <- desc[!duplicated(desc$annotation_key), , drop = FALSE]
  desc
}

annotation_id_lookup_columns <- function(desc) {
  if (is.null(desc) || ncol(desc) == 0) return(character())
  id_patterns <- paste(c(
    "^gene_id$",
    "gene.?id",
    "refseq",
    "genbank",
    "locus",
    "tag",
    "protein.?id",
    "transcript.?id",
    "uniprot",
    "geneid",
    "entrez",
    "ensembl",
    "alias",
    "symbol"
  ), collapse = "|")
  cols <- names(desc)[grepl(id_patterns, names(desc), ignore.case = TRUE)]
  setdiff(unique(c("gene_id", cols)), "annotation_key")
}

make_annotation_lookup <- function(desc, id_cols = NULL) {
  desc <- normalize_annotation_table(desc)
  if (is.null(id_cols)) id_cols <- annotation_id_lookup_columns(desc)
  id_cols <- intersect(id_cols, names(desc))
  if (!length(id_cols)) id_cols <- "gene_id"

  rows <- lapply(id_cols, function(col) {
    values <- as.character(desc[[col]])
    values[is.na(values)] <- ""
    pieces <- strsplit(values, "[;,|[:space:]]+")
    keys <- unlist(pieces, use.names = FALSE)
    idx <- rep(seq_len(nrow(desc)), lengths(pieces))
    keys <- trimws(keys)
    keep <- !is.na(keys) & nzchar(keys) & keys != "-"
    if (!any(keep)) return(NULL)
    data.frame(
      lookup_key = gene_join_key(keys[keep]),
      annotation_key = desc$annotation_key[idx[keep]],
      annotation_gene_id = desc$gene_id[idx[keep]],
      annotation_id_source_col = col,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0) {
    return(data.frame(
      lookup_key = character(),
      annotation_key = character(),
      annotation_gene_id = character(),
      annotation_id_source_col = character()
    ))
  }
  out <- out[!is.na(out$lookup_key) & out$lookup_key != "", , drop = FALSE]
  out <- out[!duplicated(out$lookup_key), , drop = FALSE]
  rownames(out) <- NULL
  out
}

load_description_file <- function(
    path = DEFAULT_ARABIDOPSIS_DESCRIPTION_URL
) {
  if (!is.null(.desc_cache) && identical(.desc_cache_path, path)) return(.desc_cache)
  if (!is_remote_path(path) && !file.exists(path)) {
    message("Description file not found: ", path, " — skipping annotation merge.")
    return(NULL)
  }
  tryCatch({
    desc <- normalize_annotation_table(read_any_table(path))
    .desc_cache <<- desc
    .desc_cache_path <<- path
    desc
  }, error = function(e) {
    message("Could not load description file: ", e$message)
    NULL
  })
}

# Left-join DE table with description file.
# Adds annotation columns (GO terms, Symbol, KEGG, etc.) that are NOT already present.
merge_with_description <- function(de_df, desc = NULL, replace_gene_id = FALSE) {
  if (is.null(desc)) desc <- load_description_file()
  if (is.null(desc)) return(de_df)
  desc <- normalize_annotation_table(desc)
  de_df$.original_gene_id_for_annotation <- as.character(de_df$gene_id)
  de_df$annotation_key <- gene_join_key(de_df$.original_gene_id_for_annotation)
  if (isTRUE(replace_gene_id)) {
    de_df$.annotation_lookup_row <- seq_len(nrow(de_df))
    lookup <- make_annotation_lookup(desc)
    de_df <- merge(
      de_df,
      lookup[, c("lookup_key", "annotation_key", "annotation_gene_id", "annotation_id_source_col"), drop = FALSE],
      by.x = "annotation_key",
      by.y = "lookup_key",
      all.x = TRUE,
      sort = FALSE,
      suffixes = c("", ".matched")
    )
    de_df <- de_df[order(de_df$.annotation_lookup_row), , drop = FALSE]
    de_df$.annotation_lookup_row <- NULL
    de_df$annotation_key <- ifelse(
      !is.na(de_df$annotation_key.matched) & de_df$annotation_key.matched != "",
      de_df$annotation_key.matched,
      de_df$annotation_key
    )
    de_df$annotation_key.matched <- NULL
  }
  new_cols <- setdiff(names(desc), c(names(de_df), "gene_id", "annotation_key"))
  if (length(new_cols) == 0) {
    if (isTRUE(replace_gene_id) && "annotation_gene_id" %in% names(de_df)) {
      if (!"original_gene_id" %in% names(de_df)) de_df$original_gene_id <- de_df$.original_gene_id_for_annotation
      replace_hit <- !is.na(de_df$annotation_gene_id) & de_df$annotation_gene_id != ""
      de_df$gene_id[replace_hit] <- de_df$annotation_gene_id[replace_hit]
      de_df$annotation_gene_id <- NULL
    }
    de_df$.original_gene_id_for_annotation <- NULL
    de_df$annotation_key <- NULL
    return(de_df)
  }
  desc_sub <- desc[, c("annotation_key", new_cols), drop = FALSE]
  de_df$annotation_row <- seq_len(nrow(de_df))
  merged <- merge(de_df, desc_sub, by = "annotation_key", all.x = TRUE, sort = FALSE)
  merged <- merged[order(merged$annotation_row), , drop = FALSE]
  if (isTRUE(replace_gene_id) && "annotation_gene_id" %in% names(merged)) {
    if (!"original_gene_id" %in% names(merged)) merged$original_gene_id <- merged$.original_gene_id_for_annotation
    replace_hit <- !is.na(merged$annotation_gene_id) & merged$annotation_gene_id != ""
    merged$gene_id[replace_hit] <- merged$annotation_gene_id[replace_hit]
    merged$annotation_gene_id <- NULL
    merged$annotation_id_source_col <- NULL
  }
  merged$annotation_key <- NULL
  merged$annotation_row <- NULL
  merged$.original_gene_id_for_annotation <- NULL
  rownames(merged) <- NULL
  merged
}

standardize_de_table <- function(df, description_df = NULL, merge_default_description = TRUE) {
  df <- as.data.frame(df, check.names = FALSE)

  gene_col <- first_existing_col(df, c("gene_id", "GeneID", "TAIR", "locus_tag"))
  lfc_col  <- first_existing_col(df, c("log2FoldChange", "log2FC", "logFC", "LFC", "fold_change_log2"))
  padj_col <- first_existing_col(df, c("padj", "FDR", "adj.P.Val", "p.adjust", "p_adj"))
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

  df$gene_id <- clean_gene_id(df$gene_id)
  df$log2FoldChange <- suppressWarnings(as.numeric(df$log2FoldChange))
  df$padj <- as_numeric_p(df$padj)
  if (!"pValue" %in% names(df)) df$pValue <- df$padj
  df$pValue <- as_numeric_p(df$pValue)
  if ("baseMean" %in% names(df)) df$baseMean <- suppressWarnings(as.numeric(df$baseMean))
  df <- df[!is.na(df$gene_id) & df$gene_id != "", , drop = FALSE]
  # Merge with selected description file to add GO terms, Symbol, KEGG, etc.
  if (!is.null(description_df)) {
    df <- merge_with_description(df, description_df)
  } else if (isTRUE(merge_default_description)) {
    df <- merge_with_description(df)
  }
  df
}

classify_de <- function(df, alpha = 0.05, lfc_cutoff = 1) {
  if (nrow(df) == 0) {
    df$DE_class <- character(0)
    return(df)
  }
  df$DE_class <- "not_significant"
  df$DE_class[!is.na(df$padj) & df$padj < alpha & df$log2FoldChange >= lfc_cutoff] <- "up"
  df$DE_class[!is.na(df$padj) & df$padj < alpha & df$log2FoldChange <= -lfc_cutoff] <- "down"
  df
}

plot_theme_choice <- function(plot_theme = "classic", base_size = 12, font_family = "serif") {
  if (is.null(plot_theme) || !nzchar(plot_theme)) plot_theme <- "classic"
  if (is.null(font_family) || !nzchar(font_family)) font_family <- "serif"
  switch(tolower(plot_theme),
    bw = ggplot2::theme_bw(base_size = base_size, base_family = font_family),
    minimal = ggplot2::theme_minimal(base_size = base_size, base_family = font_family),
    linedraw = ggplot2::theme_linedraw(base_size = base_size, base_family = font_family),
    light = ggplot2::theme_light(base_size = base_size, base_family = font_family),
    gray = ggplot2::theme_gray(base_size = base_size, base_family = font_family),
    dark = ggplot2::theme_dark(base_size = base_size, base_family = font_family),
    void = ggplot2::theme_void(base_size = base_size, base_family = font_family),
    classic = ggplot2::theme_classic(base_size = base_size, base_family = font_family),
    ggplot2::theme_classic(base_size = base_size, base_family = font_family)
  )
}

make_volcano_plot <- function(df, alpha = 0.05, lfc_cutoff = 1, title = "Volcano plot",
                              point_size = 1, point_alpha = 0.65,
                              color_up = "#B2182B", color_down = "#2166AC", color_ns = "grey70",
                              plot_theme = "classic", font_family = "serif") {
  df <- classify_de(df, alpha, lfc_cutoff)
  df$neg_log10_padj <- -log10(pmax(df$padj, .Machine$double.xmin))
  ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange, y = neg_log10_padj, color = DE_class)) +
    ggplot2::geom_point(alpha = point_alpha, size = point_size) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed", linewidth = 0.3) +
    ggplot2::scale_color_manual(values = c(up = color_up, down = color_down, not_significant = color_ns)) +
    plot_theme_choice(plot_theme, base_size = 12, font_family = font_family) +
    ggplot2::labs(x = "log2 fold change", y = "-log10 adjusted p-value", color = "Class", title = title)
}

make_ma_plot <- function(df, alpha = 0.05, lfc_cutoff = 1, title = "MA plot",
                         point_size = 1, point_alpha = 0.65,
                         color_up = "#B2182B", color_down = "#2166AC", color_ns = "grey70",
                         plot_theme = "classic", font_family = "serif") {
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
    plot_theme_choice(plot_theme, base_size = 12, font_family = font_family) +
    ggplot2::labs(x = "Mean of normalized counts", y = "log2 fold change", color = "Class", title = title)
}

make_gene_norm_counts_boxplot <- function(norm_counts, coldata, gene_id,
                                          treatment = NULL, control = NULL,
                                          plot_theme = "classic", font_family = "serif",
                                          color_trnt = "#ac783e", color_ctrl = "#505050",
                                          point_size = 2.4, point_alpha = 0.88,
                                          jitter_width = 0.12, box_width = 0.55,
                                          box_alpha = 0.28) {
  if (is.null(norm_counts) || is.null(coldata)) stop("Normalized counts and colData are required.")
  if (!"gene_id" %in% names(norm_counts)) stop("Normalized counts table is missing gene_id.")
  if (!all(c("sample_id", "condition") %in% names(coldata))) stop("colData must include sample_id and condition.")

  gene_key <- clean_gene_id(gene_id)
  gene_match <- clean_gene_id(norm_counts$gene_id) == gene_key
  if (!any(gene_match, na.rm = TRUE)) {
    gene_match <- gene_join_key(norm_counts$gene_id) == gene_join_key(gene_id)
  }
  if (!any(gene_match, na.rm = TRUE)) stop("Gene was not found in the normalized counts table: ", gene_id)

  cd <- as.data.frame(coldata, check.names = FALSE)
  groups <- c(control, treatment)
  groups <- groups[!is.na(groups) & nzchar(groups)]
  groups <- unique(groups)
  if (length(groups) > 0) {
    cd <- cd[as.character(cd$condition) %in% groups, , drop = FALSE]
    cd$condition <- factor(as.character(cd$condition), levels = groups)
  }

  sample_cols <- as.character(cd$sample_id)
  sample_cols <- sample_cols[sample_cols %in% names(norm_counts)]
  if (length(sample_cols) == 0) stop("No selected samples were found in the normalized counts table.")

  counts <- as.numeric(norm_counts[which(gene_match)[1], sample_cols, drop = TRUE])
  plot_df <- data.frame(
    sample_id = sample_cols,
    normalized_count = counts,
    stringsAsFactors = FALSE
  )
  plot_df <- merge(plot_df, cd, by = "sample_id", all.x = TRUE, sort = FALSE)
  plot_df <- plot_df[!is.na(plot_df$condition) & !is.na(plot_df$normalized_count), , drop = FALSE]
  if (nrow(plot_df) == 0) stop("No plottable normalized counts for this gene and comparison.")

  plot_df$plot_count <- plot_df$normalized_count
  condition_levels <- levels(factor(plot_df$condition))
  cols <- if (length(condition_levels) == 2) c(color_ctrl, color_trnt) else pca_palette_values(length(condition_levels), "set2")
  names(cols) <- condition_levels

  ggplot2::ggplot(plot_df, ggplot2::aes(x = condition, y = plot_count)) +
    ggplot2::geom_boxplot(width = box_width, alpha = box_alpha, outlier.shape = NA) +
    ggplot2::geom_jitter(width = jitter_width, height = 0, size = point_size, alpha = point_alpha, ggplot2::aes(color = condition)) +
    ggplot2::scale_color_manual(values = cols, guide = "none") +
    plot_theme_choice(plot_theme, base_size = 12, font_family = font_family) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 14, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 13),
    ) +
    ggplot2::labs(
      title = gene_id,
      x = NULL,
      y = "Normalized count\n"
    )
}

pca_palette_values <- function(n, palette = "default") {
  if (is.null(palette) || !nzchar(palette) || palette == "default") return(NULL)
  base_cols <- switch(tolower(palette),
    okabe_ito = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"),
    set1 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"),
    set2 = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"),
    dark2 = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"),
    paired = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"),
    tableau = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"),
    viridis = c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725"),
    plasma = c("#0D0887", "#7E03A8", "#CC4778", "#F89540", "#F0F921"),
    pastel = c("#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC"),
    NULL
  )
  if (is.null(base_cols)) return(NULL)
  if (n <= length(base_cols)) return(base_cols[seq_len(n)])
  grDevices::colorRampPalette(base_cols)(n)
}

make_pca_plot <- function(pca_df, title = "PCA", point_size = 3.2, point_alpha = 0.9,
                          show_labels = TRUE, conditions = NULL,
                          plot_theme = "classic", font_family = "serif",
                          color_palette = "default") {
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
  condition_levels <- unique(as.character(pca_df$condition))
  palette_values <- pca_palette_values(length(condition_levels), color_palette)
  if (!is.null(palette_values)) {
    names(palette_values) <- condition_levels
    p <- p + ggplot2::scale_color_manual(values = palette_values)
  }
  if (isTRUE(show_labels)) {
    p <- p + ggrepel::geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = 30)
  }
  p +
    plot_theme_choice(plot_theme, base_size = 12, font_family = font_family) +
    ggplot2::labs(x = pc1lab, y = pc2lab, color = "Condition", title = title)
}

normalize_coldata <- function(coldata, sample_col = NULL, condition_col = NULL, label_col = NULL) {
  coldata <- as.data.frame(coldata, check.names = FALSE)
  if (is.null(sample_col)) {
    sample_candidates <- seq_len(min(2, ncol(coldata)))
    sample_is_unique <- vapply(sample_candidates, function(i) {
      vals <- trimws(as.character(coldata[[i]]))
      vals <- vals[!is.na(vals) & nzchar(vals)]
      length(vals) == nrow(coldata) && !anyDuplicated(vals)
    }, logical(1))
    sample_col <- names(coldata)[sample_candidates[sample_is_unique][1]]
  }
  if (is.null(condition_col)) {
    condition_candidates <- seq_len(min(3, ncol(coldata)))
    condition_is_group <- vapply(condition_candidates, function(i) {
      vals <- trimws(as.character(coldata[[i]]))
      vals <- vals[!is.na(vals) & nzchar(vals)]
      length(vals) == nrow(coldata) && anyDuplicated(vals) > 0
    }, logical(1))
    condition_col <- names(coldata)[condition_candidates[condition_is_group][1]]
  }
  if (is.null(label_col)) {
    label_candidates <- if (ncol(coldata) >= 2) 2:min(4, ncol(coldata)) else integer()
    if (length(label_candidates) > 0) {
      label_candidates <- label_candidates[!names(coldata)[label_candidates] %in% c(sample_col, condition_col)]
      label_is_unique <- vapply(label_candidates, function(i) {
        vals <- trimws(as.character(coldata[[i]]))
        vals <- vals[!is.na(vals) & nzchar(vals)]
        length(vals) == nrow(coldata) && !anyDuplicated(vals)
      }, logical(1))
      label_col <- names(coldata)[label_candidates[label_is_unique][1]]
    }
  }
  if (is.null(sample_col) || is.na(sample_col)) {
    stop("Could not find a sample/file column in colData. It must be column 1 or 2 and contain unique non-empty values.")
  }
  if (is.null(condition_col) || is.na(condition_col)) {
    stop("Could not find a condition/group column in colData. It must be column 1, 2, or 3 and contain repeated non-empty group values.")
  }
  out <- data.frame(
    sample_id = as.character(coldata[[sample_col]]),
    condition = as.character(coldata[[condition_col]]),
    stringsAsFactors = FALSE
  )
  if (!is.null(label_col) && !is.na(label_col) && label_col %in% names(coldata)) {
    out$sample_label <- as.character(coldata[[label_col]])
  } else {
    out$sample_label <- out$sample_id
  }
  used_cols <- unique(c(sample_col, condition_col, label_col))
  used_cols <- used_cols[!is.na(used_cols) & used_cols %in% names(coldata)]
  extra_cols <- setdiff(names(coldata), used_cols)
  extra_cols <- setdiff(extra_cols, names(out))
  if (length(extra_cols) > 0) {
    out <- cbind(out, coldata[, extra_cols, drop = FALSE])
  }
  out$sample_id <- sub("\\.(genes|transcripts)\\.results$", "", out$sample_id)
  out$sample_id <- basename(out$sample_id)
  out <- out[!is.na(out$sample_id) & out$sample_id != "", , drop = FALSE]
  out
}

scan_rsem_files <- function(folder, transcript_level = FALSE) {
  if (is.null(folder) || !dir.exists(folder)) return(data.frame())
  suffix <- if (isTRUE(transcript_level)) "transcripts" else "genes"
  pattern <- paste0("\\.", suffix, "\\.results$")
  files <- list.files(folder, pattern = pattern, full.names = TRUE)
  sample_ids <- sub("\\.(genes|transcripts)\\.results$", "", basename(files))
  data.frame(
    sample_id = sample_ids,
    condition = "condition_1",
    sample_label = sample_ids,
    file = files,
    stringsAsFactors = FALSE
  )
}

scan_tximport_quant_files <- function(folder, quant_type = c("rsem", "salmon", "kallisto"),
                                      rsem_tx_ids = FALSE) {
  quant_type <- match.arg(quant_type)
  if (identical(quant_type, "rsem")) return(scan_rsem_files(folder, transcript_level = rsem_tx_ids))
  if (is.null(folder) || !dir.exists(folder)) return(data.frame())

  pattern <- switch(
    quant_type,
    salmon = "quant\\.sf$",
    kallisto = "abundance\\.tsv$"
  )
  files <- list.files(folder, pattern = pattern, full.names = TRUE, recursive = TRUE)
  if (!length(files)) return(data.frame())
  sample_ids <- basename(dirname(files))
  duplicated_ids <- duplicated(sample_ids) | duplicated(sample_ids, fromLast = TRUE)
  sample_ids[duplicated_ids] <- tools::file_path_sans_ext(basename(files[duplicated_ids]))
  data.frame(
    sample_id = sample_ids,
    condition = "condition_1",
    sample_label = sample_ids,
    file = files,
    stringsAsFactors = FALSE
  )
}

read_tx2gene_table <- function(path) {
  if (is.null(path) || !file.exists(path)) stop("A tx2gene table is required for transcript-level quantification input.")
  tx2gene <- read_any_table(path, source_name = path)
  tx2gene <- as.data.frame(tx2gene, check.names = FALSE)
  if (ncol(tx2gene) < 2) stop("tx2gene table must contain at least two columns: transcript ID and gene ID.")
  tx2gene <- tx2gene[, 1:2, drop = FALSE]
  names(tx2gene) <- c("TXNAME", "GENEID")
  tx2gene$TXNAME <- trimws(as.character(tx2gene$TXNAME))
  tx2gene$GENEID <- trimws(as.character(tx2gene$GENEID))
  tx2gene <- tx2gene[!is.na(tx2gene$TXNAME) & nzchar(tx2gene$TXNAME) &
                       !is.na(tx2gene$GENEID) & nzchar(tx2gene$GENEID), , drop = FALSE]
  tx2gene <- tx2gene[!duplicated(tx2gene$TXNAME), , drop = FALSE]
  if (nrow(tx2gene) == 0) stop("tx2gene table has no usable transcript-to-gene mappings.")
  tx2gene
}

clean_featurecounts_sample_names <- function(x) {
  x <- basename(as.character(x))
  x <- sub("\\.sorted\\.bam$", "", x, ignore.case = TRUE)
  x <- sub("\\.bam$", "", x, ignore.case = TRUE)
  x
}

read_featurecounts_count_matrix <- function(path) {
  fc <- utils::read.delim(path, comment.char = "#", check.names = FALSE)
  if (nrow(fc) == 0) stop("featureCounts file has no rows.")
  gene_col <- first_existing_col(fc, c("Geneid", "gene_id", "GeneID", "gene"))
  if (is.null(gene_col)) gene_col <- names(fc)[1]
  count_anchor <- first_existing_col(fc, c("gene_biotype", "Length"))
  if (is.null(count_anchor)) {
    stop("Could not find gene_biotype or Length column in the featureCounts file.")
  }
  count_start <- match(count_anchor, names(fc)) + 1
  if (is.na(count_start) || count_start > ncol(fc)) {
    stop("No sample count columns were found after ", count_anchor, ".")
  }
  count_mat <- fc[, count_start:ncol(fc), drop = FALSE]
  colnames(count_mat) <- clean_featurecounts_sample_names(colnames(count_mat))
  count_mat <- as.matrix(count_mat)
  suppressWarnings(storage.mode(count_mat) <- "numeric")
  if (anyNA(count_mat)) stop("featureCounts sample columns contain non-numeric values.")
  count_mat <- round(count_mat)
  storage.mode(count_mat) <- "integer"
  rownames(count_mat) <- clean_gene_id(fc[[gene_col]])
  keep <- !is.na(rownames(count_mat)) & nzchar(rownames(count_mat))
  count_mat <- count_mat[keep, , drop = FALSE]
  if (nrow(count_mat) == 0 || ncol(count_mat) == 0) stop("No usable count matrix could be extracted from featureCounts.")
  if (anyDuplicated(rownames(count_mat))) {
    count_mat <- rowsum(count_mat, group = rownames(count_mat), reorder = FALSE)
    storage.mode(count_mat) <- "integer"
  }
  count_mat
}

read_count_matrix_file <- function(path) {
  tbl <- read_any_table(path, source_name = basename(path))
  tbl <- as.data.frame(tbl, check.names = FALSE)
  if (ncol(tbl) < 2) stop("Count matrix must contain a gene ID column and at least one sample count column.")

  gene_ids <- clean_gene_id(tbl[[1]])
  count_df <- tbl[, -1, drop = FALSE]
  if (any(!nzchar(names(count_df)))) stop("All sample count columns must have names.")
  count_mat <- as.matrix(count_df)
  suppressWarnings(storage.mode(count_mat) <- "numeric")
  if (anyNA(count_mat)) stop("Count matrix sample columns contain non-numeric values.")
  count_mat <- round(count_mat)
  storage.mode(count_mat) <- "integer"
  rownames(count_mat) <- gene_ids
  keep <- !is.na(rownames(count_mat)) & nzchar(rownames(count_mat))
  count_mat <- count_mat[keep, , drop = FALSE]
  if (nrow(count_mat) == 0 || ncol(count_mat) == 0) stop("No usable count matrix could be extracted.")
  if (anyDuplicated(rownames(count_mat))) {
    count_mat <- rowsum(count_mat, group = rownames(count_mat), reorder = FALSE)
    storage.mode(count_mat) <- "integer"
  }
  count_mat
}

scan_count_matrix_file <- function(path) {
  count_mat <- read_count_matrix_file(path)
  data.frame(
    sample_id = colnames(count_mat),
    condition = "condition_1",
    sample_label = colnames(count_mat),
    stringsAsFactors = FALSE
  )
}

scan_featurecounts_file <- function(path) {
  if (is.null(path) || !file.exists(path)) return(data.frame())
  count_mat <- read_featurecounts_count_matrix(path)
  data.frame(
    sample_id = colnames(count_mat),
    condition = "condition_1",
    sample_label = colnames(count_mat),
    stringsAsFactors = FALSE
  )
}

run_deseq2_from_rsem <- function(folder, coldata, treatment, control, lfc_shrink = FALSE, min_count = 10,
                                 effect_col = NULL, effect_level = NULL, use_interaction = FALSE,
                                 all_vs_control = TRUE, quant_type = "rsem", tx2gene_file = NULL,
                                 rsem_tx_ids = FALSE) {
  quant_type <- match.arg(quant_type, c("rsem", "salmon", "kallisto"))
  required <- c("DESeq2", "tximport", "tibble")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "))

  quant_table <- scan_tximport_quant_files(folder, quant_type, rsem_tx_ids = rsem_tx_ids)
  if (nrow(quant_table) == 0) {
    expected_file <- switch(
      quant_type,
      rsem = if (isTRUE(rsem_tx_ids)) ".transcripts.results" else ".genes.results",
      salmon = "quant.sf",
      kallisto = "abundance.tsv"
    )
    stop("No ", expected_file, " files were found in the selected folder.")
  }

  coldata <- normalize_coldata(coldata, sample_col = "sample_id", condition_col = "condition", label_col = "sample_label")
  if (identical(quant_type, "rsem")) {
    coldata$sample_id <- sub("\\.(genes|transcripts)\\.results$", "", basename(coldata$sample_id))
  }
  rownames(coldata) <- coldata$sample_id
  quant_table <- quant_table[match(coldata$sample_id, quant_table$sample_id), , drop = FALSE]
  if (any(is.na(quant_table$file))) {
    missing_samples <- coldata$sample_id[is.na(quant_table$file)]
    stop("These colData sample_id values do not match ", quant_type, " quantification files: ", paste(missing_samples, collapse = ", "))
  }

  if (!treatment %in% coldata$condition) stop("Treatment is not present in colData condition column.")
  if (!control %in% coldata$condition) stop("Control is not present in colData condition column.")

  coldata$condition <- factor(coldata$condition)
  coldata$condition <- stats::relevel(coldata$condition, ref = control)
  design_formula <- ~ condition
  results_contrast <- c("condition", treatment, control)
  contrast_label <- paste0("condition: ", treatment, " vs ", control)
  design_label <- "~ condition"

  if (!is.null(effect_col) && nzchar(effect_col) && identical(effect_col, "condition")) {
    if (isTRUE(use_interaction)) {
      stop("Interaction requires a separate effect column in colData, not the condition column itself.")
    }
    design_label <- "~ condition"
    contrast_label <- paste0("condition: ", treatment, " vs ", control, " in a multi-level condition column")
  } else if (!is.null(effect_col) && nzchar(effect_col)) {
    if (!effect_col %in% names(coldata)) stop("Selected effect column is not present in colData.")
    effect_values <- trimws(as.character(coldata[[effect_col]]))
    effect_values[!nzchar(effect_values)] <- NA_character_
    coldata$de_effect <- factor(effect_values)
    coldata <- coldata[!is.na(coldata$de_effect), , drop = FALSE]
    if (nrow(coldata) == 0) stop("No samples remain after removing missing values from the selected effect column.")
    if (!treatment %in% coldata$condition) stop("Treatment is not present after removing missing values from the selected effect column.")
    if (!control %in% coldata$condition) stop("Control is not present after removing missing values from the selected effect column.")
    if (is.null(effect_level) || !nzchar(effect_level)) {
      effect_level <- levels(coldata$de_effect)[1]
    }
    if (!effect_level %in% levels(coldata$de_effect)) stop("Selected effect level is not present in the effect column.")
    coldata$de_effect <- stats::relevel(coldata$de_effect, ref = effect_level)

    if (isTRUE(use_interaction)) {
      has_treatment_at_effect <- any(coldata$condition == treatment & coldata$de_effect == effect_level)
      has_control_at_effect <- any(coldata$condition == control & coldata$de_effect == effect_level)
      if (!has_treatment_at_effect || !has_control_at_effect) {
        stop("For the interaction test, both treatment and control must have samples at effect level '", effect_level, "'.")
      }
      design_formula <- ~ condition + de_effect + condition:de_effect
      design_label <- paste0("~ condition + ", effect_col, " + condition:", effect_col)
      contrast_label <- paste0("condition: ", treatment, " vs ", control, " at ", effect_col, " = ", effect_level, " (", effect_col, " re-leveled as reference)")
    } else {
      design_formula <- ~ condition + de_effect
      design_label <- paste0("~ condition + ", effect_col)
      contrast_label <- paste0("condition: ", treatment, " vs ", control, " adjusted for ", effect_col, " (reference: ", effect_level, ")")
    }
  }

  coldata$condition <- droplevels(coldata$condition)
  if ("de_effect" %in% names(coldata)) coldata$de_effect <- droplevels(coldata$de_effect)

  quant_table <- quant_table[match(coldata$sample_id, quant_table$sample_id), , drop = FALSE]
  files <- quant_table$file
  names(files) <- coldata$sample_id

  if (identical(quant_type, "rsem")) {
    if (isTRUE(rsem_tx_ids)) {
      tx2gene <- read_tx2gene_table(tx2gene_file)
      txi <- tximport::tximport(files, type = "rsem", txIn = TRUE, txOut = FALSE, tx2gene = tx2gene)
    } else {
      txi <- tximport::tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
    }
  } else {
    tx2gene <- read_tx2gene_table(tx2gene_file)
    txi <- tximport::tximport(files, type = quant_type, tx2gene = tx2gene)
  }
  txi$length[txi$length == 0] <- 1
  dds <- DESeq2::DESeqDataSetFromTximport(txi, colData = coldata, design = design_formula)
  keep <- rowSums(DESeq2::counts(dds)) >= min_count
  dds <- dds[keep, ]
  dds <- DESeq2::DESeq(dds)

  result_table_for_contrast <- function(contrast_vec) {
    res <- DESeq2::results(dds, contrast = contrast_vec, alpha = 0.05)
    if (isTRUE(lfc_shrink)) {
      if (requireNamespace("ashr", quietly = TRUE)) {
        res <- DESeq2::lfcShrink(dds, contrast = contrast_vec, res = res, type = "ashr")
      } else {
        warning("Package 'ashr' is not installed, so lfcShrink was skipped.")
      }
    }
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    res_df <- tibble::as_tibble(res_df)
    if ("pvalue" %in% names(res_df)) names(res_df)[names(res_df) == "pvalue"] <- "pValue"
    res_df <- res_df[, c("gene_id", setdiff(names(res_df), "gene_id"))]
    res_df <- standardize_de_table(res_df, merge_default_description = FALSE)
    res_df[order(res_df$padj), , drop = FALSE]
  }

  res <- DESeq2::results(dds, contrast = results_contrast, alpha = 0.05)
  if (isTRUE(lfc_shrink)) {
    if (requireNamespace("ashr", quietly = TRUE)) {
      res <- DESeq2::lfcShrink(dds, contrast = results_contrast, res = res, type = "ashr")
    } else {
      warning("Package 'ashr' is not installed, so lfcShrink was skipped.")
    }
  }
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- tibble::as_tibble(res_df)
  if ("pvalue" %in% names(res_df)) names(res_df)[names(res_df) == "pvalue"] <- "pValue"
  res_df <- res_df[, c("gene_id", setdiff(names(res_df), "gene_id"))]
  res_df <- standardize_de_table(res_df, merge_default_description = FALSE)
  res_df <- res_df[order(res_df$padj), , drop = FALSE]

  all_comparisons <- NULL
  if (isTRUE(all_vs_control)) {
    comparison_levels <- setdiff(levels(coldata$condition), control)
    comparison_tables <- lapply(comparison_levels, function(level) {
      result_table_for_contrast(c("condition", level, control))
    })
    names(comparison_tables) <- paste0(comparison_levels, "_vs_", control)
    all_comparisons <- list(
      control = control,
      comparisons = comparison_levels,
      tables = comparison_tables
    )
  }

  norm_counts <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
  norm_counts$gene_id <- clean_gene_id(rownames(norm_counts))
  norm_counts <- norm_counts[, c("gene_id", setdiff(names(norm_counts), "gene_id"))]

  vst_obj <- if (nrow(dds) < 1000) {
    DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
  } else {
    DESeq2::vst(dds, blind = FALSE)
  }
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
    all_comparisons = all_comparisons,
    summary = paste(c(
      paste0("Design formula: ", design_label),
      paste0("Contrast: ", contrast_label),
      if (isTRUE(all_vs_control)) paste0("All comparisons vs control: ", paste(setdiff(levels(coldata$condition), control), collapse = ", "), " vs ", control) else NULL,
      capture.output(summary(res))
    ), collapse = "\n"),
    design_formula = design_label,
    contrast = contrast_label,
    coldata = coldata
  )
}

run_deseq2_from_featurecounts <- function(counts_file, coldata, treatment, control, lfc_shrink = FALSE, min_count = 10,
                                          effect_col = NULL, effect_level = NULL, use_interaction = FALSE,
                                          all_vs_control = TRUE) {
  required <- c("DESeq2", "tibble")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "))

  count_mat <- read_featurecounts_count_matrix(counts_file)
  coldata <- normalize_coldata(coldata, sample_col = "sample_id", condition_col = "condition", label_col = "sample_label")
  coldata$sample_id <- clean_featurecounts_sample_names(coldata$sample_id)
  rownames(coldata) <- coldata$sample_id
  missing_samples <- setdiff(coldata$sample_id, colnames(count_mat))
  if (length(missing_samples) > 0) {
    stop("These colData sample_id values do not match featureCounts columns: ", paste(missing_samples, collapse = ", "))
  }
  count_mat <- count_mat[, coldata$sample_id, drop = FALSE]

  if (!treatment %in% coldata$condition) stop("Treatment is not present in colData condition column.")
  if (!control %in% coldata$condition) stop("Control is not present in colData condition column.")

  coldata$condition <- factor(coldata$condition)
  coldata$condition <- stats::relevel(coldata$condition, ref = control)
  design_formula <- ~ condition
  results_contrast <- c("condition", treatment, control)
  contrast_label <- paste0("condition: ", treatment, " vs ", control)
  design_label <- "~ condition"

  if (!is.null(effect_col) && nzchar(effect_col) && identical(effect_col, "condition")) {
    if (isTRUE(use_interaction)) {
      stop("Interaction requires a separate effect column in colData, not the condition column itself.")
    }
    design_label <- "~ condition"
    contrast_label <- paste0("condition: ", treatment, " vs ", control, " in a multi-level condition column")
  } else if (!is.null(effect_col) && nzchar(effect_col)) {
    if (!effect_col %in% names(coldata)) stop("Selected effect column is not present in colData.")
    effect_values <- trimws(as.character(coldata[[effect_col]]))
    effect_values[!nzchar(effect_values)] <- NA_character_
    coldata$de_effect <- factor(effect_values)
    coldata <- coldata[!is.na(coldata$de_effect), , drop = FALSE]
    if (nrow(coldata) == 0) stop("No samples remain after removing missing values from the selected effect column.")
    if (!treatment %in% coldata$condition) stop("Treatment is not present after removing missing values from the selected effect column.")
    if (!control %in% coldata$condition) stop("Control is not present after removing missing values from the selected effect column.")
    count_mat <- count_mat[, coldata$sample_id, drop = FALSE]
    if (is.null(effect_level) || !nzchar(effect_level)) {
      effect_level <- levels(coldata$de_effect)[1]
    }
    if (!effect_level %in% levels(coldata$de_effect)) stop("Selected effect level is not present in the effect column.")
    coldata$de_effect <- stats::relevel(coldata$de_effect, ref = effect_level)

    if (isTRUE(use_interaction)) {
      has_treatment_at_effect <- any(coldata$condition == treatment & coldata$de_effect == effect_level)
      has_control_at_effect <- any(coldata$condition == control & coldata$de_effect == effect_level)
      if (!has_treatment_at_effect || !has_control_at_effect) {
        stop("For the interaction test, both treatment and control must have samples at effect level '", effect_level, "'.")
      }
      design_formula <- ~ condition + de_effect + condition:de_effect
      design_label <- paste0("~ condition + ", effect_col, " + condition:", effect_col)
      contrast_label <- paste0("condition: ", treatment, " vs ", control, " at ", effect_col, " = ", effect_level, " (", effect_col, " re-leveled as reference)")
    } else {
      design_formula <- ~ condition + de_effect
      design_label <- paste0("~ condition + ", effect_col)
      contrast_label <- paste0("condition: ", treatment, " vs ", control, " adjusted for ", effect_col, " (reference: ", effect_level, ")")
    }
  }

  coldata$condition <- droplevels(coldata$condition)
  if ("de_effect" %in% names(coldata)) coldata$de_effect <- droplevels(coldata$de_effect)

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = design_formula)
  keep <- rowSums(DESeq2::counts(dds)) >= min_count
  dds <- dds[keep, ]
  dds <- DESeq2::DESeq(dds)

  result_table_for_contrast <- function(contrast_vec) {
    res <- DESeq2::results(dds, contrast = contrast_vec, alpha = 0.05)
    if (isTRUE(lfc_shrink)) {
      if (requireNamespace("ashr", quietly = TRUE)) {
        res <- DESeq2::lfcShrink(dds, contrast = contrast_vec, res = res, type = "ashr")
      } else {
        warning("Package 'ashr' is not installed, so lfcShrink was skipped.")
      }
    }
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    res_df <- tibble::as_tibble(res_df)
    if ("pvalue" %in% names(res_df)) names(res_df)[names(res_df) == "pvalue"] <- "pValue"
    res_df <- res_df[, c("gene_id", setdiff(names(res_df), "gene_id"))]
    res_df <- standardize_de_table(res_df, merge_default_description = FALSE)
    res_df[order(res_df$padj), , drop = FALSE]
  }

  res <- DESeq2::results(dds, contrast = results_contrast, alpha = 0.05)
  if (isTRUE(lfc_shrink)) {
    if (requireNamespace("ashr", quietly = TRUE)) {
      res <- DESeq2::lfcShrink(dds, contrast = results_contrast, res = res, type = "ashr")
    } else {
      warning("Package 'ashr' is not installed, so lfcShrink was skipped.")
    }
  }
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- tibble::as_tibble(res_df)
  if ("pvalue" %in% names(res_df)) names(res_df)[names(res_df) == "pvalue"] <- "pValue"
  res_df <- res_df[, c("gene_id", setdiff(names(res_df), "gene_id"))]
  res_df <- standardize_de_table(res_df, merge_default_description = FALSE)
  res_df <- res_df[order(res_df$padj), , drop = FALSE]

  all_comparisons <- NULL
  if (isTRUE(all_vs_control)) {
    comparison_levels <- setdiff(levels(coldata$condition), control)
    comparison_tables <- lapply(comparison_levels, function(level) {
      result_table_for_contrast(c("condition", level, control))
    })
    names(comparison_tables) <- paste0(comparison_levels, "_vs_", control)
    all_comparisons <- list(
      control = control,
      comparisons = comparison_levels,
      tables = comparison_tables
    )
  }

  norm_counts <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
  norm_counts$gene_id <- clean_gene_id(rownames(norm_counts))
  norm_counts <- norm_counts[, c("gene_id", setdiff(names(norm_counts), "gene_id"))]

  vst_obj <- if (nrow(dds) < 1000) {
    DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
  } else {
    DESeq2::vst(dds, blind = FALSE)
  }
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
    all_comparisons = all_comparisons,
    summary = paste(c(
      paste0("Input: featureCounts"),
      paste0("Design formula: ", design_label),
      paste0("Contrast: ", contrast_label),
      if (isTRUE(all_vs_control)) paste0("All comparisons vs control: ", paste(setdiff(levels(coldata$condition), control), collapse = ", "), " vs ", control) else NULL,
      capture.output(summary(res))
    ), collapse = "\n"),
    design_formula = design_label,
    contrast = contrast_label,
    coldata = coldata
  )
}

run_deseq2_from_count_matrix <- function(counts_file, coldata, treatment, control, lfc_shrink = FALSE, min_count = 10,
                                         effect_col = NULL, effect_level = NULL, use_interaction = FALSE,
                                         all_vs_control = TRUE) {
  required <- c("DESeq2", "tibble")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "))

  count_mat <- read_count_matrix_file(counts_file)
  coldata <- normalize_coldata(coldata, sample_col = "sample_id", condition_col = "condition", label_col = "sample_label")
  rownames(coldata) <- coldata$sample_id
  missing_samples <- setdiff(coldata$sample_id, colnames(count_mat))
  if (length(missing_samples) > 0) {
    stop("These colData sample_id values do not match count matrix columns: ", paste(missing_samples, collapse = ", "))
  }
  count_mat <- count_mat[, coldata$sample_id, drop = FALSE]

  if (!treatment %in% coldata$condition) stop("Treatment is not present in colData condition column.")
  if (!control %in% coldata$condition) stop("Control is not present in colData condition column.")

  coldata$condition <- factor(coldata$condition)
  coldata$condition <- stats::relevel(coldata$condition, ref = control)
  design_formula <- ~ condition
  results_contrast <- c("condition", treatment, control)
  contrast_label <- paste0("condition: ", treatment, " vs ", control)
  design_label <- "~ condition"

  if (!is.null(effect_col) && nzchar(effect_col) && identical(effect_col, "condition")) {
    if (isTRUE(use_interaction)) {
      stop("Interaction requires a separate effect column in colData, not the condition column itself.")
    }
    design_label <- "~ condition"
    contrast_label <- paste0("condition: ", treatment, " vs ", control, " in a multi-level condition column")
  } else if (!is.null(effect_col) && nzchar(effect_col)) {
    if (!effect_col %in% names(coldata)) stop("Selected effect column is not present in colData.")
    effect_values <- trimws(as.character(coldata[[effect_col]]))
    effect_values[!nzchar(effect_values)] <- NA_character_
    coldata$de_effect <- factor(effect_values)
    coldata <- coldata[!is.na(coldata$de_effect), , drop = FALSE]
    if (nrow(coldata) == 0) stop("No samples remain after removing missing values from the selected effect column.")
    if (!treatment %in% coldata$condition) stop("Treatment is not present after removing missing values from the selected effect column.")
    if (!control %in% coldata$condition) stop("Control is not present after removing missing values from the selected effect column.")
    count_mat <- count_mat[, coldata$sample_id, drop = FALSE]
    if (is.null(effect_level) || !nzchar(effect_level)) {
      effect_level <- levels(coldata$de_effect)[1]
    }
    if (!effect_level %in% levels(coldata$de_effect)) stop("Selected effect level is not present in the effect column.")
    coldata$de_effect <- stats::relevel(coldata$de_effect, ref = effect_level)

    if (isTRUE(use_interaction)) {
      has_treatment_at_effect <- any(coldata$condition == treatment & coldata$de_effect == effect_level)
      has_control_at_effect <- any(coldata$condition == control & coldata$de_effect == effect_level)
      if (!has_treatment_at_effect || !has_control_at_effect) {
        stop("For the interaction test, both treatment and control must have samples at effect level '", effect_level, "'.")
      }
      design_formula <- ~ condition + de_effect + condition:de_effect
      design_label <- paste0("~ condition + ", effect_col, " + condition:", effect_col)
      contrast_label <- paste0("condition: ", treatment, " vs ", control, " at ", effect_col, " = ", effect_level, " (", effect_col, " re-leveled as reference)")
    } else {
      design_formula <- ~ condition + de_effect
      design_label <- paste0("~ condition + ", effect_col)
      contrast_label <- paste0("condition: ", treatment, " vs ", control, " adjusted for ", effect_col, " (reference: ", effect_level, ")")
    }
  }

  coldata$condition <- droplevels(coldata$condition)
  if ("de_effect" %in% names(coldata)) coldata$de_effect <- droplevels(coldata$de_effect)

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = design_formula)
  keep <- rowSums(DESeq2::counts(dds)) >= min_count
  dds <- dds[keep, ]
  dds <- DESeq2::DESeq(dds)

  result_table_for_contrast <- function(contrast_vec) {
    res <- DESeq2::results(dds, contrast = contrast_vec, alpha = 0.05)
    if (isTRUE(lfc_shrink)) {
      if (requireNamespace("ashr", quietly = TRUE)) {
        res <- DESeq2::lfcShrink(dds, contrast = contrast_vec, res = res, type = "ashr")
      } else {
        warning("Package 'ashr' is not installed, so lfcShrink was skipped.")
      }
    }
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    res_df <- tibble::as_tibble(res_df)
    if ("pvalue" %in% names(res_df)) names(res_df)[names(res_df) == "pvalue"] <- "pValue"
    res_df <- res_df[, c("gene_id", setdiff(names(res_df), "gene_id"))]
    res_df <- standardize_de_table(res_df, merge_default_description = FALSE)
    res_df[order(res_df$padj), , drop = FALSE]
  }

  res <- DESeq2::results(dds, contrast = results_contrast, alpha = 0.05)
  if (isTRUE(lfc_shrink)) {
    if (requireNamespace("ashr", quietly = TRUE)) {
      res <- DESeq2::lfcShrink(dds, contrast = results_contrast, res = res, type = "ashr")
    } else {
      warning("Package 'ashr' is not installed, so lfcShrink was skipped.")
    }
  }
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- tibble::as_tibble(res_df)
  if ("pvalue" %in% names(res_df)) names(res_df)[names(res_df) == "pvalue"] <- "pValue"
  res_df <- res_df[, c("gene_id", setdiff(names(res_df), "gene_id"))]
  res_df <- standardize_de_table(res_df, merge_default_description = FALSE)
  res_df <- res_df[order(res_df$padj), , drop = FALSE]

  all_comparisons <- NULL
  if (isTRUE(all_vs_control)) {
    comparison_levels <- setdiff(levels(coldata$condition), control)
    comparison_tables <- lapply(comparison_levels, function(level) {
      result_table_for_contrast(c("condition", level, control))
    })
    names(comparison_tables) <- paste0(comparison_levels, "_vs_", control)
    all_comparisons <- list(
      control = control,
      comparisons = comparison_levels,
      tables = comparison_tables
    )
  }

  norm_counts <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
  norm_counts$gene_id <- clean_gene_id(rownames(norm_counts))
  norm_counts <- norm_counts[, c("gene_id", setdiff(names(norm_counts), "gene_id"))]

  vst_obj <- if (nrow(dds) < 1000) {
    DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
  } else {
    DESeq2::vst(dds, blind = FALSE)
  }
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
    all_comparisons = all_comparisons,
    summary = paste(c(
      paste0("Input: count matrix"),
      paste0("Design formula: ", design_label),
      paste0("Contrast: ", contrast_label),
      if (isTRUE(all_vs_control)) paste0("All comparisons vs control: ", paste(setdiff(levels(coldata$condition), control), collapse = ", "), " vs ", control) else NULL,
      capture.output(summary(res))
    ), collapse = "\n"),
    design_formula = design_label,
    contrast = contrast_label,
    coldata = coldata
  )
}

deseq_comparison_significant_sets <- function(all_comparisons, alpha = 0.05, lfc_cutoff = 1,
                                             direction_filters = NULL) {
  if (is.null(all_comparisons) || is.null(all_comparisons$tables) || length(all_comparisons$tables) == 0) {
    return(list())
  }
  Map(function(tbl, comparison_name) {
    direction <- direction_filters[[comparison_name]] %||% "all"
    if (!direction %in% c("all", "up", "down")) direction <- "all"
    tbl <- classify_de(tbl, alpha = alpha, lfc_cutoff = lfc_cutoff)
    classes <- if (identical(direction, "all")) c("up", "down") else direction
    sort(unique(tbl$gene_id[tbl$DE_class %in% classes]))
  }, all_comparisons$tables, names(all_comparisons$tables))
}

make_comparison_venn_plot <- function(sig_sets, max_sets = 5, title = "Shared significant genes",
                                      fill_palette = "default", line_color = "#404040",
                                      selected_comparison = NULL, font_family = "serif") {
  sig_sets <- sig_sets[!vapply(sig_sets, is.null, logical(1))]
  if (length(sig_sets) == 0) {
    return(
      ggplot2::ggplot() +
        plot_theme_choice("void", base_size = 12, font_family = font_family) +
        ggplot2::annotate("text", x = 0, y = 0, label = "No all-vs-control comparisons available.", size = 4, family = font_family)
    )
  }
  sig_sets <- lapply(sig_sets, function(x) unique(as.character(x[!is.na(x) & nzchar(x)])))
  total_sets <- length(sig_sets)
  if (!is.null(selected_comparison) && nzchar(selected_comparison) && selected_comparison %in% names(sig_sets)) {
    other_names <- setdiff(names(sig_sets), selected_comparison)
    sig_sets <- sig_sets[c(selected_comparison, other_names)]
  }
  if (total_sets > max_sets) sig_sets <- sig_sets[seq_len(max_sets)]
  n <- length(sig_sets)
  raw_set_names <- names(sig_sets)
  set_names <- raw_set_names
  set_names <- gsub("_vs_", " vs ", set_names, fixed = TRUE)
  set_names <- gsub("_", " ", set_names, fixed = TRUE)
  subtitle <- NULL

  circle_df <- function(cx, cy, r, id, n_points = 240) {
    theta <- seq(0, 2 * pi, length.out = n_points)
    data.frame(x = cx + r * cos(theta), y = cy + r * sin(theta), id = id)
  }

  if (is.null(font_family) || !nzchar(font_family)) font_family <- "serif"
  if (is.null(fill_palette) || !nzchar(fill_palette) || fill_palette == "default") {
    cols <- rep("#FFFFFF", n)
    selected_idx <- match(selected_comparison, raw_set_names)
    if (!is.na(selected_idx)) cols[selected_idx] <- "#F28E2B"
  } else {
    cols <- pca_palette_values(n, fill_palette)
    if (is.null(cols)) cols <- c("#B2182B", "#2166AC", "#4D9221")[seq_len(n)]
  }
  line_color <- if (is.null(line_color) || !nzchar(line_color)) "#404040" else line_color
  if (n > 3) {
    if (!requireNamespace("VennDiagram", quietly = TRUE)) {
      return(
        ggplot2::ggplot() +
          plot_theme_choice("void", base_size = 12, font_family = font_family) +
          ggplot2::annotate("text", x = 0, y = 0, label = "Install the VennDiagram package to show 4-5 comparisons.", size = 4, family = font_family)
      )
    }
    if (requireNamespace("futile.logger", quietly = TRUE)) {
      futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    }
    names(sig_sets) <- set_names
    venn_grob <- VennDiagram::venn.diagram(
      x = sig_sets,
      filename = NULL,
      disable.logging = TRUE,
      main = title,
      main.cex = 1.25,
      main.fontfamily = font_family,
      main.fontface = "bold",
      fill = cols,
      alpha = ifelse(cols == "#FFFFFF", 0.92, 0.58),
      col = rep(line_color, n),
      lwd = rep(1.5, n),
      cex = 0.95,
      fontfamily = font_family,
      fontface = "bold",
      cat.cex = 1.22,
      cat.fontfamily = font_family,
      cat.fontface = "plain",
      cat.col = rep("black", n),
      margin = 0.08
    )
    grid::grid.newpage()
    grid::grid.draw(venn_grob)
    return(invisible(venn_grob))
  }
  if (n == 1) {
    circles <- circle_df(0, 0, 1, set_names[1])
    ggplot2::ggplot(circles, ggplot2::aes(x, y, group = id)) +
      ggplot2::geom_polygon(fill = cols[1], alpha = 0.28, color = line_color, linewidth = 0.9) +
      ggplot2::annotate("text", x = 0, y = 0, label = length(sig_sets[[1]]), size = 5.8, fontface = "bold", family = font_family) +
      ggplot2::annotate("text", x = 0, y = 1.25, label = set_names[1], size = 4.1, family = font_family) +
      ggplot2::coord_fixed(xlim = c(-1.5, 1.5), ylim = c(-1.4, 1.5), clip = "off") +
      plot_theme_choice("void", base_size = 12, font_family = font_family) +
      ggplot2::labs(title = title, subtitle = subtitle)
  } else if (n == 2) {
    a <- sig_sets[[1]]
    b <- sig_sets[[2]]
    counts <- c(
      a_only = length(setdiff(a, b)),
      both = length(intersect(a, b)),
      b_only = length(setdiff(b, a))
    )
    circles <- rbind(circle_df(-0.55, 0, 1, set_names[1]), circle_df(0.55, 0, 1, set_names[2]))
    ggplot2::ggplot(circles, ggplot2::aes(x, y, group = id, fill = id)) +
      ggplot2::geom_polygon(alpha = 0.28, color = line_color, linewidth = 0.9) +
      ggplot2::scale_fill_manual(values = stats::setNames(cols[1:2], set_names)) +
      ggplot2::annotate("text", x = -0.82, y = 0, label = counts["a_only"], size = 4.9, fontface = "bold", family = font_family) +
      ggplot2::annotate("text", x = 0, y = 0, label = counts["both"], size = 4.9, fontface = "bold", family = font_family) +
      ggplot2::annotate("text", x = 0.82, y = 0, label = counts["b_only"], size = 4.9, fontface = "bold", family = font_family) +
      ggplot2::annotate("text", x = -0.85, y = 1.22, label = set_names[1], size = 4.0, family = font_family) +
      ggplot2::annotate("text", x = 0.85, y = 1.22, label = set_names[2], size = 4.0, family = font_family) +
      ggplot2::coord_fixed(xlim = c(-1.75, 1.75), ylim = c(-1.25, 1.45), clip = "off") +
      plot_theme_choice("void", base_size = 12, font_family = font_family) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(title = title, subtitle = subtitle)
  } else {
    a <- sig_sets[[1]]
    b <- sig_sets[[2]]
    c <- sig_sets[[3]]
    ab <- intersect(a, b)
    ac <- intersect(a, c)
    bc <- intersect(b, c)
    abc <- Reduce(intersect, sig_sets[1:3])
    counts <- c(
      a_only = length(setdiff(a, union(b, c))),
      b_only = length(setdiff(b, union(a, c))),
      c_only = length(setdiff(c, union(a, b))),
      ab_only = length(setdiff(ab, c)),
      ac_only = length(setdiff(ac, b)),
      bc_only = length(setdiff(bc, a)),
      abc = length(abc)
    )
    circles <- rbind(
      circle_df(-0.55, 0.25, 1, set_names[1]),
      circle_df(0.55, 0.25, 1, set_names[2]),
      circle_df(0, -0.58, 1, set_names[3])
    )
    ggplot2::ggplot(circles, ggplot2::aes(x, y, group = id, fill = id)) +
      ggplot2::geom_polygon(alpha = 0.26, color = line_color, linewidth = 0.9) +
      ggplot2::scale_fill_manual(values = stats::setNames(cols[1:3], set_names)) +
      ggplot2::annotate("text", x = -1.05, y = 0.35, label = counts["a_only"], size = 4.3, fontface = "bold", family = font_family) +
      ggplot2::annotate("text", x = 1.05, y = 0.35, label = counts["b_only"], size = 4.3, fontface = "bold", family = font_family) +
      ggplot2::annotate("text", x = 0, y = -1.05, label = counts["c_only"], size = 4.3, fontface = "bold", family = font_family) +
      ggplot2::annotate("text", x = 0, y = 0.55, label = counts["ab_only"], size = 4.3, fontface = "bold", family = font_family) +
      ggplot2::annotate("text", x = -0.45, y = -0.28, label = counts["ac_only"], size = 4.3, fontface = "bold", family = font_family) +
      ggplot2::annotate("text", x = 0.45, y = -0.28, label = counts["bc_only"], size = 4.3, fontface = "bold", family = font_family) +
      ggplot2::annotate("text", x = 0, y = 0.02, label = counts["abc"], size = 4.5, fontface = "bold", family = font_family) +
      ggplot2::annotate("text", x = -1.05, y = 1.35, label = set_names[1], size = 3.8, family = font_family) +
      ggplot2::annotate("text", x = 1.05, y = 1.35, label = set_names[2], size = 3.8, family = font_family) +
      ggplot2::annotate("text", x = 0, y = -1.85, label = set_names[3], size = 3.8, family = font_family) +
      ggplot2::coord_fixed(xlim = c(-1.8, 1.8), ylim = c(-1.95, 1.55), clip = "off") +
      plot_theme_choice("void", base_size = 12, font_family = font_family) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(title = title, subtitle = subtitle)
  }
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

add_go_bp_column <- function(df, orgdb = "org.At.tair.db", keytype = "TAIR") {
  bp_col <- detect_go_bp_col(df)
  if (!is.null(bp_col)) {
    df$GO_BP_terms <- as.character(df[[bp_col]])
    return(df)
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace(orgdb, quietly = TRUE)) {
    df$GO_BP_terms <- NA_character_
    return(df)
  }
  orgdb_obj <- get(orgdb, envir = asNamespace(orgdb))
  df$.go_lookup_key <- normalize_orgdb_gene_keys(df$gene_id, orgdb = orgdb, keytype = keytype)
  genes <- unique(df$.go_lookup_key)
  genes <- genes[!is.na(genes) & nzchar(genes)]
  annot <- AnnotationDbi::select(orgdb_obj, keys = genes, keytype = keytype, columns = c("GO", "ONTOLOGY"))
  annot <- annot[!is.na(annot$GO) & annot$ONTOLOGY == "BP", , drop = FALSE]
  collapsed <- stats::aggregate(stats::as.formula(paste("GO ~", keytype)), data = annot, FUN = function(x) paste(unique(x), collapse = "; "))
  names(collapsed) <- c(".go_lookup_key", "GO_BP_terms")
  df$.annotation_row <- seq_len(nrow(df))
  df <- merge(df, collapsed, by = ".go_lookup_key", all.x = TRUE, sort = FALSE)
  df <- df[order(df$.annotation_row), , drop = FALSE]
  df$.annotation_row <- NULL
  df$.go_lookup_key <- NULL
  rownames(df) <- NULL
  df
}

run_topgo_enrichment <- function(de_df, direction = c("up", "down", "all"), ontology = "BP", alpha = 0.05, lfc_cutoff = 1,
                                 p_cutoff = 0.01, algorithm = "weight01", statistic = "fisher",
                                 orgdb = "org.At.tair.db", topgo_id = "entrez", keytype = NULL) {
  direction <- match.arg(direction)
  required <- c("topGO", "GO.db", "AnnotationDbi", orgdb)
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages(library(topGO))
  suppressPackageStartupMessages(library(GO.db))
  de_df <- classify_de(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  keytype <- toupper(keytype %||% switch(tolower(topgo_id %||% "entrez"),
    entrez = "ENTREZID",
    symbol = "SYMBOL",
    ensembl = "ENSEMBL",
    refseq = "REFSEQ",
    alias = "ALIAS",
    genename = "GENENAME",
    "ENTREZID"
  ))
  orgdb_obj <- getExportedValue(orgdb, orgdb)
  available_keytypes <- AnnotationDbi::keytypes(orgdb_obj)
  if (!keytype %in% available_keytypes) {
    stop("Gene ID type '", keytype, "' is not available in ", orgdb,
         ". Available key types: ", paste(available_keytypes, collapse = ", "))
  }
  available_columns <- AnnotationDbi::columns(orgdb_obj)
  go_col <- if ("GOALL" %in% available_columns) "GOALL" else "GO"
  ontology_col <- if ("ONTOLOGYALL" %in% available_columns) "ONTOLOGYALL" else "ONTOLOGY"
  bg <- unique(normalize_orgdb_gene_keys(de_df$gene_id, orgdb = orgdb, keytype = keytype))
  if (direction == "up") interesting <- unique(de_df$gene_id[de_df$DE_class == "up"])
  if (direction == "down") interesting <- unique(de_df$gene_id[de_df$DE_class == "down"])
  if (direction == "all") interesting <- unique(de_df$gene_id[de_df$DE_class %in% c("up", "down")])
  interesting <- intersect(normalize_orgdb_gene_keys(interesting, orgdb = orgdb, keytype = keytype), bg)
  if (length(interesting) < 2) stop("Too few significant genes for GO enrichment in this direction.")
  annot <- AnnotationDbi::select(orgdb_obj, keys = bg, keytype = keytype, columns = c(go_col, ontology_col))
  annot <- annot[!is.na(annot[[go_col]]) & annot[[ontology_col]] == ontology, , drop = FALSE]
  if (nrow(annot) == 0) {
    stop("No ", ontology, " GO annotations found for the loaded genes using ", orgdb, " key type ", keytype, ".")
  }
  annot <- annot[!duplicated(annot[, c(keytype, go_col)]), c(keytype, go_col), drop = FALSE]
  gene2go <- split(as.character(annot[[go_col]]), as.character(annot[[keytype]]))

  geneList <- factor(as.integer(bg %in% interesting))
  names(geneList) <- bg
  GOdata <- methods::new("topGOdata", ontology = ontology, allGenes = geneList,
                         geneSelectionFun = function(x) x == 1,
                         annot = get("annFUN.gene2GO", envir = asNamespace("topGO")), gene2GO = gene2go)
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

parse_go_ids <- function(go_ids) {
  go_ids <- trimws(unlist(strsplit(as.character(go_ids %||% ""), ",|;|\\s+")))
  go_ids <- toupper(unique(go_ids[!is.na(go_ids) & nzchar(go_ids)]))
  go_ids <- go_ids[grepl("^GO:[0-9]{7}$", go_ids)]
  if (length(go_ids) == 0) stop("Enter at least one GO ID, for example GO:0008150.")
  go_ids
}

make_go_gene_table <- function(de_df, go_ids, ontology = "BP", orgdb = "org.At.tair.db",
                               keytype = "TAIR", alpha = 0.05, lfc_cutoff = 1) {
  go_ids <- parse_go_ids(go_ids)
  required <- c("AnnotationDbi", "GO.db", orgdb)
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "))

  de_df <- classify_de(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  keytype <- toupper(keytype %||% "TAIR")
  orgdb_obj <- getExportedValue(orgdb, orgdb)
  available_keytypes <- AnnotationDbi::keytypes(orgdb_obj)
  if (!keytype %in% available_keytypes) {
    stop("Gene ID type '", keytype, "' is not available in ", orgdb,
         ". Available key types: ", paste(available_keytypes, collapse = ", "))
  }

  available_columns <- AnnotationDbi::columns(orgdb_obj)
  go_col <- if ("GOALL" %in% available_columns) "GOALL" else "GO"
  ontology_col <- if ("ONTOLOGYALL" %in% available_columns) "ONTOLOGYALL" else "ONTOLOGY"
  de_df$.go_lookup_key <- normalize_orgdb_gene_keys(de_df$gene_id, orgdb = orgdb, keytype = keytype)
  lookup_keys <- unique(de_df$.go_lookup_key)
  lookup_keys <- lookup_keys[!is.na(lookup_keys) & nzchar(lookup_keys)]
  if (length(lookup_keys) == 0) stop("No usable gene IDs were found for GO matching.")

  annot <- AnnotationDbi::select(orgdb_obj, keys = lookup_keys, keytype = keytype, columns = c(go_col, ontology_col))
  annot <- annot[!is.na(annot[[go_col]]) & annot[[go_col]] %in% go_ids & annot[[ontology_col]] == ontology, , drop = FALSE]
  if (nrow(annot) == 0) return(data.frame())
  annot <- annot[!duplicated(annot[, c(keytype, go_col)]), c(keytype, go_col), drop = FALSE]
  names(annot) <- c(".go_lookup_key", "GO_ID")

  term_map <- tryCatch({
    terms <- AnnotationDbi::select(GO.db::GO.db, keys = unique(annot$GO_ID), keytype = "GOID", columns = "TERM")
    terms <- terms[!is.na(terms$GOID), , drop = FALSE]
    stats::setNames(as.character(terms$TERM), as.character(terms$GOID))
  }, error = function(e) character())

  out <- merge(de_df, annot, by = ".go_lookup_key", all.x = FALSE, sort = FALSE)
  out$GO_term <- unname(term_map[as.character(out$GO_ID)])
  out$GO_term[is.na(out$GO_term) | !nzchar(out$GO_term)] <- out$GO_ID[is.na(out$GO_term) | !nzchar(out$GO_term)]
  out$.go_lookup_key <- NULL
  out <- out[order(out$GO_ID, out$padj, out$pValue, na.last = TRUE), , drop = FALSE]
  out <- out[!duplicated(out[, intersect(c("gene_id", "GO_ID"), names(out)), drop = FALSE]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

make_go_bubble_plot <- function(go_df, title = "GO enrichment", top_n = 20, direction = NULL, point_alpha = 0.82,
                                color_up = "#B2182B", color_down = "#2166AC", color_all = "#e6ac6a",
                                plot_theme = "classic", font_family = "serif") {
  if (is.null(go_df) || nrow(go_df) == 0) stop("No GO terms to plot.")
  go_df <- go_df[order(go_df$pValue_num), , drop = FALSE]
  go_df <- head(go_df, top_n)
  go_df$Term_short <- stringr::str_trunc(go_df$Term, 75)
  go_df$neg_log10_p <- -log10(pmax(go_df$pValue_num, .Machine$double.xmin))
  high_col <- color_up
  if (!is.null(direction)) {
    dir_low <- tolower(direction)
    if (dir_low == "down") high_col <- color_down
    if (dir_low == "all") high_col <- color_all
  }
  ggplot2::ggplot(go_df, ggplot2::aes(x = FoldEnrichment, y = stats::reorder(Term_short, FoldEnrichment))) +
    ggplot2::geom_point(ggplot2::aes(size = Significant, color = neg_log10_p), alpha = point_alpha) +
    plot_theme_choice(plot_theme, base_size = 12, font_family = font_family) +
    ggplot2::scale_color_gradient(low = "#979797", high = high_col) +
    ggplot2::labs(x = "Fold enrichment", y = NULL, size = "Genes", color = "-log10(p)", title = title)
}

normalize_msigdb_gene_ids <- function(x, keytype = "SYMBOL") {
  x <- trimws(as.character(x))
  x <- sub("\\.0$", "", x)
  if (is.null(keytype)) keytype <- "SYMBOL"
  keytype <- toupper(keytype)
  if (keytype %in% c("SYMBOL", "TAIR", "ENSEMBL")) x <- toupper(x)
  x[!is.na(x) & nzchar(x)]
}

msigdb_gene_column <- function(msig_df, keytype = "SYMBOL") {
  if (is.null(keytype)) keytype <- "SYMBOL"
  keytype <- toupper(keytype)
  candidates <- switch(keytype,
    ENTREZID = c("ncbi_gene", "entrez_gene", "db_ncbi_gene", "db_entrez_gene"),
    ENSEMBL = c("ensembl_gene", "db_ensembl_gene"),
    TAIR = c("gene_symbol", "db_gene_symbol", "human_gene_symbol"),
    c("gene_symbol", "db_gene_symbol", "human_gene_symbol")
  )
  hit <- intersect(candidates, names(msig_df))
  if (length(hit) == 0) {
    stop("Could not find a matching MSigDB gene ID column for Gene ID type '", keytype, "'. Available columns: ", paste(names(msig_df), collapse = ", "))
  }
  hit[1]
}

load_msigdb_hallmark_sets <- function(species = "Homo sapiens", keytype = "SYMBOL") {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Package 'msigdbr' is required for MSigDB/Hallmark analysis. Run install.bat or install the msigdbr package.")
  }
  if (is.null(species)) species <- "Homo sapiens"
  species <- trimws(species)
  if (!nzchar(species)) species <- "Homo sapiens"
  msig <- tryCatch(
    msigdbr::msigdbr(species = species, collection = "H"),
    error = function(e) suppressWarnings(msigdbr::msigdbr(species = species, category = "H"))
  )
  if (is.null(msig) || nrow(msig) == 0) stop("No MSigDB Hallmark sets were found for species: ", species)
  id_col <- msigdb_gene_column(msig, keytype)
  msig$analysis_gene_id <- normalize_msigdb_gene_ids(msig[[id_col]], keytype)
  msig <- msig[!is.na(msig$analysis_gene_id) & msig$analysis_gene_id != "", , drop = FALSE]
  split(msig$analysis_gene_id, msig$gs_name)
}

run_msigdb_hallmark_enrichment <- function(de_df, direction = c("up", "down", "all"),
                                           species = "Homo sapiens", keytype = "SYMBOL",
                                           alpha = 0.05, lfc_cutoff = 1,
                                           min_set_size = 5, max_set_size = 500,
                                           p_adjust_method = "BH") {
  direction <- match.arg(direction)
  de_df <- classify_de(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  universe <- unique(normalize_msigdb_gene_ids(de_df$gene_id, keytype))
  universe <- universe[!is.na(universe) & nzchar(universe)]
  if (length(universe) < 5) stop("Too few background genes after applying Gene ID type '", keytype, "'.")

  if (direction == "up") interesting <- de_df$gene_id[de_df$DE_class == "up"]
  if (direction == "down") interesting <- de_df$gene_id[de_df$DE_class == "down"]
  if (direction == "all") interesting <- de_df$gene_id[de_df$DE_class %in% c("up", "down")]
  interesting <- intersect(unique(normalize_msigdb_gene_ids(interesting, keytype)), universe)
  if (length(interesting) < 2) stop("Too few significant genes for MSigDB/Hallmark enrichment in this direction.")

  sets <- load_msigdb_hallmark_sets(species = species, keytype = keytype)
  rows <- lapply(names(sets), function(gs_name) {
    set_genes <- intersect(unique(sets[[gs_name]]), universe)
    set_size <- length(set_genes)
    if (set_size < min_set_size || set_size > max_set_size) return(NULL)
    a <- length(intersect(set_genes, interesting))
    b <- set_size - a
    c <- length(setdiff(interesting, set_genes))
    d <- max(length(universe) - a - b - c, 0)
    ft <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = "greater")
    data.frame(
      Hallmark = gs_name,
      Term = gsub("_", " ", sub("^HALLMARK_", "", gs_name)),
      Direction = direction,
      Gene_ID_type = toupper(keytype),
      Species = species,
      Significant_in_set = a,
      Set_size = set_size,
      Significant_total = length(interesting),
      Background_total = length(universe),
      FoldEnrichment = (a / max(length(interesting), 1)) / (set_size / max(length(universe), 1)),
      pValue = ft$p.value,
      Genes = paste(sort(intersect(set_genes, interesting)), collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out) || nrow(out) == 0) stop("No Hallmark gene sets overlap the loaded DE table after size filtering.")
  out$pAdjusted <- stats::p.adjust(out$pValue, method = p_adjust_method)
  out <- out[order(out$pAdjusted, out$pValue), , drop = FALSE]
  rownames(out) <- NULL
  out
}

parse_hallmark_ids <- function(hallmark_ids) {
  hallmark_ids <- trimws(unlist(strsplit(as.character(hallmark_ids %||% ""), ",|;")))
  hallmark_ids <- toupper(unique(hallmark_ids[!is.na(hallmark_ids) & nzchar(hallmark_ids)]))
  hallmark_ids <- gsub("[^A-Z0-9]+", "_", hallmark_ids)
  hallmark_ids <- gsub("^_+|_+$", "", hallmark_ids)
  hallmark_ids <- ifelse(grepl("^HALLMARK_", hallmark_ids), hallmark_ids, paste0("HALLMARK_", hallmark_ids))
  if (length(hallmark_ids) == 0) stop("Enter at least one Hallmark code.")
  hallmark_ids
}

make_msigdb_hallmark_gene_table <- function(de_df, hallmark_ids, species = "Homo sapiens",
                                            keytype = "SYMBOL", alpha = 0.05, lfc_cutoff = 1) {
  hallmark_ids <- parse_hallmark_ids(hallmark_ids)
  de_df <- classify_de(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  de_df$.hallmark_lookup_key <- normalize_msigdb_gene_ids(de_df$gene_id, keytype)
  de_df <- de_df[!is.na(de_df$.hallmark_lookup_key) & nzchar(de_df$.hallmark_lookup_key), , drop = FALSE]
  if (nrow(de_df) == 0) stop("No usable gene IDs were found for Hallmark matching.")

  sets <- load_msigdb_hallmark_sets(species = species, keytype = keytype)
  known_codes <- names(sets)
  rows <- lapply(hallmark_ids, function(code) {
    hit <- known_codes[toupper(known_codes) == toupper(code)]
    if (length(hit) == 0) return(NULL)
    hallmark_id <- hit[1]
    set_genes <- unique(sets[[hallmark_id]])
    sub <- de_df[de_df$.hallmark_lookup_key %in% set_genes, , drop = FALSE]
    if (nrow(sub) == 0) return(NULL)
    sub$Hallmark <- hallmark_id
    sub$Hallmark_term <- gsub("_", " ", sub("^HALLMARK_", "", hallmark_id))
    sub$Hallmark_set_size <- length(set_genes)
    sub
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out) || nrow(out) == 0) return(data.frame())
  out$.hallmark_lookup_key <- NULL
  out <- out[order(out$Hallmark, out$padj, out$pValue, na.last = TRUE), , drop = FALSE]
  out <- out[!duplicated(out[, intersect(c("gene_id", "Hallmark"), names(out)), drop = FALSE]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

make_msigdb_hallmark_plot <- function(msig_df, title = "MSigDB Hallmark enrichment", top_n = 20, point_alpha = 0.82,
                                      color_up = "#B2182B", color_down = "#2166AC", color_all = "#e6ac6a",
                                      plot_theme = "classic", font_family = "serif") {
  if (is.null(msig_df) || nrow(msig_df) == 0) stop("No Hallmark terms to plot.")
  msig_df <- msig_df[order(msig_df$pAdjusted, msig_df$pValue), , drop = FALSE]
  msig_df <- head(msig_df, top_n)
  msig_df$Term_short <- stringr::str_trunc(msig_df$Term, 70)
  msig_df$neg_log10_padj <- -log10(pmax(msig_df$pAdjusted, .Machine$double.xmin))
  direction <- unique(tolower(msig_df$Direction))
  high_col <- color_up
  if (length(direction) == 1 && direction == "down") high_col <- color_down
  if (length(direction) == 1 && direction == "all") high_col <- color_all
  ggplot2::ggplot(msig_df, ggplot2::aes(x = FoldEnrichment, y = stats::reorder(Term_short, FoldEnrichment))) +
    ggplot2::geom_point(ggplot2::aes(size = Significant_in_set, color = neg_log10_padj), alpha = point_alpha) +
    plot_theme_choice(plot_theme, base_size = 12, font_family = font_family) +
    ggplot2::scale_color_gradient(low = "#979797", high = high_col) +
    ggplot2::labs(x = "Fold enrichment", y = NULL, size = "Genes", color = "-log10(FDR)", title = title)
}

pmn_builtin_catalog <- function() {
  data.frame(
    tax_id = c(
      3702, 4577, 3847, 39947, 4530, 4081, 4113, 3694, 3880, 29760,
      15368, 4558, 4555, 4556, 3659, 3661, 3983, 29780, 4615, 3711,
      3712, 3708, 4097, 42345, 3641, 3635, 29730, 29729, 4553, 4565,
      3760, 49390, 214687, 29760, 4087, 29760
    ),
    species = c(
      "Arabidopsis thaliana", "Zea mays", "Glycine max", "Oryza sativa japonica",
      "Oryza sativa", "Solanum lycopersicum", "Solanum tuberosum",
      "Populus trichocarpa", "Medicago truncatula", "Vitis vinifera",
      "Brachypodium distachyon", "Sorghum bicolor", "Setaria italica",
      "Setaria viridis", "Cucumis sativus", "Cucumis melo", "Ricinus communis",
      "Manihot esculenta", "Musa acuminata", "Brassica rapa", "Brassica oleracea",
      "Brassica napus", "Nicotiana tabacum", "Helianthus annuus",
      "Theobroma cacao", "Gossypium raimondii", "Fragaria vesca",
      "Fragaria x ananassa", "Secale cereale", "Triticum aestivum",
      "Prunus persica", "Malus domestica", "Citrus sinensis", "Vitis vinifera",
      "Capsicum annuum", "Vitis vinifera"
    ),
    cyc_db = c(
      "AraCyc", "CornCyc", "SoyCyc", "OryzaCyc", "OryzaCyc", "TomatoCyc",
      "PotatoCyc", "PoplarCyc", "MtruncatulaCyc", "GrapeCyc",
      "BrachypodiumCyc", "SorghumbicolorCyc", "SetariaCyc", "SviridisCyc",
      "CucumberCyc", "MuskmelonCyc", "CastorbeanCyc", "CassavaCyc",
      "BananaCyc", "Brapa_fpscCyc", "Boleracea_oleraceaCyc", "OilseedrapeCyc",
      "Ntabacum_tn90Cyc", "SunflowerCyc", "CocoaCyc", "GraimondiiCyc",
      "Fvesca_vescaCyc", "StrawberryCyc", "RyeCyc", "BreadwheatCyc",
      "PpersicaCyc", "MdomesticaCyc", "SweetorangeCyc", "GrapeCyc",
      "CannuumCyc", "GrapeCyc"
    ),
    stringsAsFactors = FALSE
  )
}

pmn_catalog_choices <- function() {
  catalog <- pmn_builtin_catalog()
  catalog <- catalog[!duplicated(catalog$cyc_db), , drop = FALSE]
  labels <- paste0(catalog$cyc_db, " - ", catalog$species)
  stats::setNames(catalog$cyc_db, labels)
}

pmn_database_for_tax <- function(tax_id, organism_label = NULL) {
  catalog <- pmn_builtin_catalog()
  tax_id <- suppressWarnings(as.integer(tax_id))
  if (!is.na(tax_id)) {
    hit <- catalog[catalog$tax_id == tax_id, , drop = FALSE]
    if (nrow(hit) > 0) return(hit$cyc_db[1])
  }
  organism_label <- trimws(as.character(organism_label %||% ""))
  if (nzchar(organism_label)) {
    label_base <- sub("\\s*\\([^)]*\\)\\s*$", "", organism_label)
    hit <- catalog[tolower(catalog$species) == tolower(label_base), , drop = FALSE]
    if (nrow(hit) > 0) return(hit$cyc_db[1])
    label_low <- tolower(label_base)
    species_low <- tolower(catalog$species)
    partial_match <- grepl(label_low, species_low, fixed = TRUE) |
      vapply(species_low, function(sp) grepl(sp, label_low, fixed = TRUE), logical(1))
    hit <- catalog[partial_match, , drop = FALSE]
    if (nrow(hit) > 0) return(hit$cyc_db[1])
  }
  ""
}

normalize_pmn_gene_ids <- function(x) {
  x <- gene_join_key(x)
  x[!is.na(x) & nzchar(x)]
}

pmn_cyc_file_stem <- function(cyc_db) {
  cyc_db <- trimws(as.character(cyc_db %||% ""))
  stem <- tolower(gsub("[^A-Za-z0-9_]+", "", cyc_db))
  if (!grepl("cyc$", stem)) stem <- paste0(stem, "cyc")
  stem
}

pmn_pathway_file_url <- function(cyc_db) {
  paste0(PMN_PATHWAYS_BASE_URL, "/", pmn_cyc_file_stem(cyc_db), "_pathways.", PMN_PATHWAYS_DATE)
}

load_pmn_pathways_from_tsv <- function(cyc_db) {
  url <- pmn_pathway_file_url(cyc_db)
  tbl <- suppressWarnings(read_any_table(url, delim = "\t"))
  tbl <- as.data.frame(tbl, check.names = FALSE)
  pathway_col <- first_existing_col(tbl, c("Pathway-id", "pathway.id", "PathwayID", "pathway_code"))
  name_col <- first_existing_col(tbl, c("Pathway-name", "pathway.name", "PathwayName", "pathway_name"))
  gene_col <- first_existing_col(tbl, c("Gene-id", "gene_id", "GeneID", "Gene.ID", "Gene"))
  if (is.null(pathway_col) || is.null(name_col) || is.null(gene_col)) {
    stop("PMN pathway file for ", cyc_db, " is missing required columns Pathway-id, Pathway-name, and Gene-id.")
  }
  tbl <- tbl[!is.na(tbl[[pathway_col]]) & nzchar(tbl[[pathway_col]]) &
               !is.na(tbl[[gene_col]]) & nzchar(tbl[[gene_col]]), , drop = FALSE]
  if (nrow(tbl) == 0) stop("PMN pathway file for ", cyc_db, " contains no pathway-gene rows.")
  tbl$pmn_gene_id <- normalize_pmn_gene_ids(tbl[[gene_col]])
  tbl <- tbl[!is.na(tbl$pmn_gene_id) & nzchar(tbl$pmn_gene_id), , drop = FALSE]
  if (nrow(tbl) == 0) stop("No usable gene IDs were found in the PMN pathway file for ", cyc_db, ".")
  genes_by_pathway <- split(tbl$pmn_gene_id, tbl[[pathway_col]])
  genes_by_pathway <- lapply(genes_by_pathway, function(x) sort(unique(x[!is.na(x) & nzchar(x)])))
  genes_by_pathway <- genes_by_pathway[lengths(genes_by_pathway) > 0]
  pathway_names <- stats::aggregate(tbl[[name_col]], by = list(pathway_id = tbl[[pathway_col]]), FUN = function(x) {
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0) "" else x[1]
  })
  pathway_names <- stats::setNames(as.character(pathway_names$x), as.character(pathway_names$pathway_id))
  list(
    genes_by_pathway = genes_by_pathway,
    pathway_names = pathway_names[names(genes_by_pathway)],
    source = url
  )
}

pmn_org_candidates <- function(cyc_db) {
  cyc_db <- trimws(as.character(cyc_db %||% ""))
  no_cyc <- sub("Cyc$", "", cyc_db, ignore.case = FALSE)
  unique(c(cyc_db, no_cyc, toupper(no_cyc), tolower(no_cyc), tolower(cyc_db)))
}

parse_biocyc_object_ids <- function(xml_text, classes = c("Gene", "Pathway")) {
  if (is.null(xml_text) || !nzchar(xml_text)) return(character())
  class_pattern <- paste(classes, collapse = "|")
  patterns <- c(
    paste0("<(?:", class_pattern, ")[^>]*(?:frameid|ID)=\"([^\"]+)\""),
    paste0("<(?:", class_pattern, ")[^>]*rdf:about=\"[^\"]*#([^\"]+)\"")
  )
  ids <- character()
  for (pat in patterns) {
    m <- gregexpr(pat, xml_text, perl = TRUE, ignore.case = TRUE)
    vals <- regmatches(xml_text, m)[[1]]
    if (length(vals) > 0 && !identical(vals, character(0))) {
      ids <- c(ids, sub(pat, "\\1", vals, perl = TRUE, ignore.case = TRUE))
    }
  }
  unique(ids[!is.na(ids) & nzchar(ids)])
}

parse_biocyc_pathway_names <- function(xml_text, pathway_ids) {
  out <- stats::setNames(pathway_ids, pathway_ids)
  for (pid in pathway_ids) {
    block_pat <- paste0("<[^>]*(?:frameid|ID)=\"", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", pid), "\"[^>]*>.*?</[^>]+>")
    block_match <- regmatches(xml_text, regexpr(block_pat, xml_text, perl = TRUE, ignore.case = TRUE))
    block <- if (length(block_match) > 0) block_match[[1]] else NA_character_
    if (!is.na(block) && nzchar(block)) {
      name_match <- regmatches(block, regexpr("<common-name[^>]*>.*?</common-name>", block, perl = TRUE, ignore.case = TRUE))
      name <- if (length(name_match) > 0) name_match[[1]] else NA_character_
      if (!is.na(name) && nzchar(name)) {
        name <- gsub("<[^>]+>", "", name)
        name <- trimws(utils::URLdecode(name))
        name <- gsub("&amp;", "&", name, fixed = TRUE)
        if (nzchar(name)) out[pid] <- name
      }
    }
  }
  out
}

pmn_read_url <- function(url) {
  con <- url(url, open = "rt")
  on.exit(close(con), add = TRUE)
  paste(readLines(con, warn = FALSE), collapse = "\n")
}

load_pmn_pathways_from_webservice <- function(cyc_db) {
  hosts <- c("https://pmn.plantcyc.org", "https://plantcyc.org", "https://websvc.biocyc.org")
  last_error <- NULL
  for (host in hosts) {
    for (org in pmn_org_candidates(cyc_db)) {
      query <- paste0("[x:x<-", tolower(org), "^^pathways]")
      url_all <- paste0(host, "/xmlquery?query=", utils::URLencode(query, reserved = TRUE), "&detail=low")
      xml <- tryCatch(pmn_read_url(url_all), error = function(e) { last_error <<- e$message; NULL })
      pathway_ids <- parse_biocyc_object_ids(xml, classes = "Pathway")
      if (length(pathway_ids) == 0) next
      pathway_names <- parse_biocyc_pathway_names(xml, pathway_ids)
      genes_by_pathway <- list()
      for (pid in pathway_ids) {
        obj_id <- if (grepl(":", pid, fixed = TRUE)) pid else paste0(org, ":", pid)
        gene_url <- paste0(host, "/apixml?fn=genes-of-pathway&id=", utils::URLencode(obj_id, reserved = TRUE), "&detail=none")
        gene_xml <- tryCatch(pmn_read_url(gene_url), error = function(e) NULL)
        genes <- normalize_pmn_gene_ids(parse_biocyc_object_ids(gene_xml, classes = "Gene"))
        if (length(genes) > 0) genes_by_pathway[[pid]] <- genes
      }
      genes_by_pathway <- genes_by_pathway[lengths(genes_by_pathway) > 0]
      if (length(genes_by_pathway) > 0) {
        return(list(
          genes_by_pathway = genes_by_pathway,
          pathway_names = pathway_names[names(genes_by_pathway)],
          source = paste0(host, " / ", org)
        ))
      }
    }
  }
  stop("Could not load PMN/BioCyc pathway genes for ", cyc_db, ". Last error: ", last_error %||% "no pathway genes returned")
}

get_pmn_genes <- function(cyc_db = "AraCyc") {
  cyc_db <- trimws(as.character(cyc_db %||% ""))
  if (!nzchar(cyc_db)) stop("No PMN Cyc database selected.")
  tryCatch(
    load_pmn_pathways_from_tsv(cyc_db),
    error = function(e) {
      stop(
        "Could not load PMN pathway table for ", cyc_db, ". Expected file: ",
        pmn_pathway_file_url(cyc_db), ". Original error: ", e$message
      )
    }
  )
}

run_pmn_enrichment <- function(de_df, direction = c("up", "down", "all"), cyc_db = "AraCyc",
                               alpha = 0.05, lfc_cutoff = 1, min_set_size = 3,
                               p_adjust_method = "BH") {
  direction <- match.arg(direction)
  de_df <- classify_de(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  universe <- unique(normalize_pmn_gene_ids(de_df$gene_id))
  if (length(universe) < 5) stop("Too few background genes after PMN ID normalization.")

  if (direction == "up") interesting <- de_df$gene_id[de_df$DE_class == "up"]
  if (direction == "down") interesting <- de_df$gene_id[de_df$DE_class == "down"]
  if (direction == "all") interesting <- de_df$gene_id[de_df$DE_class %in% c("up", "down")]
  interesting <- intersect(unique(normalize_pmn_gene_ids(interesting)), universe)
  if (length(interesting) < 2) stop("Too few significant genes for PMN enrichment in this direction.")

  pmn_data <- get_pmn_genes(cyc_db = cyc_db)
  sets <- pmn_data$genes_by_pathway
  pathway_names <- pmn_data$pathway_names
  pmn_universe <- unique(unlist(sets, use.names = FALSE))
  overlap_n <- length(intersect(universe, pmn_universe))
  if (overlap_n == 0) {
    stop("No overlap between loaded gene IDs and ", cyc_db, " pathway genes. PMN usually expects the organism genome locus IDs used by that Cyc database.")
  }

  rows <- lapply(names(sets), function(pathway_id) {
    set_genes <- intersect(unique(sets[[pathway_id]]), universe)
    set_size <- length(set_genes)
    if (set_size < min_set_size) return(NULL)
    a <- length(intersect(set_genes, interesting))
    if (a == 0) return(NULL)
    b <- set_size - a
    c <- length(setdiff(interesting, set_genes))
    d <- max(length(universe) - a - b - c, 0)
    ft <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = "greater")
    pname <- unname(pathway_names[pathway_id])
    if (length(pname) == 0 || is.na(pname) || !nzchar(pname)) pname <- pathway_id
    data.frame(
      pathway.code = pathway_id,
      pathway.name = pname,
      PMN_DB = cyc_db,
      Direction = direction,
      Significant = a,
      Annotated = set_size,
      Significant_total = length(interesting),
      Background_total = length(universe),
      PMN_overlap_gene_ids = overlap_n,
      FoldEnrichment = (a / max(length(interesting), 1)) / (set_size / max(length(universe), 1)),
      p.value = ft$p.value,
      Genes = paste(sort(intersect(set_genes, interesting)), collapse = "; "),
      Data_source = pmn_data$source %||% "PMN/BioCyc web services",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out) || nrow(out) == 0) stop("No PMN pathways overlap the selected significant genes after size filtering.")
  out$p.adjusted <- stats::p.adjust(out$p.value, method = p_adjust_method)
  out <- out[order(out$p.adjusted, out$p.value), , drop = FALSE]
  rownames(out) <- NULL
  out
}

plot_pmn_bubble <- function(pmn_df, p_value_threshold = 0.05, top_n = 20, point_alpha = 0.82,
                            color_up = "#B2182B", color_down = "#2166AC", color_all = "#e6ac6a",
                            plot_theme = "classic", font_family = "serif") {
  if (is.null(pmn_df) || nrow(pmn_df) == 0) return(NULL)
  plot_df <- pmn_df[!is.na(pmn_df$p.adjusted) & pmn_df$p.adjusted <= p_value_threshold, , drop = FALSE]
  if (nrow(plot_df) == 0) return(NULL)
  plot_df <- head(plot_df[order(plot_df$p.adjusted, plot_df$p.value), , drop = FALSE], top_n)
  plot_df$pathway.short <- stringr::str_trunc(plot_df$pathway.name, 70)
  plot_df$neg_log10_padj <- -log10(pmax(plot_df$p.adjusted, .Machine$double.xmin))
  direction <- unique(tolower(plot_df$Direction))
  high_col <- color_up
  if (length(direction) == 1 && direction == "down") high_col <- color_down
  if (length(direction) == 1 && direction == "all") high_col <- color_all
  ggplot2::ggplot(plot_df, ggplot2::aes(x = FoldEnrichment, y = stats::reorder(pathway.short, FoldEnrichment))) +
    ggplot2::geom_point(ggplot2::aes(size = Significant, color = neg_log10_padj), alpha = point_alpha) +
    ggplot2::scale_color_gradient(low = "#979797", high = high_col) +
    plot_theme_choice(plot_theme, base_size = 12, font_family = font_family) +
    ggplot2::labs(x = "Fold enrichment", y = NULL, size = "Genes", color = "-log10(FDR)", title = "PMN pathway enrichment")
}

read_gene_family_table_utf8 <- function(path) {
  read_raw_text <- function(path) {
    if (is_remote_path(path)) {
      con <- url(path, open = "rb")
      on.exit(close(con), add = TRUE)
      raw <- readBin(con, what = "raw", n = 50 * 1024 * 1024)
    } else {
      con <- file(path, open = "rb")
      on.exit(close(con), add = TRUE)
      raw <- readBin(con, what = "raw", n = file.info(path)$size)
    }
    rawToChar(raw)
  }

  tbl <- tryCatch(
    suppressWarnings(read_any_table(path, delim = "\t")),
    error = function(e) NULL
  )
  if (is.null(tbl)) {
    txt <- read_raw_text(path)
    txt <- iconv(txt, from = "latin1", to = "UTF-8", sub = "")
    tbl <- utils::read.table(text = txt, sep = "\t", header = TRUE, check.names = FALSE,
                             quote = "", comment.char = "", fill = TRUE)
  }
  tbl <- as.data.frame(tbl, check.names = FALSE)
  names(tbl) <- iconv(names(tbl), from = "", to = "UTF-8", sub = "")
  chr_cols <- vapply(tbl, is.character, logical(1))
  tbl[chr_cols] <- lapply(tbl[chr_cols], function(x) iconv(as.character(x), from = "", to = "UTF-8", sub = ""))
  tbl
}

load_at_gene_families <- function(path = DEFAULT_AT_GENE_FAMILIES_URL) {
  tbl <- read_gene_family_table_utf8(path)
  locus_col <- first_existing_col(tbl, c("Genomic_Locus_Tag", "gene_id", "GeneID", "Gene.ID", "Locus"))
  family_col <- first_existing_col(tbl, c("Gene_Family", "family", "GeneFamily"))
  subfamily_col <- first_existing_col(tbl, c("Sub_Family", "subfamily", "SubFamily"))
  if (is.null(locus_col) || is.null(family_col)) {
    stop("Gene family file must contain Genomic_Locus_Tag and Gene_Family columns.")
  }
  if (is.null(subfamily_col)) tbl$Sub_Family <- "" else names(tbl)[names(tbl) == subfamily_col] <- "Sub_Family"
  names(tbl)[names(tbl) == locus_col] <- "gene_id"
  names(tbl)[names(tbl) == family_col] <- "Gene_Family"
  tbl$gene_id <- toupper(trimws(as.character(tbl$gene_id)))
  tbl$Gene_Family <- trimws(as.character(tbl$Gene_Family))
  tbl$Sub_Family <- trimws(as.character(tbl$Sub_Family))
  tbl <- tbl[!is.na(tbl$gene_id) & nzchar(tbl$gene_id) &
               !is.na(tbl$Gene_Family) & nzchar(tbl$Gene_Family), , drop = FALSE]
  tbl <- tbl[!duplicated(tbl[, c("gene_id", "Gene_Family", "Sub_Family"), drop = FALSE]), , drop = FALSE]
  rownames(tbl) <- NULL
  tbl
}

setup_at_gene_family_analysis <- function(de_df, alpha = 0.05, lfc_cutoff = 1,
                                          family_file = DEFAULT_AT_GENE_FAMILIES_URL) {
  fam <- load_at_gene_families(family_file)
  de_df <- classify_de(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  de_df$family_gene_id <- toupper(trimws(as.character(de_df$gene_id)))
  de_df$gene_id <- de_df$family_gene_id
  fam$family_gene_id <- fam$gene_id
  merged <- merge(de_df, fam[, c("family_gene_id", "Gene_Family", "Sub_Family"), drop = FALSE],
                  by = "family_gene_id", all.x = FALSE, sort = FALSE)
  merged$family_gene_id <- NULL
  rownames(merged) <- NULL
  families <- sort(unique(fam$Gene_Family))
  list(
    families = fam,
    merged = merged,
    choices = families,
    alpha = alpha,
    lfc_cutoff = lfc_cutoff,
    source = family_file
  )
}

compute_at_gene_family <- function(family_name, ctx) {
  if (is.null(ctx) || is.null(ctx$merged)) stop("Gene family context has not been built.")
  family_name <- trimws(as.character(family_name %||% ""))
  family_name <- family_name[!is.na(family_name) & nzchar(family_name)]
  if (length(family_name) == 0) stop("No gene family selected.")
  out <- ctx$merged[ctx$merged$Gene_Family %in% family_name, , drop = FALSE]
  out <- out[order(out$padj, out$pValue, na.last = TRUE), , drop = FALSE]
  out <- out[!duplicated(out$gene_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

make_at_gene_family_volcano_plot <- function(family_df, family_name, alpha = 0.05, lfc_cutoff = 1,
                                             color_up = "#B2182B", color_down = "#2166AC", color_ns = "grey70",
                                             plot_theme = "classic", font_family = "serif") {
  if (is.null(family_df) || nrow(family_df) == 0) stop("No genes found for selected family.")
  make_gene_group_volcano_plot(
    family_df,
    family_name,
    alpha = alpha,
    lfc_cutoff = lfc_cutoff,
    color_up = color_up,
    color_down = color_down,
    color_ns = color_ns,
    plot_theme = plot_theme,
    font_family = font_family
  )
}

run_at_gene_family_enrichment <- function(de_df, direction = c("up", "down", "all"),
                                          alpha = 0.05, lfc_cutoff = 1,
                                          family_file = DEFAULT_AT_GENE_FAMILIES_URL,
                                          min_set_size = 3, p_adjust_method = "BH") {
  direction <- match.arg(direction)
  de_df <- classify_de(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  de_df$family_gene_id <- toupper(trimws(as.character(de_df$gene_id)))
  universe <- unique(de_df$family_gene_id[!is.na(de_df$family_gene_id) & nzchar(de_df$family_gene_id)])
  if (direction == "up") interesting <- de_df$family_gene_id[de_df$DE_class == "up"]
  if (direction == "down") interesting <- de_df$family_gene_id[de_df$DE_class == "down"]
  if (direction == "all") interesting <- de_df$family_gene_id[de_df$DE_class %in% c("up", "down")]
  interesting <- intersect(unique(interesting[!is.na(interesting) & nzchar(interesting)]), universe)
  if (length(interesting) < 2) stop("Too few significant genes for gene-family enrichment.")

  fam <- load_at_gene_families(family_file)
  sets <- split(fam$gene_id, fam$Gene_Family)
  rows <- lapply(names(sets), function(family_name) {
    set_genes <- intersect(unique(sets[[family_name]]), universe)
    set_size <- length(set_genes)
    if (set_size < min_set_size) return(NULL)
    a <- length(intersect(set_genes, interesting))
    if (a == 0) return(NULL)
    b <- set_size - a
    c <- length(setdiff(interesting, set_genes))
    d <- max(length(universe) - a - b - c, 0)
    ft <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = "greater")
    data.frame(
      Gene_Family = family_name,
      Direction = direction,
      Significant = a,
      Annotated = set_size,
      Significant_total = length(interesting),
      Background_total = length(universe),
      FoldEnrichment = (a / max(length(interesting), 1)) / (set_size / max(length(universe), 1)),
      p.value = ft$p.value,
      Genes = paste(sort(intersect(set_genes, interesting)), collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out) || nrow(out) == 0) stop("No gene families overlap the selected significant genes after size filtering.")
  out$p.adjusted <- stats::p.adjust(out$p.value, method = p_adjust_method)
  out <- out[order(out$p.adjusted, out$p.value), , drop = FALSE]
  rownames(out) <- NULL
  out
}

plot_at_gene_family_enrichment <- function(enrichment_df, p_value_threshold = 0.05, top_n = 10,
                                           point_alpha = 0.82, color_up = "#B2182B",
                                           color_down = "#2166AC", color_all = "#e6ac6a",
                                           plot_theme = "classic", font_family = "serif") {
  if (is.null(enrichment_df) || nrow(enrichment_df) == 0) return(NULL)
  if (!"p.adjusted" %in% names(enrichment_df)) {
    if (!"p.value" %in% names(enrichment_df)) return(NULL)
    enrichment_df$p.adjusted <- enrichment_df$p.value
  }
  plot_df <- enrichment_df[!is.na(enrichment_df$p.adjusted) & enrichment_df$p.adjusted <= p_value_threshold, , drop = FALSE]
  title <- "Gene family enrichment"
  if (nrow(plot_df) == 0) {
    plot_df <- enrichment_df[!is.na(enrichment_df$p.adjusted), , drop = FALSE]
    if (nrow(plot_df) == 0) return(NULL)
    title <- paste0("Top gene family enrichment (none pass FDR <= ", p_value_threshold, ")")
  }
  plot_df <- head(plot_df[order(plot_df$p.adjusted, plot_df$p.value), , drop = FALSE], top_n)
  plot_df$Family_short <- stringr::str_trunc(plot_df$Gene_Family, 70)
  plot_df$neg_log10_padj <- -log10(pmax(plot_df$p.adjusted, .Machine$double.xmin))
  direction <- unique(tolower(plot_df$Direction))
  high_col <- color_up
  if (length(direction) == 1 && direction == "down") high_col <- color_down
  if (length(direction) == 1 && direction == "all") high_col <- color_all
  ggplot2::ggplot(plot_df, ggplot2::aes(x = FoldEnrichment, y = stats::reorder(Family_short, FoldEnrichment))) +
    ggplot2::geom_point(ggplot2::aes(size = Significant, color = neg_log10_padj), alpha = point_alpha) +
    ggplot2::scale_color_gradient(low = "#979797", high = high_col) +
    plot_theme_choice(plot_theme, base_size = 12, font_family = font_family) +
    ggplot2::labs(x = "Fold enrichment", y = NULL, size = "Genes", color = "-log10(FDR)", title = title)
}

normalize_hgnc_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.0$", "", x)
  x <- ifelse(!is.na(x) & nzchar(x) & !grepl("^HGNC:", x, ignore.case = TRUE), paste0("HGNC:", x), x)
  toupper(x)
}

load_hgnc_gene_families <- function(family_url = DEFAULT_HGNC_FAMILY_URL,
                                    gene_has_family_url = DEFAULT_HGNC_GENE_HAS_FAMILY_URL) {
  fam <- as.data.frame(suppressWarnings(read_any_table(family_url, delim = ",")), check.names = FALSE)
  rel <- as.data.frame(suppressWarnings(read_any_table(gene_has_family_url, delim = ",")), check.names = FALSE)
  names(fam) <- trimws(names(fam))
  names(rel) <- trimws(names(rel))
  family_id_col <- first_existing_col(fam, c("family_id", "id", "Family ID", "family.id"))
  rel_family_id_col <- first_existing_col(rel, c("family_id", "Family ID", "family.id"))
  rel_hgnc_col <- first_existing_col(rel, c("hgnc_id", "HGNC ID", "hgnc.id", "gene_id"))
  if (!"name" %in% names(fam)) {
    stop("HGNC family.csv is missing the required family name column: name.")
  }
  if (is.null(family_id_col) || is.null(rel_family_id_col) || is.null(rel_hgnc_col)) {
    stop("HGNC family tables are missing required family_id/name or hgnc_id/family_id columns.")
  }
  fam_sub <- data.frame(
    family_id = trimws(as.character(fam[[family_id_col]])),
    Gene_Family = trimws(as.character(fam[["name"]])),
    stringsAsFactors = FALSE
  )
  rel_sub <- data.frame(
    hgnc_id = normalize_hgnc_id(rel[[rel_hgnc_col]]),
    family_id = trimws(as.character(rel[[rel_family_id_col]])),
    stringsAsFactors = FALSE
  )
  out <- merge(rel_sub, fam_sub, by = "family_id", all.x = TRUE, sort = FALSE)
  out <- out[!is.na(out$hgnc_id) & nzchar(out$hgnc_id) &
               !is.na(out$Gene_Family) & nzchar(out$Gene_Family), , drop = FALSE]
  out$Sub_Family <- ""
  out$gene_id <- out$hgnc_id
  out <- out[!duplicated(out[, c("hgnc_id", "family_id", "Gene_Family"), drop = FALSE]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

map_gene_ids_to_hgnc_from_complete_set <- function(gene_ids, gene_id_type = "ENSEMBL",
                                                   complete_set_url = DEFAULT_HGNC_COMPLETE_SET_URL) {
  hgnc <- as.data.frame(suppressWarnings(read_any_table(complete_set_url, delim = "\t")), check.names = FALSE)
  names(hgnc) <- trimws(names(hgnc))
  keytype <- toupper(gene_id_type %||% "ENSEMBL")
  original <- as.character(gene_ids)
  lookup <- trimws(original)
  if (keytype == "ENSEMBL") lookup <- sub("\\.[0-9]+$", "", lookup)
  key_col <- switch(keytype,
    ENSEMBL = first_existing_col(hgnc, c("ensembl_gene_id", "ensembl_id", "ensembl")),
    SYMBOL = first_existing_col(hgnc, c("symbol", "approved_symbol")),
    ENTREZID = first_existing_col(hgnc, c("entrez_id", "entrezgene_id", "entrez")),
    REFSEQ = first_existing_col(hgnc, c("refseq_accession", "refseq_ids", "refseq_id", "refseq")),
    UNIPROT = first_existing_col(hgnc, c("uniprot_ids", "uniprot_id")),
    NULL
  )
  hgnc_col <- first_existing_col(hgnc, c("hgnc_id", "HGNC ID", "hgnc.id"))
  if (is.null(key_col) || is.null(hgnc_col)) {
    return(data.frame(gene_id = original, hgnc_id = NA_character_, stringsAsFactors = FALSE))
  }
  map <- data.frame(
    lookup_id = as.character(hgnc[[key_col]]),
    hgnc_id = normalize_hgnc_id(hgnc[[hgnc_col]]),
    stringsAsFactors = FALSE
  )
  map <- map[!is.na(map$lookup_id) & nzchar(map$lookup_id) &
               !is.na(map$hgnc_id) & nzchar(map$hgnc_id), , drop = FALSE]
  if (keytype %in% c("UNIPROT", "REFSEQ")) {
    map <- do.call(rbind, lapply(seq_len(nrow(map)), function(i) {
      ids <- trimws(unlist(strsplit(map$lookup_id[i], "\\||,|;")))
      ids <- ids[nzchar(ids)]
      data.frame(lookup_id = ids, hgnc_id = map$hgnc_id[i], stringsAsFactors = FALSE)
    }))
  }
  map$lookup_id <- trimws(as.character(map$lookup_id))
  if (keytype %in% c("ENSEMBL", "SYMBOL")) map$lookup_id <- toupper(map$lookup_id)
  lookup_key <- lookup
  if (keytype %in% c("ENSEMBL", "SYMBOL")) lookup_key <- toupper(lookup_key)
  map <- map[!duplicated(map$lookup_id), , drop = FALSE]
  out <- merge(data.frame(gene_id = original, lookup_id = lookup_key, stringsAsFactors = FALSE),
               map, by = "lookup_id", all.x = TRUE, sort = FALSE)
  out[, c("gene_id", "hgnc_id"), drop = FALSE]
}

map_gene_ids_to_hgnc <- function(gene_ids, gene_id_type = "ENSEMBL", orgdb = "org.Hs.eg.db") {
  original <- as.character(gene_ids)
  lookup <- trimws(original)
  keytype <- toupper(gene_id_type %||% "ENSEMBL")
  if (keytype == "ENSEMBL") lookup <- sub("\\.[0-9]+$", "", lookup)

  if (keytype %in% c("HGNC", "HGNC_ID")) {
    return(data.frame(gene_id = original, hgnc_id = normalize_hgnc_id(lookup), stringsAsFactors = FALSE))
  }

  hgnc_complete <- tryCatch(
    map_gene_ids_to_hgnc_from_complete_set(original, gene_id_type = keytype),
    error = function(e) NULL
  )
  if (!is.null(hgnc_complete) && any(!is.na(hgnc_complete$hgnc_id) & nzchar(hgnc_complete$hgnc_id))) {
    return(hgnc_complete)
  }

  out <- data.frame(gene_id = original, lookup_id = lookup, hgnc_id = NA_character_, stringsAsFactors = FALSE)
  if (requireNamespace("biomaRt", quietly = TRUE)) {
    mart_map <- tryCatch({
      mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
      filter <- switch(keytype,
        ENSEMBL = "ensembl_gene_id",
        SYMBOL = "hgnc_symbol",
        ENTREZID = "entrezgene_id",
        REFSEQ = NA_character_,
        UNIPROT = NA_character_,
        "ensembl_gene_id"
      )
      if (is.na(filter)) stop("biomaRt fallback is not configured for ", keytype)
      vals <- unique(out$lookup_id[!is.na(out$lookup_id) & nzchar(out$lookup_id)])
      bm <- biomaRt::getBM(attributes = c(filter, "hgnc_id"), filters = filter, values = vals, mart = mart)
      names(bm)[names(bm) == filter] <- "lookup_id"
      bm$lookup_id <- as.character(bm$lookup_id)
      bm$hgnc_id <- normalize_hgnc_id(bm$hgnc_id)
      bm <- bm[!is.na(bm$hgnc_id) & nzchar(bm$hgnc_id), c("lookup_id", "hgnc_id"), drop = FALSE]
      bm[!duplicated(bm$lookup_id), , drop = FALSE]
    }, error = function(e) NULL)
    if (!is.null(mart_map) && nrow(mart_map) > 0) {
      out <- merge(out, mart_map, by = "lookup_id", all.x = TRUE, sort = FALSE, suffixes = c("", ".bm"))
      out$hgnc_id <- out$hgnc_id.bm
      out$hgnc_id.bm <- NULL
    }
  }

  if (all(is.na(out$hgnc_id)) && requireNamespace("AnnotationDbi", quietly = TRUE) && requireNamespace(orgdb, quietly = TRUE)) {
    orgdb_obj <- getExportedValue(orgdb, orgdb)
    db_keytype <- keytype
    if (db_keytype == "HGNC_ID") db_keytype <- "HGNC"
    if (db_keytype %in% AnnotationDbi::keytypes(orgdb_obj) && "HGNC" %in% AnnotationDbi::columns(orgdb_obj)) {
      db_map <- tryCatch({
        AnnotationDbi::select(orgdb_obj, keys = unique(out$lookup_id), keytype = db_keytype, columns = "HGNC")
      }, error = function(e) NULL)
      if (!is.null(db_map) && nrow(db_map) > 0) {
        names(db_map)[names(db_map) == db_keytype] <- "lookup_id"
        names(db_map)[names(db_map) == "HGNC"] <- "hgnc_id"
        db_map$lookup_id <- as.character(db_map$lookup_id)
        db_map$hgnc_id <- normalize_hgnc_id(db_map$hgnc_id)
        db_map <- db_map[!is.na(db_map$hgnc_id) & nzchar(db_map$hgnc_id), c("lookup_id", "hgnc_id"), drop = FALSE]
        db_map <- db_map[!duplicated(db_map$lookup_id), , drop = FALSE]
        out <- merge(out[, c("gene_id", "lookup_id"), drop = FALSE], db_map, by = "lookup_id", all.x = TRUE, sort = FALSE)
      }
    }
  }
  out[, c("gene_id", "hgnc_id"), drop = FALSE]
}

setup_hgnc_gene_family_analysis <- function(de_df, alpha = 0.05, lfc_cutoff = 1,
                                            gene_id_type = "ENSEMBL", orgdb = "org.Hs.eg.db",
                                            family_url = DEFAULT_HGNC_FAMILY_URL,
                                            gene_has_family_url = DEFAULT_HGNC_GENE_HAS_FAMILY_URL) {
  hgnc_fam <- load_hgnc_gene_families(family_url, gene_has_family_url)
  de_df <- classify_de(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  de_df$original_gene_id <- as.character(de_df$gene_id)
  id_map <- map_gene_ids_to_hgnc(de_df$gene_id, gene_id_type = gene_id_type, orgdb = orgdb)
  de_df$join_row <- seq_len(nrow(de_df))
  mapped <- merge(de_df, id_map, by.x = "original_gene_id", by.y = "gene_id", all.x = TRUE, sort = FALSE)
  mapped <- mapped[order(mapped$join_row), , drop = FALSE]
  mapped$join_row <- NULL
  merged <- merge(mapped, hgnc_fam[, c("hgnc_id", "family_id", "Gene_Family", "Sub_Family"), drop = FALSE],
                  by = "hgnc_id", all.x = FALSE, sort = FALSE)
  if (nrow(merged) == 0) {
    mapped_n <- sum(!is.na(mapped$hgnc_id) & nzchar(mapped$hgnc_id))
    stop("No loaded genes could be matched to HGNC gene families. Mapped ",
         mapped_n, " input IDs to hgnc_id. Check that Gene ID type is correct for the current human DE table.")
  }
  merged$gene_id <- merged$original_gene_id
  rownames(merged) <- NULL
  list(
    families = hgnc_fam,
    merged = merged,
    choices = sort(unique(hgnc_fam$Gene_Family)),
    alpha = alpha,
    lfc_cutoff = lfc_cutoff,
    source = paste(family_url, gene_has_family_url, sep = " + "),
    backend = "human_hgnc"
  )
}

gene_family_backend_for_tax <- function(tax_id) {
  tax_id <- suppressWarnings(as.integer(tax_id))
  if (identical(tax_id, 3702L)) return("arabidopsis")
  if (identical(tax_id, 9606L)) return("human_hgnc")
  NA_character_
}

setup_gene_family_analysis <- function(de_df, tax_id = 3702, alpha = 0.05, lfc_cutoff = 1,
                                       gene_id_type = "TAIR", orgdb = "org.At.tair.db") {
  backend <- gene_family_backend_for_tax(tax_id)
  if (identical(backend, "arabidopsis")) {
    ctx <- setup_at_gene_family_analysis(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
    ctx$backend <- backend
    return(ctx)
  }
  if (identical(backend, "human_hgnc")) {
    return(setup_hgnc_gene_family_analysis(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff,
                                           gene_id_type = gene_id_type, orgdb = orgdb))
  }
  stop("Gene-family analysis is available only for Arabidopsis thaliana and Homo sapiens.")
}

run_gene_family_enrichment <- function(de_df, tax_id = 3702, direction = c("up", "down", "all"),
                                       alpha = 0.05, lfc_cutoff = 1, min_set_size = 3,
                                       gene_id_type = "TAIR", orgdb = "org.At.tair.db",
                                       p_adjust_method = "BH", ctx = NULL) {
  backend <- gene_family_backend_for_tax(tax_id)
  if (identical(backend, "arabidopsis")) {
    if (!is.null(ctx) && identical(ctx$backend %||% "arabidopsis", "arabidopsis")) {
      de_df <- ctx$merged
      de_df$family_gene_id <- toupper(trimws(as.character(de_df$gene_id)))
      de_df <- classify_de(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
      direction <- match.arg(direction)
      universe <- unique(de_df$family_gene_id[!is.na(de_df$family_gene_id) & nzchar(de_df$family_gene_id)])
      if (direction == "up") interesting <- de_df$family_gene_id[de_df$DE_class == "up"]
      if (direction == "down") interesting <- de_df$family_gene_id[de_df$DE_class == "down"]
      if (direction == "all") interesting <- de_df$family_gene_id[de_df$DE_class %in% c("up", "down")]
      interesting <- intersect(unique(interesting[!is.na(interesting) & nzchar(interesting)]), universe)
      fam <- ctx$families
      sets <- split(fam$gene_id, fam$Gene_Family)
    } else {
      return(run_at_gene_family_enrichment(de_df, direction = direction, alpha = alpha,
                                           lfc_cutoff = lfc_cutoff, min_set_size = min_set_size,
                                           p_adjust_method = p_adjust_method))
    }
    if (length(interesting) < 2) stop("Too few significant genes for gene-family enrichment.")
    rows <- lapply(names(sets), function(family_name) {
      set_genes <- intersect(unique(sets[[family_name]]), universe)
      set_size <- length(set_genes)
      if (set_size < min_set_size) return(NULL)
      a <- length(intersect(set_genes, interesting))
      if (a == 0) return(NULL)
      b <- set_size - a
      c <- length(setdiff(interesting, set_genes))
      d <- max(length(universe) - a - b - c, 0)
      ft <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = "greater")
      data.frame(Gene_Family = family_name, Direction = direction, Significant = a,
                 Annotated = set_size, Significant_total = length(interesting),
                 Background_total = length(universe),
                 FoldEnrichment = (a / max(length(interesting), 1)) / (set_size / max(length(universe), 1)),
                 p.value = ft$p.value,
                 Genes = paste(sort(intersect(set_genes, interesting)), collapse = "; "),
                 stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
    if (is.null(out) || nrow(out) == 0) stop("No gene families overlap the selected significant genes after size filtering.")
    out$p.adjusted <- stats::p.adjust(out$p.value, method = p_adjust_method)
    out <- out[order(out$p.adjusted, out$p.value), , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  if (identical(backend, "human_hgnc")) {
    if (is.null(ctx) || !identical(ctx$backend, "human_hgnc")) {
      ctx <- setup_hgnc_gene_family_analysis(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff,
                                             gene_id_type = gene_id_type, orgdb = orgdb)
    }
    universe <- unique(ctx$merged$hgnc_id[!is.na(ctx$merged$hgnc_id) & nzchar(ctx$merged$hgnc_id)])
    df <- classify_de(ctx$merged, alpha = alpha, lfc_cutoff = lfc_cutoff)
    direction <- match.arg(direction)
    if (direction == "up") interesting <- df$hgnc_id[df$DE_class == "up"]
    if (direction == "down") interesting <- df$hgnc_id[df$DE_class == "down"]
    if (direction == "all") interesting <- df$hgnc_id[df$DE_class %in% c("up", "down")]
    interesting <- intersect(unique(interesting[!is.na(interesting) & nzchar(interesting)]), universe)
    if (length(interesting) < 2) stop("Too few significant genes with HGNC IDs for gene-family enrichment.")
    sets <- split(ctx$families$hgnc_id, ctx$families$Gene_Family)
    rows <- lapply(names(sets), function(family_name) {
      set_genes <- intersect(unique(sets[[family_name]]), universe)
      set_size <- length(set_genes)
      if (set_size < min_set_size) return(NULL)
      a <- length(intersect(set_genes, interesting))
      if (a == 0) return(NULL)
      b <- set_size - a
      c <- length(setdiff(interesting, set_genes))
      d <- max(length(universe) - a - b - c, 0)
      ft <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = "greater")
      data.frame(
        Gene_Family = family_name,
        Direction = direction,
        Significant = a,
        Annotated = set_size,
        Significant_total = length(interesting),
        Background_total = length(universe),
        FoldEnrichment = (a / max(length(interesting), 1)) / (set_size / max(length(universe), 1)),
        p.value = ft$p.value,
        Genes = paste(sort(unique(df$gene_id[df$hgnc_id %in% intersect(set_genes, interesting)])), collapse = "; "),
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
    if (is.null(out) || nrow(out) == 0) stop("No HGNC gene families overlap the selected significant genes after size filtering.")
    out$p.adjusted <- stats::p.adjust(out$p.value, method = p_adjust_method)
    out <- out[order(out$p.adjusted, out$p.value), , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  stop("Gene-family enrichment is available only for Arabidopsis thaliana and Homo sapiens.")
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

make_go_offspring_summary <- function(de_df, parent_go_ids, alpha = 0.05, lfc_cutoff = 1,
                                      orgdb = "org.At.tair.db", keytype = "TAIR") {
  de_df <- add_go_bp_column(de_df, orgdb = orgdb, keytype = keytype)
  de_df <- classify_de(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  parent_go_ids <- trimws(unlist(strsplit(parent_go_ids, ",|;|\\s+")))
  parent_go_ids <- parent_go_ids[parent_go_ids != ""]
  rows <- lapply(parent_go_ids, function(go_id) {
    terms <- go_offspring_terms(go_id)
    genes <- genes_matching_go_terms(de_df, terms)
    sub <- de_df[de_df$gene_id %in% genes, , drop = FALSE]
    sig <- sub[sub$DE_class %in% c("up", "down"), , drop = FALSE]
    up_ids <- sort(unique(sig$gene_id[sig$DE_class == "up"]))
    down_ids <- sort(unique(sig$gene_id[sig$DE_class == "down"]))
    sig_ids <- sort(unique(sig$gene_id))
    data.frame(
      Parent_GO_ID = go_id,
      Category = go_term_title(go_id),
      Offspring_terms = length(terms),
      Total = length(unique(sub$gene_id)),
      Upregulated = length(up_ids),
      Downregulated = length(down_ids),
      Significant = length(sig_ids),
      Upregulated_gene_ids = paste(up_ids, collapse = "; "),
      Downregulated_gene_ids = paste(down_ids, collapse = "; "),
      Significant_gene_ids = paste(sig_ids, collapse = "; "),
      Percentage = ifelse(length(unique(sub$gene_id)) > 0, paste0(round(100 * length(unique(sig$gene_id)) / length(unique(sub$gene_id)), 2), "%"), "NA"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

parse_pmn_pathway_codes <- function(pathway_codes) {
  pathway_codes <- trimws(unlist(strsplit(as.character(pathway_codes %||% ""), ",|;|\\s+")))
  pathway_codes <- unique(pathway_codes[!is.na(pathway_codes) & nzchar(pathway_codes)])
  if (length(pathway_codes) == 0) stop("Enter at least one PMN pathway code.")
  pathway_codes
}

make_pmn_pathway_summary <- function(de_df, pathway_codes, cyc_db = "AraCyc",
                                     alpha = 0.05, lfc_cutoff = 1, pmn_data = NULL) {
  pathway_codes <- parse_pmn_pathway_codes(pathway_codes)
  de_df <- classify_de(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  de_df$pmn_gene_id <- normalize_pmn_gene_ids(de_df$gene_id)
  de_df <- de_df[!is.na(de_df$pmn_gene_id) & nzchar(de_df$pmn_gene_id), , drop = FALSE]
  if (nrow(de_df) == 0) stop("No usable gene IDs were found in the loaded DE table for PMN matching.")

  if (is.null(pmn_data)) pmn_data <- get_pmn_genes(cyc_db = cyc_db)
  sets <- pmn_data$genes_by_pathway
  pathway_names <- pmn_data$pathway_names
  known_codes <- names(sets)
  rows <- lapply(pathway_codes, function(code) {
    hit <- known_codes[toupper(known_codes) == toupper(code)]
    if (length(hit) == 0) {
      return(data.frame(
        pathway.code = code,
        pathway.name = NA_character_,
        PMN_DB = cyc_db,
        Status = "Pathway code not found",
        Total_pathway_genes = 0,
        Matched_loaded_genes = 0,
        Upregulated = 0,
        Downregulated = 0,
        Significant = 0,
        Upregulated_gene_ids = "",
        Downregulated_gene_ids = "",
        Significant_gene_ids = "",
        stringsAsFactors = FALSE
      ))
    }
    pathway_id <- hit[1]
    pathway_genes <- unique(sets[[pathway_id]])
    sub <- de_df[de_df$pmn_gene_id %in% pathway_genes, , drop = FALSE]
    sig <- sub[sub$DE_class %in% c("up", "down"), , drop = FALSE]
    up_ids <- sort(unique(sig$gene_id[sig$DE_class == "up"]))
    down_ids <- sort(unique(sig$gene_id[sig$DE_class == "down"]))
    sig_ids <- sort(unique(sig$gene_id))
    matched_ids <- sort(unique(sub$gene_id))
    pname <- unname(pathway_names[pathway_id])
    if (length(pname) == 0 || is.na(pname) || !nzchar(pname)) pname <- pathway_id
    data.frame(
      pathway.code = pathway_id,
      pathway.name = pname,
      PMN_DB = cyc_db,
      Status = "Found",
      Total_pathway_genes = length(pathway_genes),
      Matched_loaded_genes = length(matched_ids),
      Upregulated = length(up_ids),
      Downregulated = length(down_ids),
      Significant = length(sig_ids),
      Upregulated_gene_ids = paste(up_ids, collapse = "; "),
      Downregulated_gene_ids = paste(down_ids, collapse = "; "),
      Significant_gene_ids = paste(sig_ids, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

make_pmn_pathway_gene_table <- function(de_df, pathway_codes, cyc_db = "AraCyc",
                                        alpha = 0.05, lfc_cutoff = 1, pmn_data = NULL) {
  pathway_codes <- parse_pmn_pathway_codes(pathway_codes)
  de_df <- classify_de(de_df, alpha = alpha, lfc_cutoff = lfc_cutoff)
  de_df$pmn_gene_id <- normalize_pmn_gene_ids(de_df$gene_id)
  de_df <- de_df[!is.na(de_df$pmn_gene_id) & nzchar(de_df$pmn_gene_id), , drop = FALSE]
  if (nrow(de_df) == 0) stop("No usable gene IDs were found in the loaded DE table for PMN matching.")

  if (is.null(pmn_data)) pmn_data <- get_pmn_genes(cyc_db = cyc_db)
  sets <- pmn_data$genes_by_pathway
  pathway_names <- pmn_data$pathway_names
  known_codes <- names(sets)

  rows <- lapply(pathway_codes, function(code) {
    hit <- known_codes[toupper(known_codes) == toupper(code)]
    if (length(hit) == 0) return(NULL)
    pathway_id <- hit[1]
    pathway_genes <- unique(sets[[pathway_id]])
    sub <- de_df[de_df$pmn_gene_id %in% pathway_genes, , drop = FALSE]
    if (nrow(sub) == 0) return(NULL)
    pname <- unname(pathway_names[pathway_id])
    if (length(pname) == 0 || is.na(pname) || !nzchar(pname)) pname <- pathway_id
    sub$pathway.code <- pathway_id
    sub$pathway.name <- pname
    sub$PMN_DB <- cyc_db
    sub$Total_pathway_genes <- length(pathway_genes)
    sub
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out) || nrow(out) == 0) {
    return(data.frame())
  }
  out <- out[order(out$pathway.code, out$padj, out$pValue, na.last = TRUE), , drop = FALSE]
  out <- out[!duplicated(out[, intersect(c("pathway.code", "gene_id"), names(out)), drop = FALSE]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

make_abiotic_stress_table <- function(de_df, dataset = c("all", "up", "down"), alpha = 0.05, lfc_cutoff = 1,
                                      orgdb = "org.At.tair.db", keytype = "TAIR") {
  dataset <- match.arg(dataset)
  de_df <- add_go_bp_column(de_df, orgdb = orgdb, keytype = keytype)
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

make_abiotic_stress_plot <- function(stress_df, title = "Abiotic stress enrichment",
                                     plot_theme = "classic", font_family = "serif") {
  if (is.null(stress_df) || nrow(stress_df) == 0) stop("No stress enrichment table to plot.")
  stress_df$Label <- paste0("p=", signif(stress_df$P_value, 2), "\n", stress_df$Sig_in_term, "/", stress_df$Total_in_term)
  ggplot2::ggplot(stress_df, ggplot2::aes(x = stats::reorder(Test_name, Fold_enrichment), y = Fold_enrichment, fill = Significant)) +
    ggplot2::geom_col(alpha = 0.85) +
    ggplot2::geom_text(ggplot2::aes(label = Label), hjust = -0.05, size = 3) +
    ggplot2::coord_flip() +
    plot_theme_choice(plot_theme, base_size = 12, font_family = font_family) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey70")) +
    ggplot2::labs(x = NULL, y = "Fold enrichment", fill = "p < 0.05", title = title) +
    ggplot2::expand_limits(y = max(stress_df$Fold_enrichment, na.rm = TRUE) * 1.25)
}

plot_revigo_data <- function(plot_df, title = "REVIGO-like semantic reduction",
                             color_palette = "set1", show_labels = TRUE,
                             max_parent_labels = 8, plot_theme = "classic", font_family = "serif") {
  if (is.null(plot_df) || nrow(plot_df) == 0) stop("No REVIGO-like plot data available.")
  plot_df <- plot_df[!is.na(plot_df$plot_X) & !is.na(plot_df$plot_Y), , drop = FALSE]
  if (nrow(plot_df) == 0) stop("REVIGO-like reduction produced no plottable terms.")

  parent_terms <- unique(plot_df$parentTerm)
  cols <- pca_palette_values(length(parent_terms), color_palette)
  if (is.null(cols)) cols <- grDevices::rainbow(length(parent_terms))
  names(cols) <- parent_terms
  plot_df$hex_col <- cols[plot_df$parentTerm]

  label_df <- data.frame()
  if (isTRUE(show_labels)) {
    n_labels <- suppressWarnings(as.integer(max_parent_labels %||% 8))
    if (!is.finite(n_labels) || is.na(n_labels)) n_labels <- 8
    n_labels <- max(0, n_labels)
    parent_counts <- sort(table(plot_df$parentTerm), decreasing = TRUE)
    label_parents <- names(parent_counts)[seq_len(min(n_labels, length(parent_counts)))]
    label_df <- plot_df[plot_df$parentTerm %in% label_parents, , drop = FALSE]
    label_df <- label_df[!duplicated(label_df$parentTerm), , drop = FALSE]
  }

  x_range <- diff(range(plot_df$plot_X, na.rm = TRUE))
  y_range <- diff(range(plot_df$plot_Y, na.rm = TRUE))
  if (!is.finite(x_range) || x_range == 0) x_range <- 1
  if (!is.finite(y_range) || y_range == 0) y_range <- 1
  label_margin <- if (isTRUE(show_labels) && nrow(label_df) > 0) {
    max(nchar(label_df$parentTerm), na.rm = TRUE) * 0.15
  } else {
    0.1
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = plot_X, y = plot_Y, color = parentTerm)) +
    ggplot2::geom_point(ggplot2::aes(size = value), alpha = 0.5, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = cols, guide = "none") +
    ggplot2::geom_point(ggplot2::aes(size = value),
                        shape = 21, fill = "transparent", colour = ggplot2::alpha("black", 0.6),
                        show.legend = FALSE) +
    plot_theme_choice(plot_theme, base_size = 12, font_family = font_family) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed") +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.75),
      axis.text.y = ggplot2::element_text(size = 9),
      axis.text.x = ggplot2::element_text(size = 9),
      axis.line = ggplot2::element_blank(),
      plot.margin = grid::unit(c(0.1, label_margin, 0.1, 0.1), "cm")
    ) +
    ggplot2::labs(y = "semantic space y", x = "semantic space x", title = title) +
    ggplot2::coord_cartesian(
      clip = "off",
      xlim = c(min(plot_df$plot_X, na.rm = TRUE) - x_range / 10, max(plot_df$plot_X, na.rm = TRUE) + x_range / 10),
      ylim = c(min(plot_df$plot_Y, na.rm = TRUE) - y_range / 10, max(plot_df$plot_Y, na.rm = TRUE) + y_range / 10)
    )
  if (isTRUE(show_labels) && nrow(label_df) > 0) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = parentTerm),
      data = label_df,
      direction = "y",
      hjust = 0,
      fontface = "bold",
      size = 3,
      colour = "black",
      xlim = c(max(plot_df$plot_X, na.rm = TRUE) + x_range / 7.5, Inf),
      show.legend = FALSE
    )
  }
  list(plot = p, plot_data = plot_df, legend = data.frame(parentTerm = parent_terms, hex_col = unname(cols), stringsAsFactors = FALSE))
}

make_revigo_scatter_plot <- function(simMatrix, reducedTerms, title = "REVIGO-like semantic reduction",
                                     algorithm = "umap", color_palette = "set1",
                                     show_labels = TRUE, max_parent_labels = 8,
                                     plot_theme = "classic", font_family = "serif") {
  algorithm <- match.arg(tolower(algorithm), c("umap", "pca"))
  scatter_data <- tryCatch({
    rrvgo::scatterPlot(simMatrix, reducedTerms, algorithm = algorithm)$data
  }, error = function(e) {
    if (algorithm == "umap") {
      rrvgo::scatterPlot(simMatrix, reducedTerms, algorithm = "pca")$data
    } else {
      stop(e)
    }
  })
  plot_df <- data.frame(
    term_ID = row.names(scatter_data),
    description = as.character(scatter_data$term),
    plot_X = as.numeric(scatter_data$V1),
    plot_Y = as.numeric(scatter_data$V2),
    value = as.numeric(scatter_data$score),
    parent = as.character(scatter_data$parent),
    parentTerm = as.character(scatter_data$parentTerm),
    stringsAsFactors = FALSE
  )
  plot_df <- plot_df[!is.na(plot_df$plot_X) & !is.na(plot_df$plot_Y), , drop = FALSE]
  if (nrow(plot_df) == 0) stop("REVIGO-like reduction produced no plottable terms.")

  parent_terms <- unique(plot_df$parentTerm)
  cols <- pca_palette_values(length(parent_terms), color_palette)
  if (is.null(cols)) cols <- grDevices::rainbow(length(parent_terms))
  names(cols) <- parent_terms
  plot_df$hex_col <- cols[plot_df$parentTerm]

  label_df <- data.frame()
  if (isTRUE(show_labels)) {
    n_labels <- suppressWarnings(as.integer(max_parent_labels %||% 8))
    if (!is.finite(n_labels) || is.na(n_labels)) n_labels <- 8
    n_labels <- max(0, n_labels)
    parent_counts <- sort(table(plot_df$parentTerm), decreasing = TRUE)
    label_parents <- names(parent_counts)[seq_len(min(n_labels, length(parent_counts)))]
    label_df <- plot_df[plot_df$parentTerm %in% label_parents, , drop = FALSE]
    label_df <- label_df[!duplicated(label_df$parentTerm), , drop = FALSE]
  }

  x_range <- diff(range(plot_df$plot_X, na.rm = TRUE))
  y_range <- diff(range(plot_df$plot_Y, na.rm = TRUE))
  if (!is.finite(x_range) || x_range == 0) x_range <- 1
  if (!is.finite(y_range) || y_range == 0) y_range <- 1
  label_margin <- if (isTRUE(show_labels) && nrow(label_df) > 0) {
    max(nchar(label_df$parentTerm), na.rm = TRUE) * 0.15
  } else {
    0.1
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = plot_X, y = plot_Y, color = parentTerm)) +
    ggplot2::geom_point(ggplot2::aes(size = value), alpha = 0.5, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = cols, guide = "none") +
    ggplot2::geom_point(ggplot2::aes(size = value),
                        shape = 21, fill = "transparent", colour = ggplot2::alpha("black", 0.6),
                        show.legend = FALSE) +
    plot_theme_choice(plot_theme, base_size = 12, font_family = font_family) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed") +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.75),
      axis.text.y = ggplot2::element_text(size = 9),
      axis.text.x = ggplot2::element_text(size = 9),
      axis.line = ggplot2::element_blank(),
      plot.margin = grid::unit(c(0.1, label_margin, 0.1, 0.1), "cm")
    ) +
    ggplot2::labs(y = "semantic space y", x = "semantic space x", title = title) +
    ggplot2::coord_cartesian(
      clip = "off",
      xlim = c(min(plot_df$plot_X, na.rm = TRUE) - x_range / 10, max(plot_df$plot_X, na.rm = TRUE) + x_range / 10),
      ylim = c(min(plot_df$plot_Y, na.rm = TRUE) - y_range / 10, max(plot_df$plot_Y, na.rm = TRUE) + y_range / 10)
    )
  if (isTRUE(show_labels) && nrow(label_df) > 0) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = parentTerm),
      data = label_df,
      direction = "y",
      hjust = 0,
      fontface = "bold",
      size = 3,
      colour = "black",
      xlim = c(max(plot_df$plot_X, na.rm = TRUE) + x_range / 7.5, Inf),
      show.legend = FALSE
    )
  }
  list(plot = p, plot_data = plot_df, legend = data.frame(parentTerm = parent_terms, hex_col = unname(cols), stringsAsFactors = FALSE))
}

run_rrvgo_reduce <- function(go_df, ontology = "BP", top_n = 80, threshold = 0.7, title = "REVIGO-like semantic reduction",
                             orgdb = "org.At.tair.db", plot_theme = "classic", font_family = "serif",
                             algorithm = "umap", color_palette = "set1", show_labels = TRUE,
                             max_parent_labels = 8) {
  required <- c("rrvgo", orgdb, "ggplot2", "ggrepel")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "))
  if (is.null(go_df) || nrow(go_df) < 2) stop("At least two GO terms are required for semantic reduction.")
  go_df <- go_df[order(go_df$pValue_num), , drop = FALSE]
  go_df <- head(go_df, top_n)
  scores <- stats::setNames(-log10(pmax(go_df$pValue_num, .Machine$double.xmin)), go_df$GO.ID)
  simMatrix <- rrvgo::calculateSimMatrix(names(scores), orgdb = orgdb, ont = ontology, method = "Rel")
  reducedTerms <- rrvgo::reduceSimMatrix(simMatrix, scores, threshold = threshold, orgdb = orgdb)
  scatter <- make_revigo_scatter_plot(
    simMatrix, reducedTerms,
    title = title,
    algorithm = algorithm,
    color_palette = color_palette,
    show_labels = show_labels,
    max_parent_labels = max_parent_labels,
    plot_theme = plot_theme,
    font_family = font_family
  )
  list(plot = scatter$plot, table = as.data.frame(reducedTerms), plot_data = scatter$plot_data, legend = scatter$legend)
}

make_revigo_treemap_plot <- function(reduced_terms, title = "REVIGO-like treemap") {
  if (is.null(reduced_terms) || nrow(reduced_terms) == 0) stop("No REVIGO-like reduced terms available.")
  rrvgo::treemapPlot(reduced_terms, size = "score", title = title)
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
