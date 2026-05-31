############################################################
# RefSeq GTF annotation builder
# Builds a description table from an NCBI RefSeq GTF file.
############################################################

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

refseq_gtf_required_packages <- function() {
  required <- c("readr")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required packages: ", paste(missing, collapse = ", "))
  }
}

refseq_gtf_collapse_unique <- function(x, sep = "; ") {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (!length(x)) return(NA_character_)
  x <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  x <- trimws(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  paste(unique(x), collapse = sep)
}

refseq_gtf_first_nonempty <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  x[1]
}

refseq_gtf_clean_db_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "unknown")
}

refseq_gtf_split_ids <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (!length(x)) return(character())
  ids <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  ids <- trimws(ids)
  unique(ids[!is.na(ids) & nzchar(ids)])
}

refseq_gtf_escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

refseq_gtf_extract_attribute <- function(attributes, key) {
  pattern <- paste0("(?:^|;)\\s*", refseq_gtf_escape_regex(key), "\\s+\"([^\"]*)\"")
  if (requireNamespace("stringi", quietly = TRUE)) {
    hit <- stringi::stri_match_first_regex(attributes, pattern)
    return(hit[, 2])
  }
  m <- regexec(pattern, attributes, perl = TRUE)
  r <- regmatches(attributes, m)
  vapply(r, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
}

refseq_gtf_discover_attribute_keys <- function(attributes) {
  pattern <- "(?:^|;)\\s*([A-Za-z0-9_.:-]+)\\s+\""
  if (requireNamespace("stringi", quietly = TRUE)) {
    keys <- unlist(lapply(stringi::stri_match_all_regex(attributes, pattern), function(x) x[, 2]), use.names = FALSE)
  } else {
    matches <- regmatches(attributes, gregexpr(pattern, attributes, perl = TRUE))
    keys <- unlist(lapply(matches, function(x) {
      sub(".*;\\s*", "", sub("\\s+\"$", "", x, perl = TRUE), perl = TRUE)
    }), use.names = FALSE)
  }
  keys <- keys[!is.na(keys) & nzchar(keys)]
  unique(keys)
}

refseq_gtf_extract_db_xrefs <- function(attributes) {
  pattern <- "(?:^|;)\\s*db_xref\\s+\"([^\"]*)\""
  if (requireNamespace("stringi", quietly = TRUE)) {
    matches <- stringi::stri_match_all_regex(attributes, pattern)
    return(lapply(matches, function(x) {
      if (is.null(x) || nrow(x) == 0) character() else x[, 2]
    }))
  }
  raw <- regmatches(attributes, gregexpr(pattern, attributes, perl = TRUE))
  lapply(raw, function(x) {
    if (!length(x)) return(character())
    sub(".*db_xref\\s+\"", "", sub("\"$", "", x, perl = TRUE), perl = TRUE)
  })
}

read_refseq_gtf <- function(gtf_file, n_max = Inf) {
  refseq_gtf_required_packages()
  if (is.null(gtf_file) || !nzchar(gtf_file) || !file.exists(gtf_file)) {
    stop("GTF file not found: ", gtf_file %||% "")
  }

  cols <- c("seqid", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  gtf <- readr::read_tsv(
    gtf_file,
    comment = "#",
    col_names = cols,
    col_types = readr::cols(.default = readr::col_character()),
    n_max = n_max,
    progress = FALSE
  )

  gtf <- as.data.frame(gtf, check.names = FALSE)
  if (nrow(gtf) == 0) stop("No annotation rows were found in the GTF file.")
  gtf$start <- suppressWarnings(as.integer(gtf$start))
  gtf$end <- suppressWarnings(as.integer(gtf$end))
  gtf
}

parse_refseq_gtf_attribute_one <- function(attribute) {
  if (is.na(attribute) || !nzchar(trimws(attribute))) return(list())
  matches <- regmatches(attribute, gregexpr('([A-Za-z0-9_.:-]+) "([^"]*)"', attribute, perl = TRUE))[[1]]
  if (!length(matches)) {
    matches <- regmatches(attribute, gregexpr("([A-Za-z0-9_.:-]+)=([^;]+)", attribute, perl = TRUE))[[1]]
  }
  if (!length(matches) || identical(matches, character(0))) return(list())

  out <- list()
  for (part in matches) {
    if (grepl("=", part, fixed = TRUE)) {
      key <- sub("=.*$", "", part, perl = TRUE)
      value <- sub("^[^=]+=", "", part, perl = TRUE)
    } else {
      key <- sub("\\s.*$", "", part, perl = TRUE)
      value <- sub("^[^[:space:]]+[[:space:]]+\"", "", part, perl = TRUE)
      value <- sub("\"$", "", value, perl = TRUE)
    }
    key <- trimws(key)
    value <- trimws(value)
    if (!nzchar(key) || !nzchar(value)) next
    out[[key]] <- c(out[[key]], value)

    if (identical(key, "db_xref") && grepl(":", value, fixed = TRUE)) {
      db <- sub(":.*$", "", value)
      db_value <- sub("^[^:]+:", "", value)
      db_col <- paste0("db_xref_", refseq_gtf_clean_db_name(db))
      out[[db_col]] <- c(out[[db_col]], db_value)
    }
  }

  if (!length(out)) return(out)
  lapply(out, refseq_gtf_collapse_unique)
}

parse_refseq_gtf_attributes <- function(attributes) {
  parsed <- lapply(attributes, parse_refseq_gtf_attribute_one)
  attr_names <- unique(unlist(lapply(parsed, names), use.names = FALSE))
  if (!length(attr_names)) {
    return(data.frame(.row_id = seq_along(attributes)))
  }

  out <- as.data.frame(
    stats::setNames(lapply(attr_names, function(col) {
      vapply(parsed, function(row) {
        val <- row[[col]]
        if (is.null(val) || is.na(val)) NA_character_ else as.character(val)
      }, character(1))
    }), attr_names),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out$.row_id <- seq_along(attributes)
  out
}

fast_parse_refseq_gtf_attributes <- function(attributes, attribute_keys = NULL) {
  preferred_keys <- c(
    "gene_id", "gene", "locus_tag", "protein_id", "transcript_id",
    "gbkey", "gene_biotype", "gene_synonym", "product", "description",
    "note", "Ontology_term", "transl_table", "exon_number",
    "transcript_biotype", "pseudo", "orig_protein_id", "exception",
    "protein", "transl_except"
  )
  if (is.null(attribute_keys)) {
    discovered_keys <- refseq_gtf_discover_attribute_keys(attributes)
    keys <- unique(c(preferred_keys[preferred_keys %in% discovered_keys], setdiff(discovered_keys, c(preferred_keys, "db_xref"))))
  } else {
    keys <- unique(setdiff(attribute_keys, c("db_xref", grep("^db_xref_", attribute_keys, value = TRUE))))
  }

  out_list <- stats::setNames(lapply(keys, function(key) {
    refseq_gtf_extract_attribute(attributes, key)
  }), keys)

  db_xrefs <- refseq_gtf_extract_db_xrefs(attributes)
  if (length(db_xrefs)) {
    out_list$db_xref <- vapply(db_xrefs, function(x) {
      x <- trimws(as.character(x))
      x <- x[!is.na(x) & nzchar(x)]
      if (!length(x)) return(NA_character_)
      paste(unique(x), collapse = "; ")
    }, character(1))

    db_names <- unique(unlist(lapply(db_xrefs, function(x) {
      x <- x[grepl(":", x, fixed = TRUE)]
      if (!length(x)) return(character())
      refseq_gtf_clean_db_name(sub(":.*$", "", x, perl = TRUE))
    }), use.names = FALSE))
    db_names <- db_names[!is.na(db_names) & nzchar(db_names)]

    for (db in db_names) {
      col <- paste0("db_xref_", db)
      out_list[[col]] <- vapply(db_xrefs, function(x) {
        x <- x[grepl(":", x, fixed = TRUE)]
        if (!length(x)) return(NA_character_)
        db_raw_name <- refseq_gtf_clean_db_name(sub(":.*$", "", x, perl = TRUE))
        vals <- sub("^[^:]+:", "", x[db_raw_name == db], perl = TRUE)
        vals <- trimws(vals)
        vals <- vals[!is.na(vals) & nzchar(vals)]
        if (!length(vals)) NA_character_ else paste(unique(vals), collapse = "; ")
      }, character(1))
    }
  }

  if (!length(out_list)) return(data.frame(.row_id = seq_along(attributes)))
  out <- as.data.frame(out_list, stringsAsFactors = FALSE, check.names = FALSE)
  out$.row_id <- seq_along(attributes)
  out
}

refseq_gtf_coalesce_rows <- function(df, cols) {
  cols <- intersect(cols, names(df))
  out <- rep(NA_character_, nrow(df))
  if (!length(cols)) return(out)
  for (col in cols) {
    vals <- trimws(as.character(df[[col]]))
    hit <- (is.na(out) | !nzchar(out)) & !is.na(vals) & nzchar(vals)
    out[hit] <- vals[hit]
  }
  out
}

refseq_gtf_id_source_table <- function(gtf_file, max_rows = 50000) {
  gtf <- read_refseq_gtf(gtf_file, n_max = max_rows)
  attrs <- fast_parse_refseq_gtf_attributes(
    gtf$attribute,
    attribute_keys = c("gene_id", "gene", "locus_tag", "protein_id", "transcript_id")
  )
  attrs$.row_id <- NULL
  if (!ncol(attrs)) {
    return(data.frame(source = character(), label = character(), non_empty = integer(), unique_values = integer()))
  }

  rows <- lapply(names(attrs), function(col) {
    vals <- refseq_gtf_split_ids(attrs[[col]])
    data.frame(
      source = col,
      label = if (grepl("^db_xref_", col)) paste0("db_xref:", sub("^db_xref_", "", col)) else col,
      non_empty = sum(!is.na(attrs[[col]]) & nzchar(trimws(attrs[[col]]))),
      unique_values = length(vals),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[out$non_empty > 0, , drop = FALSE]

  preferred <- c(
    "gene_id", "gene", "locus_tag", "protein_id", "transcript_id",
    "db_xref_GeneID", "db_xref_GenBank", "db_xref_RefSeq"
  )
  out$rank <- match(out$source, preferred)
  out$rank[is.na(out$rank)] <- length(preferred) + seq_len(sum(is.na(match(out$source, preferred))))
  out <- out[order(out$rank, tolower(out$label)), c("source", "label", "non_empty", "unique_values"), drop = FALSE]
  rownames(out) <- NULL
  out
}

refseq_gtf_id_source_choices <- function(gtf_file, max_rows = 50000) {
  sources <- refseq_gtf_id_source_table(gtf_file, max_rows = max_rows)
  stats::setNames(sources$source, paste0(sources$label, " (", sources$unique_values, " IDs)"))
}

resolve_refseq_gtf_id_source <- function(gene_id_source, available_cols) {
  if (is.null(gene_id_source) || !nzchar(gene_id_source)) {
    preferred <- c(
      "gene_id", "locus_tag", "gene", "db_xref_GeneID",
      "protein_id", "transcript_id", "db_xref_GenBank", "db_xref_RefSeq"
    )
    hit <- preferred[preferred %in% available_cols]
    if (length(hit)) return(hit[1])
    return(available_cols[1])
  }

  if (gene_id_source %in% available_cols) return(gene_id_source)
  normalized <- gsub("^db_xref:", "db_xref_", gene_id_source, ignore.case = TRUE)
  normalized <- paste0("db_xref_", refseq_gtf_clean_db_name(sub("^db_xref_", "", normalized)))
  hit <- available_cols[tolower(available_cols) == tolower(normalized)]
  if (length(hit)) return(hit[1])

  db_hit <- available_cols[tolower(available_cols) == tolower(paste0("db_xref_", refseq_gtf_clean_db_name(gene_id_source)))]
  if (length(db_hit)) return(db_hit[1])

  stop(
    "ID source '", gene_id_source, "' was not found in the GTF attributes. Available sources: ",
    paste(available_cols, collapse = ", ")
  )
}

summarise_refseq_gtf_by_gene <- function(gtf, attrs) {
  attrs$.row_id <- NULL
  row_df <- cbind(gtf[, c("seqid", "source", "feature", "start", "end", "score", "strand", "frame"), drop = FALSE], attrs)

  group_candidates <- c(
    "gene_id", "db_xref_GeneID", "locus_tag", "gene", "protein_id", "transcript_id"
  )
  row_df$gene_group_key <- refseq_gtf_coalesce_rows(row_df, group_candidates)
  missing_key <- is.na(row_df$gene_group_key) | !nzchar(row_df$gene_group_key)
  row_df$gene_group_key[missing_key] <- paste0("gtf_row_", which(missing_key))

  char_cols <- setdiff(names(row_df), c("start", "end", "gene_group_key"))
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::as.data.table(row_df)
    out <- dt[
      ,
      c(
        lapply(.SD, refseq_gtf_collapse_unique),
        list(
          start = {
            value <- suppressWarnings(min(start, na.rm = TRUE))
            if (is.finite(value)) value else NA_integer_
          },
          end = {
            value <- suppressWarnings(max(end, na.rm = TRUE))
            if (is.finite(value)) value else NA_integer_
          }
        )
      ),
      by = gene_group_key,
      .SDcols = char_cols
    ]
    return(as.data.frame(out, check.names = FALSE))
  }

  char_summary <- stats::aggregate(
    row_df[, char_cols, drop = FALSE],
    by = list(gene_group_key = row_df$gene_group_key),
    FUN = refseq_gtf_collapse_unique
  )
  start_summary <- stats::aggregate(
    row_df$start,
    by = list(gene_group_key = row_df$gene_group_key),
    FUN = function(x) {
      out <- suppressWarnings(min(x, na.rm = TRUE))
      if (is.finite(out)) out else NA_integer_
    }
  )
  names(start_summary)[2] <- "start"
  end_summary <- stats::aggregate(
    row_df$end,
    by = list(gene_group_key = row_df$gene_group_key),
    FUN = function(x) {
      out <- suppressWarnings(max(x, na.rm = TRUE))
      if (is.finite(out)) out else NA_integer_
    }
  )
  names(end_summary)[2] <- "end"

  out <- merge(char_summary, start_summary, by = "gene_group_key", all.x = TRUE, sort = FALSE)
  out <- merge(out, end_summary, by = "gene_group_key", all.x = TRUE, sort = FALSE)
  rownames(out) <- NULL
  out
}

build_refseq_gtf_description_file <- function(gtf_file,
                                              gene_id_source = NULL,
                                              output_dir = NULL,
                                              output_file = NULL,
                                              one_row_per_id = TRUE,
                                              input_data = NULL) {
  gtf <- read_refseq_gtf(gtf_file)
  selected_attr_key <- if (!is.null(gene_id_source) && nzchar(gene_id_source) && !grepl("^db_xref[_:]", gene_id_source)) {
    gene_id_source
  } else {
    character()
  }
  attrs <- fast_parse_refseq_gtf_attributes(
    gtf$attribute,
    attribute_keys = unique(c(
      "gene_id", "gene", "locus_tag", "protein_id", "transcript_id",
      "product", "description", "note", "Ontology_term", selected_attr_key
    ))
  )
  attrs_for_choices <- attrs
  attrs_for_choices$.row_id <- NULL
  available_sources <- names(attrs_for_choices)
  if (!length(available_sources)) stop("No parseable GTF attributes were found.")

  id_source <- resolve_refseq_gtf_id_source(gene_id_source, available_sources)
  gene_summary <- summarise_refseq_gtf_by_gene(gtf, attrs)
  gene_summary$selected_gene_id <- gene_summary[[id_source]]
  gene_summary <- gene_summary[!is.na(gene_summary$selected_gene_id) & nzchar(trimws(gene_summary$selected_gene_id)), , drop = FALSE]
  if (nrow(gene_summary) == 0) stop("The selected ID source has no usable IDs: ", id_source)

  if (isTRUE(one_row_per_id)) {
    expanded <- lapply(seq_len(nrow(gene_summary)), function(i) {
      ids <- refseq_gtf_split_ids(gene_summary$selected_gene_id[i])
      if (!length(ids)) return(NULL)
      row <- gene_summary[rep(i, length(ids)), , drop = FALSE]
      row$selected_gene_id <- ids
      row
    })
    gene_summary <- do.call(rbind, expanded)
    rownames(gene_summary) <- NULL
  }

  final_df <- data.frame(
    gene_id = gene_summary$selected_gene_id,
    Symbol = gene_summary$gene %||% NA_character_,
    Protein_name = gene_summary$product %||% NA_character_,
    Function_description = if ("description" %in% names(gene_summary)) {
      gene_summary$description
    } else if ("note" %in% names(gene_summary)) {
      gene_summary$note
    } else {
      gene_summary$product %||% NA_character_
    },
    RefSeq_gene_id = gene_summary$gene_id %||% NA_character_,
    RefSeq_gene_group_key = gene_summary$gene_group_key,
    RefSeq_selected_id_source = id_source,
    locus_tag = gene_summary$locus_tag %||% NA_character_,
    transcript_id = gene_summary$transcript_id %||% NA_character_,
    protein_id = gene_summary$protein_id %||% NA_character_,
    db_xref = gene_summary$db_xref %||% NA_character_,
    seqid = gene_summary$seqid %||% NA_character_,
    start = gene_summary$start,
    end = gene_summary$end,
    strand = gene_summary$strand %||% NA_character_,
    feature = gene_summary$feature %||% NA_character_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if ("db_xref_GO" %in% names(gene_summary)) {
    final_df$GO_terms <- gene_summary$db_xref_GO
  } else if ("Ontology_term" %in% names(gene_summary)) {
    final_df$GO_terms <- gene_summary$Ontology_term
  }

  extra_db_cols <- grep("^db_xref_", names(gene_summary), value = TRUE)
  extra_db_cols <- setdiff(extra_db_cols, c("db_xref_GO", "db_xref_KEGG"))
  for (col in extra_db_cols) {
    out_col <- sub("^db_xref_", "", col)
    if (!out_col %in% names(final_df)) final_df[[out_col]] <- gene_summary[[col]]
  }

  key_fun <- if (exists("gene_join_key", mode = "function")) {
    gene_join_key
  } else {
    function(x) toupper(sub("^.*:", "", sub("\\.[0-9]+$", "", trimws(as.character(x)))))
  }

  final_df$annotation_key <- key_fun(final_df$gene_id)
  final_df <- final_df[!is.na(final_df$annotation_key) & final_df$annotation_key != "", , drop = FALSE]
  final_df <- final_df[!duplicated(final_df$annotation_key), , drop = FALSE]
  final_df$annotation_key <- NULL

  if (!is.null(input_data)) {
    input_data <- as.data.frame(input_data, check.names = FALSE)
    names(input_data)[1] <- "gene_id"
    keep_keys <- key_fun(input_data$gene_id)
    keep_keys <- keep_keys[!is.na(keep_keys) & keep_keys != ""]
    final_df <- final_df[key_fun(final_df$gene_id) %in% keep_keys, , drop = FALSE]
  }

  if (!is.null(output_dir) || !is.null(output_file)) {
    if (is.null(output_file)) {
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      safe_source <- gsub("[^A-Za-z0-9_]+", "_", id_source)
      output_file <- file.path(output_dir, paste0("refseq_gtf_description_", safe_source, "_", Sys.Date(), ".csv"))
    } else {
      dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    }
    utils::write.csv(final_df, output_file, row.names = FALSE)
    attr(final_df, "path") <- output_file
  }

  final_df
}

# Backward-compatible alias for the requested file name typo.
build_refseq_gftf_description_file <- build_refseq_gtf_description_file
