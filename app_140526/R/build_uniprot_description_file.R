############################################################
# UniProt annotation builder
# Downloads one organism's UniProt annotation and maps it to input gene IDs.
############################################################

uniprot_required_packages <- function(extra = character()) {
  required <- unique(c("dplyr", "readr", "tidyr", "stringr", extra))
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required packages: ", paste(missing, collapse = ", "))
  }
}

uniprot_annotation_field_ids <- function() {
  c(
    "accession",
    "id",
    "gene_names",
    "gene_primary",
    "gene_synonym",
    "gene_oln",
    "gene_orf",
    "xref_geneid",
    "xref_refseq",
    "xref_ensembl",
    "organism_name",
    "organism_id",
    "protein_name",
    "cc_function",
    "go_p",
    "go_f",
    "go_c",
    "xref_kegg"
  )
}

download_uniprot_annotation <- function(tax_id = 3702, reviewed_only = FALSE, size = NULL) {
  uniprot_required_packages("readr")
  query <- paste0("organism_id:", tax_id)
  if (isTRUE(reviewed_only)) query <- paste0(query, " AND reviewed:true")

  fields <- paste(uniprot_annotation_field_ids(), collapse = ",")
  endpoint <- if (is.null(size)) "stream" else "search"
  if (!is.null(size)) size <- min(as.integer(size), 500L)
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/", endpoint, "?",
    "query=", utils::URLencode(query, reserved = TRUE),
    "&fields=", fields,
    "&format=tsv"
  )
  if (!is.null(size)) url <- paste0(url, "&size=", as.integer(size))

  message("Downloading UniProt annotation for tax_id = ", tax_id)
  uniprot_df <- readr::read_tsv(url, show_col_types = FALSE, progress = FALSE)
  if (nrow(uniprot_df) == 0) stop("UniProt returned no records for tax_id ", tax_id, ".")
  colnames(uniprot_df) <- make.names(colnames(uniprot_df))

  for (cc in c("Entry", "Entry.Name", "Gene.Names", "Gene.Names..primary.",
               "Gene.Names..synonym.", "Gene.Names..ordered.locus.",
               "Gene.Names..ORF.", "Organism", "Organism..ID.",
               "GeneID", "RefSeq", "Ensembl",
               "Protein.names", "Function..CC.",
               "Gene.Ontology..biological.process.",
               "Gene.Ontology..molecular.function.",
               "Gene.Ontology..cellular.component.",
               "Cross.reference..KEGG.")) {
    if (!cc %in% colnames(uniprot_df)) uniprot_df[[cc]] <- NA_character_
  }

  copy_uniprot_alias <- function(target, candidates) {
    hit <- candidates[candidates %in% colnames(uniprot_df)]
    if (length(hit) == 0) return(invisible(NULL))
    if (!target %in% colnames(uniprot_df) || all(is.na(uniprot_df[[target]]) | uniprot_df[[target]] == "")) {
      uniprot_df[[target]] <<- uniprot_df[[hit[1]]]
    }
    invisible(NULL)
  }
  copy_uniprot_alias("Cross.reference..KEGG.", c("KEGG", "Cross.reference..KEGG."))

  uniprot_df
}

uniprot_id_columns <- function(uniprot_df) {
  annotation_cols <- c(
    "Organism", "Organism..ID.", "Protein.names", "Function..CC.",
    "Gene.Ontology..biological.process.",
    "Gene.Ontology..molecular.function.",
    "Gene.Ontology..cellular.component.",
    "Cross.reference..KEGG."
  )
  setdiff(colnames(uniprot_df), annotation_cols)
}

uniprot_split_ids <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (!length(x)) return(character())
  ids <- unlist(strsplit(x, "[;, ]+"), use.names = FALSE)
  ids <- trimws(ids)
  unique(ids[!is.na(ids) & nzchar(ids) & ids != "-"])
}

uniprot_id_source_table <- function(tax_id = 3702, reviewed_only = FALSE, max_records = 500) {
  uniprot_df <- download_uniprot_annotation(tax_id = tax_id, reviewed_only = reviewed_only, size = max_records)
  id_cols <- uniprot_id_columns(uniprot_df)

  rows <- lapply(id_cols, function(col) {
    vals <- uniprot_split_ids(uniprot_df[[col]])
    data.frame(
      source = col,
      label = col,
      non_empty = sum(!is.na(uniprot_df[[col]]) & nzchar(trimws(as.character(uniprot_df[[col]])))),
      unique_values = length(vals),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[out$non_empty > 0, , drop = FALSE]

  preferred <- c("Entry", "Gene.Names..ordered.locus.", "Gene.Names..primary.", "GeneID", "RefSeq", "Ensembl", "Gene.Names", "Entry.Name")
  out$rank <- match(out$source, preferred)
  out$rank[is.na(out$rank)] <- length(preferred) + seq_len(sum(is.na(match(out$source, preferred))))
  out <- out[order(out$rank, tolower(out$label)), c("source", "label", "non_empty", "unique_values"), drop = FALSE]
  rownames(out) <- NULL
  out
}

uniprot_id_source_choices <- function(tax_id = 3702, reviewed_only = FALSE, max_records = 500) {
  sources <- uniprot_id_source_table(tax_id = tax_id, reviewed_only = reviewed_only, max_records = max_records)
  stats::setNames(sources$source, paste0(sources$label, " (", sources$unique_values, " IDs)"))
}

resolve_uniprot_id_source <- function(gene_id_source, available_cols) {
  if (is.null(gene_id_source) || !nzchar(gene_id_source)) {
    preferred <- c("Gene.Names..ordered.locus.", "Gene.Names..primary.", "GeneID", "RefSeq", "Entry", "Gene.Names")
    hit <- preferred[preferred %in% available_cols]
    if (length(hit)) return(hit[1])
    return(available_cols[1])
  }

  if (gene_id_source %in% available_cols) return(gene_id_source)
  hit <- available_cols[tolower(available_cols) == tolower(make.names(gene_id_source))]
  if (length(hit)) return(hit[1])
  stop(
    "ID source '", gene_id_source, "' was not found in the UniProt table. Available sources: ",
    paste(available_cols, collapse = ", ")
  )
}

search_uniprot_taxonomy <- function(query, size = 50) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Missing required package: readr")
  }

  query <- trimws(as.character(query))
  if (!nzchar(query)) stop("Enter an organism name or taxonomy ID to search.")

  url <- paste0(
    "https://rest.uniprot.org/taxonomy/search?",
    "query=", utils::URLencode(query, reserved = TRUE),
    "&format=tsv",
    "&size=", as.integer(size)
  )

  taxonomy_df <- readr::read_tsv(url, show_col_types = FALSE, progress = FALSE)
  if (nrow(taxonomy_df) == 0) {
    stop("No UniProt taxonomy matches found for: ", query)
  }

  original_names <- names(taxonomy_df)
  names(taxonomy_df) <- make.names(names(taxonomy_df))

  pick_col <- function(candidates) {
    hit <- match(tolower(candidates), tolower(names(taxonomy_df)))
    hit <- hit[!is.na(hit)]
    if (length(hit)) names(taxonomy_df)[hit[1]] else NULL
  }

  id_col <- pick_col(c("Taxon.Id", "Taxonomy.ID", "Organism.ID", "id"))
  scientific_col <- pick_col(c("Scientific.name", "Scientific.Name", "scientificName", "scientific_name"))
  common_col <- pick_col(c("Common.name", "Common.Name", "commonName", "common_name"))
  rank_col <- pick_col(c("Rank", "rank"))

  if (is.null(id_col)) {
    stop("UniProt taxonomy result did not include a taxonomy ID column. Columns returned: ", paste(original_names, collapse = ", "))
  }

  out <- data.frame(
    tax_id = suppressWarnings(as.integer(taxonomy_df[[id_col]])),
    scientific_name = if (!is.null(scientific_col)) as.character(taxonomy_df[[scientific_col]]) else NA_character_,
    common_name = if (!is.null(common_col)) as.character(taxonomy_df[[common_col]]) else NA_character_,
    rank = if (!is.null(rank_col)) as.character(taxonomy_df[[rank_col]]) else NA_character_,
    stringsAsFactors = FALSE
  )

  out <- out[!is.na(out$tax_id), , drop = FALSE]
  out$label <- paste0(
    ifelse(is.na(out$scientific_name) | out$scientific_name == "", "Unknown organism", out$scientific_name),
    " (", out$tax_id, ")",
    ifelse(!is.na(out$common_name) & out$common_name != "", paste0(" - ", out$common_name), ""),
    ifelse(!is.na(out$rank) & out$rank != "", paste0(" [", out$rank, "]"), "")
  )
  out[!duplicated(out$tax_id), , drop = FALSE]
}

build_orgdb_uniprot_bridge <- function(input_data, key_fun, uniprot_ids,
                                       gene_id_type = NULL, orgdb = NULL) {
  empty <- uniprot_ids[FALSE, , drop = FALSE]
  if (is.null(gene_id_type) || is.null(orgdb) || !nzchar(gene_id_type) || !nzchar(orgdb)) return(empty)
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace(orgdb, quietly = TRUE)) return(empty)

  orgdb_obj <- get(orgdb, envir = asNamespace(orgdb))
  keytype <- toupper(gene_id_type)
  if (!keytype %in% AnnotationDbi::keytypes(orgdb_obj)) return(empty)

  candidate_cols <- intersect(c("UNIPROT", "SYMBOL", "ENTREZID", "ENSEMBL"), AnnotationDbi::columns(orgdb_obj))
  candidate_cols <- setdiff(candidate_cols, keytype)
  if (length(candidate_cols) == 0) return(empty)

  lookup <- data.frame(
    gene_id = input_data$gene_id,
    annotation_key = input_data$annotation_key,
    lookup_id = trimws(as.character(input_data$gene_id)),
    stringsAsFactors = FALSE
  )
  if (keytype == "ENSEMBL") lookup$lookup_id <- sub("\\.[0-9]+$", "", lookup$lookup_id)
  lookup <- lookup[!is.na(lookup$lookup_id) & lookup$lookup_id != "", , drop = FALSE]
  if (nrow(lookup) == 0) return(empty)

  mapped <- tryCatch(
    suppressMessages(AnnotationDbi::select(
      orgdb_obj,
      keys = unique(lookup$lookup_id),
      keytype = keytype,
      columns = candidate_cols
    )),
    error = function(e) NULL
  )
  if (is.null(mapped) || nrow(mapped) == 0 || !keytype %in% names(mapped)) return(empty)

  mapped <- dplyr::inner_join(
    lookup,
    mapped,
    by = stats::setNames(keytype, "lookup_id")
  )
  if (nrow(mapped) == 0) return(empty)

  long <- mapped |>
    dplyr::select(gene_id, input_annotation_key = annotation_key, dplyr::all_of(candidate_cols)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(candidate_cols),
      names_to = "mapped_source",
      values_to = "mapped_id"
    ) |>
    dplyr::filter(!is.na(.data$mapped_id), .data$mapped_id != "") |>
    tidyr::separate_rows(mapped_id, sep = "[;, ]+") |>
    dplyr::mutate(mapped_key = key_fun(.data$mapped_id)) |>
    dplyr::filter(!is.na(.data$mapped_key), .data$mapped_key != "")

  matched <- long |>
    dplyr::inner_join(
      uniprot_ids |>
        dplyr::select(mapped_key = annotation_key, UniProt, UniProt_entry),
      by = "mapped_key"
    ) |>
    dplyr::transmute(
      annotation_key = .data$input_annotation_key,
      UniProt = .data$UniProt,
      UniProt_entry = .data$UniProt_entry,
      id_source = paste0("OrgDb_", keytype, "_to_", .data$mapped_source),
      gene_id_for_merge = .data$gene_id
    ) |>
    dplyr::filter(!is.na(.data$UniProt), .data$UniProt != "") |>
    dplyr::distinct(.data$annotation_key, .data$UniProt, .keep_all = TRUE)

  if (nrow(matched) == 0) return(empty)
  matched
}

build_uniprot_description_file <- function(input_data = NULL,
                                           tax_id = 3702,
                                           reviewed_only = FALSE,
                                           gene_id_source = NULL,
                                           gene_id_type = NULL,
                                           orgdb = NULL,
                                           output_dir = NULL) {
  uniprot_required_packages()

  if (!is.null(input_data)) {
    input_data <- as.data.frame(input_data, check.names = FALSE)
    if (nrow(input_data) == 0) stop("Input data has no rows.")
    names(input_data)[1] <- "gene_id"
  }

  key_fun <- if (exists("gene_join_key", mode = "function")) {
    gene_join_key
  } else {
    function(x) toupper(sub("^.*:", "", sub("\\.[0-9]+$", "", trimws(as.character(x)))))
  }

  if (!is.null(input_data)) {
    input_data$gene_id <- trimws(as.character(input_data$gene_id))
    input_data$annotation_key <- key_fun(input_data$gene_id)
    input_data <- input_data[!is.na(input_data$annotation_key) & input_data$annotation_key != "", , drop = FALSE]
  }

  uniprot_df <- download_uniprot_annotation(tax_id = tax_id, reviewed_only = reviewed_only)
  id_cols <- uniprot_id_columns(uniprot_df)

  all_ids_df <- uniprot_df |>
    dplyr::mutate(
      UniProt = .data$Entry,
      UniProt_entry = .data$Entry.Name
    ) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(id_cols), as.character)) |>
    dplyr::select(UniProt, UniProt_entry, dplyr::all_of(id_cols)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(id_cols),
      names_to = "id_source",
      values_to = "gene_id_for_merge"
    ) |>
    dplyr::filter(!is.na(.data$gene_id_for_merge)) |>
    tidyr::separate_rows(gene_id_for_merge, sep = "[;, ]+") |>
    dplyr::mutate(gene_id_for_merge = stringr::str_trim(.data$gene_id_for_merge)) |>
    dplyr::filter(!is.na(.data$gene_id_for_merge), .data$gene_id_for_merge != "", .data$gene_id_for_merge != "-") |>
    dplyr::mutate(annotation_key = key_fun(.data$gene_id_for_merge)) |>
    dplyr::filter(!is.na(.data$annotation_key), .data$annotation_key != "") |>
    dplyr::distinct(.data$annotation_key, .data$UniProt, .keep_all = TRUE)

  without_prefix <- all_ids_df |>
    dplyr::filter(stringr::str_detect(.data$gene_id_for_merge, ":")) |>
    dplyr::mutate(
      gene_id_for_merge = stringr::str_replace(.data$gene_id_for_merge, "^.*:", ""),
      annotation_key = key_fun(.data$gene_id_for_merge),
      id_source = paste0(.data$id_source, "_without_prefix")
    )

  without_version <- all_ids_df |>
    dplyr::filter(stringr::str_detect(.data$gene_id_for_merge, "\\.[0-9]+$")) |>
    dplyr::mutate(
      gene_id_for_merge = stringr::str_replace(.data$gene_id_for_merge, "\\.[0-9]+$", ""),
      annotation_key = key_fun(.data$gene_id_for_merge),
      id_source = paste0(.data$id_source, "_without_version")
    )

  all_ids_df_one_match <- dplyr::bind_rows(all_ids_df, without_prefix, without_version) |>
    dplyr::filter(!is.na(.data$annotation_key), .data$annotation_key != "") |>
    dplyr::group_by(.data$annotation_key) |>
    dplyr::slice(1) |>
    dplyr::ungroup()

  if (!is.null(input_data)) {
    orgdb_bridge <- build_orgdb_uniprot_bridge(
      input_data = input_data,
      key_fun = key_fun,
      uniprot_ids = all_ids_df_one_match,
      gene_id_type = gene_id_type,
      orgdb = orgdb
    )
    if (nrow(orgdb_bridge) > 0) {
      all_ids_df_one_match <- dplyr::bind_rows(orgdb_bridge, all_ids_df_one_match) |>
        dplyr::filter(!is.na(.data$annotation_key), .data$annotation_key != "") |>
        dplyr::group_by(.data$annotation_key) |>
        dplyr::slice(1) |>
        dplyr::ungroup()
    }
  }

  uniprot_annotation <- uniprot_df |>
    dplyr::transmute(
      UniProt = .data$Entry,
      UniProt_entry = .data$Entry.Name,
      UniProt_gene_names = .data$Gene.Names,
      UniProt_primary_gene = .data$Gene.Names..primary.,
      UniProt_gene_synonyms = .data$Gene.Names..synonym.,
      UniProt_ordered_locus = .data$Gene.Names..ordered.locus.,
      UniProt_ORF = .data$Gene.Names..ORF.,
      GeneID = .data$GeneID,
      RefSeq = .data$RefSeq,
      Ensembl = .data$Ensembl,
      Organism = .data$Organism,
      Organism_ID = .data$Organism..ID.,
      Protein_name = .data$Protein.names,
      Function_description = .data$Function..CC.,
      GO_biological_process = .data$Gene.Ontology..biological.process.,
      GO_molecular_function = .data$Gene.Ontology..molecular.function.,
      GO_cellular_component = .data$Gene.Ontology..cellular.component.,
      KEGG_pathway = .data$Cross.reference..KEGG.
    )

  if (!is.null(gene_id_source) || is.null(input_data)) {
    selected_source <- resolve_uniprot_id_source(gene_id_source, id_cols)
    selected_ids <- all_ids_df |>
      dplyr::filter(.data$id_source == selected_source) |>
      dplyr::transmute(
        gene_id = .data$gene_id_for_merge,
        annotation_key = .data$annotation_key,
        UniProt = .data$UniProt,
        UniProt_entry = .data$UniProt_entry,
        id_source = .data$id_source
      ) |>
      dplyr::filter(!is.na(.data$gene_id), .data$gene_id != "") |>
      dplyr::group_by(.data$annotation_key) |>
      dplyr::slice(1) |>
      dplyr::ungroup()

    if (!is.null(input_data)) {
      keep_keys <- unique(input_data$annotation_key)
      selected_ids <- selected_ids |>
        dplyr::filter(.data$annotation_key %in% keep_keys)
    }

    final_df <- selected_ids |>
      dplyr::left_join(uniprot_annotation, by = c("UniProt", "UniProt_entry")) |>
      dplyr::mutate(
        Symbol = sub(" .*", "", .data$UniProt_gene_names),
        Symbol = ifelse(is.na(.data$Symbol) | toupper(.data$Symbol) == toupper(.data$gene_id), NA, .data$Symbol)
      ) |>
      dplyr::select(-annotation_key) |>
      dplyr::relocate(Symbol, .before = dplyr::any_of("Protein_name")) |>
      dplyr::relocate(
        dplyr::any_of(c("id_source", "UniProt_gene_names", "UniProt_primary_gene", "UniProt_gene_synonyms",
                        "UniProt_ordered_locus", "UniProt_ORF", "GeneID", "RefSeq", "Ensembl",
                        "UniProt", "UniProt_entry", "Organism", "Organism_ID")),
        .after = dplyr::last_col()
      )

    if (!is.null(output_dir)) {
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      safe_source <- gsub("[^A-Za-z0-9_]+", "_", selected_source)
      out_path <- file.path(output_dir, paste0("uniprot_description_taxid_", tax_id, "_", safe_source, "_", Sys.Date(), ".csv"))
      utils::write.csv(final_df, out_path, row.names = FALSE)
      attr(final_df, "path") <- out_path
    }

    return(final_df)
  }

  final_df <- input_data |>
    dplyr::left_join(all_ids_df_one_match, by = "annotation_key") |>
    dplyr::left_join(uniprot_annotation, by = c("UniProt", "UniProt_entry")) |>
    dplyr::mutate(
      Symbol = sub(" .*", "", .data$UniProt_gene_names),
      Symbol = ifelse(is.na(.data$Symbol) | toupper(.data$Symbol) == toupper(.data$gene_id), NA, .data$Symbol)
    ) |>
    dplyr::select(-annotation_key, -dplyr::any_of("gene_id_for_merge")) |>
    dplyr::relocate(Symbol, .before = dplyr::any_of("Protein_name")) |>
    dplyr::relocate(
      dplyr::any_of(c("id_source", "UniProt_gene_names", "UniProt_primary_gene", "UniProt_gene_synonyms",
                      "UniProt_ordered_locus", "UniProt_ORF", "GeneID", "RefSeq", "Ensembl",
                      "UniProt", "UniProt_entry", "Organism", "Organism_ID")),
      .after = dplyr::last_col()
    )

  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    out_path <- file.path(output_dir, paste0("uniprot_description_taxid_", tax_id, "_", Sys.Date(), ".csv"))
    utils::write.csv(final_df, out_path, row.names = FALSE)
    attr(final_df, "path") <- out_path
  }

  final_df
}

build_uniprot_fun <- function(input_file, tax_id = 3702) {
  build_uniprot_description_file(input_file, tax_id = tax_id)
}
