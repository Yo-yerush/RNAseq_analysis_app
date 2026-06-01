# Local RNA-seq dashboard app
# Run with: source('launch_app.R') or double-click RA_RNAseq_analysis_app.bat on Windows.

suppressPackageStartupMessages({
  suppressWarnings({
    library(shiny)
    library(shinyFiles)
    library(DT)
    library(ggplot2)
    library(dplyr)
    library(readr)
    library(stringr)
    if (requireNamespace("shinythemes", quietly = TRUE)) library(shinythemes)
    if (requireNamespace("bslib", quietly = TRUE)) library(bslib)
  })
})

# Allow large annotation uploads, including compressed RefSeq GTF files.
options(shiny.maxRequestSize = 500 * 1024^2)

source(file.path("R", "opening_yo.R"), local = TRUE)
source(file.path("R", "helpers.R"), local = TRUE)
source(file.path("R", "build_uniprot_description_file.R"), local = TRUE)
source(file.path("R", "build_refseq_gftf_description_file.R"), local = TRUE)
source(file.path("legacy_scripts", "volcano_TEG_overlap_with_TE_families_RNAseq.R"), local = TRUE)
source(file.path("legacy_scripts", "genes_into_groups.R"), local = TRUE)
source(file.path("legacy_scripts", "kegg_analysis.R"), local = TRUE)

wrap_html <- function(x, width = 45) {
  HTML(paste(strwrap(x, width = width), collapse = "<br>"))
}

# Keep raw analysis objects and downloads unchanged; only round selected
# DE-result columns when DataTables render in the browser.
datatable <- function(data, ..., numeric_digits = 4) {
  tbl <- DT::datatable(data, ...)
  df <- tryCatch(as.data.frame(data, check.names = FALSE), error = function(e) NULL)
  if (is.null(df) || ncol(df) == 0) return(tbl)

  round_cols <- names(df)[tolower(names(df)) %in% c("log2foldchange", "basemean")]
  round_cols <- round_cols[vapply(df[round_cols], is.numeric, logical(1))]
  if (length(round_cols) > 0) {
    tbl <- DT::formatRound(tbl, columns = round_cols, digits = numeric_digits, mark = "")
  }

  pvalue_cols <- names(df)[tolower(names(df)) %in% c("padj", "pvalue")]
  pvalue_cols <- pvalue_cols[vapply(df[pvalue_cols], is.numeric, logical(1))]
  if (length(pvalue_cols) > 0) {
    tbl <- DT::formatSignif(tbl, columns = pvalue_cols, digits = numeric_digits - 1, mark = "")
  }

  tbl
}

format_row_detail_value <- function(value, column_name = "") {
  if (length(value) == 0 || is.na(value[1])) return("")
  if (is.numeric(value)) {
    if (tolower(column_name) %in% c("padj", "pvalue")) {
      return(format(signif(value[1], 4), scientific = TRUE, trim = TRUE))
    }
    return(format(round(value[1], 4), scientific = FALSE, trim = TRUE, big.mark = ""))
  }
  as.character(value[1])
}

show_row_details_modal <- function(row, title = "Row details") {
  row <- as.data.frame(row, check.names = FALSE)
  if (nrow(row) == 0) return(NULL)
  values <- row[1, , drop = FALSE]
  showModal(modalDialog(
    title = title,
    div(class = "row-detail-modal",
      tags$table(class = "table table-condensed row-detail-table",
        tags$tbody(lapply(names(values), function(col) {
          tags$tr(
            tags$th(col),
            tags$td(format_row_detail_value(values[[col]], col))
          )
        }))
      )
    ),
    size = "l",
    easyClose = TRUE,
    footer = modalButton("Close")
  ))
}

add_row_detail_buttons <- function(d, table_id) {
  d <- as.data.frame(d, check.names = FALSE)
  detail_buttons <- vapply(seq_len(nrow(d)), function(row_i) {
    paste0(
      "<button type=\"button\" class=\"btn btn-link btn-xs row-detail-btn\" style=\"font-size: 22px; padding: 0; line-height: 1;\" ",
      "data-table=\"", htmltools::htmlEscape(table_id, attribute = TRUE), "\" ",
      "data-row=\"", row_i, "\">⋮</button>"
    )
  }, character(1))

  first_col <- names(d)[1] %||% ""
  first_values <- if (ncol(d) > 0) as.character(d[[1]]) else character()
  has_control_col <- ncol(d) > 0 &&
    !is.na(first_col) &&
    trimws(first_col) == "" &&
    any(grepl("gene-counts-btn", first_values, fixed = TRUE))

  if (isTRUE(has_control_col)) {
    d[[1]] <- paste(detail_buttons, first_values)
  } else {
    d <- data.frame(" " = detail_buttons, d, check.names = FALSE)
  }
  d
}

row_detail_button_callback <- JS(
  "table.on('click', 'button.row-detail-btn', function(e) {",
  "  e.stopPropagation();",
  "  Shiny.setInputValue('row_detail_clicked', {",
  "    table: this.getAttribute('data-table'),",
  "    row: parseInt(this.getAttribute('data-row'), 10)",
  "  }, {priority: 'event'});",
  "});"
)

row_detail_button_defs <- function(target) {
  list(list(targets = target, orderable = FALSE, searchable = FALSE, width = "110px", className = "dt-left"))
}

comparison_display_label <- function(x) {
  x <- gsub("_vs_", " vs ", x, fixed = TRUE)
  gsub("_", " ", x, fixed = TRUE)
}

gene_id_type_choices <- c(
  "TAIR", "ENTREZID", "SYMBOL", "ALIAS", "GENENAME",
  "ENSEMBL", "REFSEQ", "ACCNUM", "UNIPROT"
)

gene_id_type_choices_for_orgdb <- function(orgdb = NULL) {
  choices <- gene_id_type_choices
  if (!is.null(orgdb) && nzchar(orgdb) && requireNamespace("AnnotationDbi", quietly = TRUE) && requireNamespace(orgdb, quietly = TRUE)) {
    orgdb_obj <- tryCatch(getExportedValue(orgdb, orgdb), error = function(e) NULL)
    orgdb_keytypes <- if (is.null(orgdb_obj)) character() else tryCatch(AnnotationDbi::keytypes(orgdb_obj), error = function(e) character())
    orgdb_keytypes <- setdiff(orgdb_keytypes, c("GO", "GOALL", "ONTOLOGY", "ONTOLOGYALL", "EVIDENCE", "EVIDENCEALL"))
    choices <- unique(c(choices[choices %in% orgdb_keytypes], orgdb_keytypes, choices))
  }
  stats::setNames(choices, choices)
}

detect_gene_id_type_from_values <- function(gene_ids) {
  ids <- trimws(as.character(gene_ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0) return(NULL)
  ids <- unique(ids)
  ids <- head(ids, 5000)
  ids <- sub("^.*:", "", ids)
  ids_no_version <- sub("\\.[0-9]+$", "", ids)
  ids_upper <- toupper(ids_no_version)

  scores <- c(
    TAIR = mean(grepl("^AT[1-5CM]G[0-9]{5}$", ids_upper)),
    ALIAS = mean(grepl("^B[0-9]{4}$", ids_upper)),
    ENSEMBL = mean(grepl("^ENS[A-Z]*G[0-9]+$", ids_upper)),
    REFSEQ = mean(grepl("^(NM|NR|XM|XR|NP|XP|YP|WP|NC|NG|NT|NW)_[0-9]+$", ids_upper)),
    ENTREZID = mean(grepl("^[0-9]+$", ids_upper)),
    UNIPROT = mean(grepl("^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])(-[0-9]+)?$", ids_upper)),
    SYMBOL = mean(grepl("^[A-Z][A-Z0-9_.-]{1,30}$", ids_upper))
  )
  scores <- scores[is.finite(scores)]
  if (length(scores) == 0) return(NULL)

  # SYMBOL is intentionally a fallback because many structured IDs also look like valid symbols.
  structured <- scores[setdiff(names(scores), "SYMBOL")]
  if (length(structured) && max(structured) >= 0.35) {
    keytype <- names(which.max(structured))
    confidence <- unname(max(structured))
  } else {
    keytype <- "SYMBOL"
    confidence <- unname(scores[["SYMBOL"]])
    if (is.na(confidence) || confidence < 0.7) return(NULL)
  }
  list(keytype = keytype, confidence = confidence, tested = length(ids))
}

te_file_path <- DEFAULT_TAIR_TE_URL
te_super_family_choices <- tryCatch({
  te_tbl <- load_te_file(te_file_path)
  sort(unique(te_tbl$Transposon_Super_Family))
}, error = function(e) character())
default_taxon_choices <- c(
  "Anopheles gambiae (7165)" = "7165",
  "Arabidopsis thaliana (3702)" = "3702",
  "Bos taurus (9913)" = "9913",
  "Caenorhabditis elegans (6239)" = "6239",
  "Canis lupus familiaris (9615)" = "9615",
  "Danio rerio (7955)" = "7955",
  "Drosophila melanogaster (7227)" = "7227",
  "Escherichia coli K-12 MG1655 (511145)" = "511145",
  "Escherichia coli O157:H7 Sakai (386585)" = "386585",
  "Gallus gallus (9031)" = "9031",
  "Homo sapiens (9606)" = "9606",
  "Mus musculus (10090)" = "10090",
  "Pan troglodytes (9598)" = "9598",
  "Plasmodium falciparum 3D7 (36329)" = "36329",
  "Rattus norvegicus (10116)" = "10116",
  "Oryza sativa japonica (39947)" = "39947",
  "Oryza sativa (4530)" = "4530",
  "Sus scrofa (9823)" = "9823",
  "Zea mays (4577)" = "4577",
  "Solanum lycopersicum (4081)" = "4081",
  "Triticum aestivum (4565)" = "4565",
  "Xenopus laevis (8355)" = "8355",
  "Saccharomyces cerevisiae (559292)" = "559292"
)
organism_analysis_config <- data.frame(
  tax_id = c(
    7165, 3702, 9913, 6239, 9615, 7955, 7227, 511145, 386585, 9031,
    9606, 10090, 9598, 36329, 10116, 39947, 4530, 9823, 4577, 4081,
    4565, 8355, 559292
  ),
  label = c(
    "Anopheles gambiae",
    "Arabidopsis thaliana",
    "Bos taurus",
    "Caenorhabditis elegans",
    "Canis lupus familiaris",
    "Danio rerio",
    "Drosophila melanogaster",
    "Escherichia coli K-12 MG1655",
    "Escherichia coli O157:H7 Sakai",
    "Gallus gallus",
    "Homo sapiens",
    "Mus musculus",
    "Pan troglodytes",
    "Plasmodium falciparum 3D7",
    "Rattus norvegicus",
    "Oryza sativa japonica",
    "Oryza sativa",
    "Sus scrofa",
    "Zea mays",
    "Solanum lycopersicum",
    "Triticum aestivum",
    "Xenopus laevis",
    "Saccharomyces cerevisiae"
  ),
  orgdb = c(
    "org.Ag.eg.db",
    "org.At.tair.db",
    "org.Bt.eg.db",
    "org.Ce.eg.db",
    "org.Cf.eg.db",
    "org.Dr.eg.db",
    "org.Dm.eg.db",
    "org.EcK12.eg.db",
    "org.EcSakai.eg.db",
    "org.Gg.eg.db",
    "org.Hs.eg.db",
    "org.Mm.eg.db",
    "org.Pt.eg.db",
    "org.Pf.plasmo.db",
    "org.Rn.eg.db",
    "org.Osativa.eg.db",
    "org.Osativa.eg.db",
    "org.Ss.eg.db",
    "org.Zm.eg.db",
    "org.Slycopersicum.eg.db",
    "org.Ta.eg.db",
    "org.Xl.eg.db",
    "org.Sc.sgd.db"
  ),
  go_keytype = c(
    "ENTREZID", "TAIR", "ENTREZID", "ENTREZID", "ENTREZID", "ENTREZID", "ENTREZID",
    "ALIAS", "ALIAS", "ENTREZID", "ENTREZID", "ENTREZID", "ENTREZID", "ENTREZID",
    "ENTREZID", "ENTREZID", "ENTREZID", "ENTREZID", "ENTREZID", "ENTREZID",
    "ENTREZID", "ENTREZID", "ENTREZID"
  ),
  topgo_id = c(
    "entrez", "entrez", "entrez", "entrez", "entrez", "entrez", "entrez",
    "alias", "alias", "entrez", "entrez", "entrez", "entrez", "entrez",
    "entrez", "entrez", "entrez", "entrez", "entrez", "entrez", "entrez",
    "entrez", "entrez"
  ),
  kegg_species = c(
    "aga", "ath", "bta", "cel", "cfa", "dre", "dme", "eco", "ecs", "gga",
    "hsa", "mmu", "ptr", "pfa", "rno", "osa", "osa", "ssc", "zma", "sly",
    "tae", "xla", "sce"
  ),
  stringsAsFactors = FALSE
)
known_orgdb_catalog <- data.frame(
  package = c(
    "org.Ag.eg.db",
    "org.At.tair.db",
    "org.Bt.eg.db",
    "org.Ce.eg.db",
    "org.Cf.eg.db",
    "org.Dm.eg.db",
    "org.Dr.eg.db",
    "org.EcK12.eg.db",
    "org.EcSakai.eg.db",
    "org.Gg.eg.db",
    "org.Hbacteriophora.eg.db",
    "org.Hs.eg.db",
    "org.Mm.eg.db",
    "org.Mmu.eg.db",
    "org.Mxanthus.db",
    "org.Osativa.eg.db",
    "org.Pf.plasmo.db",
    "org.Pt.eg.db",
    "org.Rn.eg.db",
    "org.Sc.sgd.db",
    "org.Slycopersicum.eg.db",
    "org.Ss.eg.db",
    "org.Ta.eg.db",
    "org.Xl.eg.db",
    "org.Zm.eg.db"
  ),
  organism = c(
    "Anopheles",
    "Arabidopsis thaliana",
    "Bovine",
    "Caenorhabditis elegans",
    "Canine",
    "Drosophila melanogaster",
    "Danio rerio",
    "E. coli K-12",
    "E. coli Sakai",
    "Chicken",
    "Heterorhabditis bacteriophora",
    "Homo sapiens",
    "Mus musculus",
    "Rhesus macaque",
    "Myxococcus xanthus DK 1622",
    "Oryza sativa",
    "Plasmodium falciparum",
    "Pan troglodytes",
    "Rattus norvegicus",
    "Saccharomyces cerevisiae",
    "Solanum lycopersicum",
    "Pig",
    "Triticum aestivum",
    "Xenopus",
    "Zea mays"
  ),
  stringsAsFactors = FALSE
)
installed_orgdb_package_names <- function() {
  tryCatch({
    pkgs <- rownames(utils::installed.packages())
    sort(pkgs[grepl("^org\\..*\\.db$", pkgs)])
  }, error = function(e) character())
}
make_go_orgdb_choices <- function() {
  installed_orgdb_packages <- installed_orgdb_package_names()
  all_orgdb_packages <- sort(unique(c(known_orgdb_catalog$package, organism_analysis_config$orgdb, installed_orgdb_packages)))
  go_orgdb_choice_labels <- vapply(all_orgdb_packages, function(pkg) {
    label <- known_orgdb_catalog$organism[match(pkg, known_orgdb_catalog$package)]
    if (is.na(label) || !nzchar(label)) label <- pkg
    installed_suffix <- if (pkg %in% installed_orgdb_packages) " [installed]" else ""
    paste0(label, " (", pkg, ")", installed_suffix)
  }, character(1))
  stats::setNames(all_orgdb_packages, go_orgdb_choice_labels)
}
go_orgdb_choices <- make_go_orgdb_choices()
msigdb_species_choices <- tryCatch({
  if (!requireNamespace("msigdbr", quietly = TRUE)) character() else {
    sp <- msigdbr::msigdbr_species()
    sort(unique(as.character(sp$species_name)))
  }
}, error = function(e) character())
default_msigdb_species <- if ("Arabidopsis thaliana" %in% msigdb_species_choices) {
  "Arabidopsis thaliana"
} else {
  ""
}
msigdb_species_select_choices <- c("Not available for selected organism" = "", msigdb_species_choices)
pmn_database_choices <- c("Not available for selected organism" = "", pmn_catalog_choices())

ui <- fluidPage(
  theme = shinythemes::shinytheme("united"),
  tags$head(
    tags$style(HTML("\n    .app-title { margin-top: 10px; margin-bottom: 4px; font-weight: 700; }\n    .muted { color: #666; font-size: 0.92em; }\n    .tab-content { padding: 16px; border: 1px solid #ddd; border-top: none; }\n    .download-row .btn { margin-right: 8px; margin-top: 6px; }\n    .annotation-input-row { display: flex; align-items: stretch; }\n    .annotation-input-row > [class*='col-'] { display: flex; }\n    .annotation-input-row .well { flex: 1; width: 100%; }\n    pre { white-space: pre-wrap; }\n    table.dataTable tbody td { padding-top: 1.5px; padding-bottom: 1.5px; line-height: 1.25; }\n    #settings_toggle_btn { margin: 4px 0 10px 0; }\n    body.settings-hidden #settings_sidebar { display: none; }\n    body.settings-hidden #main_content { width: 100%; }\n    .row-detail-modal { max-height: 70vh; overflow-y: auto; }\n    .row-detail-table th { width: 190px; vertical-align: top; white-space: nowrap; }\n    .row-detail-table td { white-space: pre-wrap; word-break: break-word; }\n  ")),
    tags$script(HTML("\n      $(document).on('click', '#settings_toggle_btn', function() {\n        var hidden = !$('body').hasClass('settings-hidden');\n        $('body').toggleClass('settings-hidden', hidden);\n        $(this).text(hidden ? '☰ ►' : '◄ ☰');\n      });\n    "))
  ),

  # Add the theme selector if shinythemes is installed (for easy testing of themes)
  # if (requireNamespace("shinythemes", quietly = TRUE)) shinythemes::themeSelector(),

  titlePanel(div(class = "app-title", "RNA-seq Analysis Dashboard")),
  div(class = "muted", "Load DE results directly, or run DESeq2 from `RSEM .genes.results` files or `featureCounts` output."),
  br(),
  div(
    style = "text-align: left;",
    tags$button(
      id = "settings_toggle_btn",
      type = "button",
      class = "btn btn-warning btn-xs",
      "◄ ☰"
    )
  ),

  sidebarLayout(
    sidebarPanel(id = "settings_sidebar", width = 3,
      conditionalPanel("input.tabs == 'Data input'",
        h4("Data input"),
        radioButtons("data_mode", NULL,
          choices = c(
            "Upload DE results (Excel/CSV/TSV)" = "csv",
            "Run DESeq2 (RSEM/featureCounts)" = "rsem"
          ), selected = "csv"),

        conditionalPanel("input.data_mode == 'csv'",
          fileInput("de_file", "DE results table", accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
          div(class = "muted", "Required columns: gene_id, log2FoldChange, padj. baseMean enables MA plot.")
        ),

        conditionalPanel("input.data_mode == 'rsem'",
          radioButtons("deseq_input_type", "DESeq2 input",
            choices = c("RSEM .genes.results folder" = "rsem", "featureCounts output table" = "featurecounts"),
            selected = "rsem"),
          conditionalPanel("input.deseq_input_type == 'rsem'",
            shinyDirButton("choose_rsem_dir", "Choose RSEM folder", "Select a folder"),
            textInput("rsem_path", "RSEM folder path", value = ""),
            actionButton("scan_rsem", "Scan folder", class = "btn-primary")
          ),
          conditionalPanel("input.deseq_input_type == 'featurecounts'",
            fileInput("featurecounts_file", "featureCounts output", accept = c(".txt", ".tsv", ".csv")),
            actionButton("scan_featurecounts", "Load samples", class = "btn-primary"),
            div(class = "muted", "Uses columns after gene_biotype when present; otherwise columns after Length.")
          ),
          tags$hr(),
          fileInput("coldata_file", "Optional colData (Excel/CSV/TSV)", accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
          div(class = "muted", "The table must include a unique `sample-ID` column and a non-unique `group/condition` column. It may also include `sample labels` and additional `effect/batch` columns. You can also edit the scanned table."),
          actionLink("show_coldata_example", "colData example"),
          tags$hr(),
          uiOutput("contrast_ui"),
          uiOutput("effect_ui"),
          numericInput("min_count", "Filter genes: min total count", value = 10, min = 0, step = 1),
          checkboxInput("lfc_shrink", "Use ashr LFC shrinkage", value = FALSE),
          actionButton("run_deseq", "Run DESeq2", class = "btn-success"),
          tags$hr(),
          actionButton("run_deseq_pca_venn", "PCA & Venn-diagram [all samples]", class = "btn-default"),
          div(class = "muted", "This DE results are for visualization only and may be different from the main separated pairwise comparison."),
          checkboxInput("pca_venn_pairwise", "Run Venn as separate pairwise DE comparisons (can filter more/less genes while analyze. takes longer time)", value = FALSE)
        ),
        tags$hr(),
      ),

      h4("Thresholds"),
      numericInput("alpha", "padj cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
      numericInput("lfc_cutoff", "|log2FC| cutoff", value = 1, min = 0, step = 0.25),
      conditionalPanel("input.tabs == 'GO'",
        tags$hr(),
        numericInput("go_pcut", "GO p-value display cutoff", value = 0.01, min = 0, max = 1, step = 0.001),
        selectInput("ontology", "GO ontology", choices = c("Biological Process" = "BP", "Molecular Function" = "MF", "Cellular Component" = "CC"), selected = "BP"),
        numericInput("top_n", "Top GO terms to display", value = 20, min = 5, max = 200, step = 5)
      ),
      conditionalPanel("input.tabs == 'Hallmark'",
        tags$hr(),
        numericInput("msigdb_pcut", "Hallmark FDR display cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput("msigdb_top_n", "Top Hallmark terms to display", value = 20, min = 5, max = 50, step = 5)
      ),
      conditionalPanel("input.tabs == 'PMN (plants)'",
        tags$hr(),
        numericInput("pmn_pcut", "PMN FDR display cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
        numericInput("pmn_top_n", "Top PMN pathways to display", value = 20, min = 5, max = 100, step = 5)
      ),

      tags$hr(),
      
      h4("Plot Settings"),
      tags$label("Point colors", style = "font-weight: 600; font-size: 0.95em;"),
      tags$div(style = "display: flex; align-items: center; gap: 8px; margin-top: 4px;",
        tags$input(type = "color", id = "color_up", value = "#B2182B",
          style = "width: 30px; height: 22px; border: none; cursor: pointer; padding: 0;",
          oninput = "Shiny.setInputValue('color_up', this.value, {priority: 'event'})"),
        tags$span("Upregulated", style = "font-size: 0.9em;")
      ),
      tags$div(style = "display: flex; align-items: center; gap: 8px; margin-top: 4px;",
        tags$input(type = "color", id = "color_down", value = "#2166AC",
          style = "width: 30px; height: 22px; border: none; cursor: pointer; padding: 0;",
          oninput = "Shiny.setInputValue('color_down', this.value, {priority: 'event'})"),
        tags$span("Downregulated", style = "font-size: 0.9em;")
      ),
      tags$div(style = "display: flex; align-items: center; gap: 8px; margin-top: 4px;",
        tags$input(type = "color", id = "color_ns", value = "#B3B3B3",
          style = "width: 30px; height: 22px; border: none; cursor: pointer; padding: 0;",
          oninput = "Shiny.setInputValue('color_ns', this.value, {priority: 'event'})"),
        tags$span("Not significant", style = "font-size: 0.9em;")
      ),

      tags$hr(),
      
      sliderInput("plot_point_size", "Point size", min = 0.1, max = 5, value = 1, step = 0.1),
      sliderInput("plot_alpha", "Point transparency", min = 0.1, max = 1, value = 0.65, step = 0.05),

      tags$hr(),

      selectInput("plot_theme", "Plot theme",
        choices = c(
          "Classic" = "classic",
          "Black and white" = "bw",
          "Minimal" = "minimal",
          "Linedraw" = "linedraw",
          "Light" = "light",
          "Gray" = "gray",
          "Dark" = "dark",
          "Void" = "void"
        ),
        selected = "classic"
      ),
      selectizeInput("plot_font_family", "Font family",
        choices = c("serif", "sans", "mono"),
        selected = "serif",
        options = list(create = TRUE)
      ),

      tags$hr(),
      
      # Per-tab width/height sliders with sensible defaults per plot type
      conditionalPanel("input.tabs == 'Data input'",
        sliderInput("venn_plot_width",  "Venn plot width (px)",  min = 300, max = 1600, value = 650, step = 50),
        sliderInput("venn_plot_height", "Venn plot height (px)", min = 250, max = 1200, value = 350, step = 50)
      ),
      conditionalPanel("input.tabs == 'DE results'",
        sliderInput("de_plot_width",  "Plot width (px)",  min = 200, max = 1600, value = 350, step = 50),
        sliderInput("de_plot_height", "Plot height (px)", min = 200, max = 1200, value = 200, step = 50)
      ),
      conditionalPanel("input.tabs == 'GO'",
        sliderInput("go_plot_width",  "Plot width (px)",  min = 200, max = 1600, value = 450, step = 50),
        sliderInput("go_plot_height", "Plot height (px)", min = 200, max = 1200, value = 400, step = 50)
      ),
      conditionalPanel("input.tabs == 'Hallmark'",
        sliderInput("msigdb_plot_width",  "Plot width (px)",  min = 200, max = 1600, value = 550, step = 50),
        sliderInput("msigdb_plot_height", "Plot height (px)", min = 200, max = 1200, value = 400, step = 50)
      ),
      conditionalPanel("input.tabs == 'TE analysis (At)'",
        sliderInput("teg_plot_width",  "Plot width (px)",  min = 200, max = 1600, value = 550, step = 50),
        sliderInput("teg_plot_height", "Plot height (px)", min = 100, max = 1200, value = 300, step = 50)
      ),
      conditionalPanel("input.tabs == 'Genes groups'",
        sliderInput("grp_plot_width",  "Plot width (px)",  min = 200, max = 1600, value = 500, step = 50),
        sliderInput("grp_plot_height", "Plot height (px)", min = 200, max = 1200, value = 200, step = 50)
      ),
      conditionalPanel("input.tabs == 'KEGG'",
        sliderInput("kegg_plot_width",  "Plot width (px)",  min = 200, max = 1600, value = 850, step = 50),
        sliderInput("kegg_plot_height", "Plot height (px)", min = 200, max = 1200, value = 350, step = 50)
      ),
      conditionalPanel("input.tabs == 'PMN (plants)'",
        sliderInput("pmn_plot_width",  "Plot width (px)",  min = 200, max = 1600, value = 550, step = 50),
        sliderInput("pmn_plot_height", "Plot height (px)", min = 200, max = 1200, value = 350, step = 50)
      )
    ),


    mainPanel(id = "main_content", width = 9,
      tabsetPanel(id = "tabs",
        tabPanel("Data input",
          conditionalPanel("input.data_mode == 'rsem'",
          div(
            style = "display: flex; align-items: center; justify-content: flex-start; gap: 8px;",
            h4(style = "margin: 0;", "Editable colData"),
            actionButton(
              "add_coldata_effect_col",
              "+",
              class = "btn-default btn-xs",
              style = "padding: 0.5px 4px; font-size: 12px;"
            )
          ),
            DTOutput("coldata_table"),
            div(class = "muted", "Edit the condition column, then select treatment/control in the sidebar and run DESeq2."),
            div(class = "muted", "Use ✚ button to add effect column."),
            tags$hr()
          ),
          uiOutput("data_summary_box"),
          uiOutput("all_comparison_venn_ui"),
          conditionalPanel("input.data_mode == 'rsem'",
            fluidRow(
              column(8,
                h4("PCA"),
                plotOutput("pca_plot", width = "auto", height = "auto"),
                textOutput("pca_message"),
                download_plot_ui("pca", "Download PCA plot"),
                div(class = "download-row", downloadButton("download_pca_table", "Download PCA table"))
              ),
              column(4,
                wellPanel(
                  h4("PCA options"),
                  checkboxInput("pca_show_labels", "Show sample labels", value = TRUE),
                  selectInput("pca_color_palette", "PCA color set",
                    choices = c(
                      "default" = "default",
                      "Okabe-Ito" = "okabe_ito",
                      "Set 1" = "set1",
                      "Set 2" = "set2",
                      "Dark 2" = "dark2",
                      "Paired" = "paired",
                      "Tableau" = "tableau",
                      "Viridis" = "viridis",
                      "Plasma" = "plasma",
                      "Pastel" = "pastel"
                    ),
                    selected = "default"
                  ),
                  uiOutput("pca_conditions_ui")
                )
              )
            )
          ),
          tags$hr(),
          h4("DE table preview"),
          DTOutput("de_preview"),
          div(class = "download-row", downloadButton("download_de", "Download DE table"), downloadButton("download_norm_counts", "Download normalized counts"))
        ),

        tabPanel("Organism annotations",
          fluidRow(
            column(8,
              wellPanel(
                h4("Organism"),
                selectizeInput("annotation_taxon_choice", "Search organism name or taxonomy ID",
                  choices = default_taxon_choices,
                  selected = "3702",
                  options = list(
                    placeholder = "Type organism name or NCBI taxonomy ID",
                    create = TRUE
                  )
                ),
                uiOutput("organism_availability_ui"),
                textOutput("taxonomy_search_status")
              )
            ),
            column(4,
              wellPanel(
                h4("Analysis ID settings"),
                selectInput("go_keytype", "Gene ID type", choices = gene_id_type_choices_for_orgdb("org.At.tair.db"), selected = "TAIR"),
                textOutput("gene_id_type_detection"),
                div(class = "muted", "Used for description-file building, GO, KEGG, Hallmark, and PMN analysis.")
              )
            )
          ),
          tags$hr(),
          fluidRow(class = "annotation-input-row",
            column(4,
              wellPanel(
                h4("Load annotation table"),
                fileInput("annotation_file", "Annotation CSV/TSV/TXT", accept = c(".csv", ".tsv", ".txt")),
                div(class = "muted", "The first column, or a column named gene_id/gene/GeneID, is used to match the loaded DE table. All other columns are kept as annotations."),
                actionButton("load_annotation_file", "Load annotations", class = "btn-primary", style = "width:100%;")
              )
            ),
            column(4,
              wellPanel(
                h4("Build from UniProt"),
                checkboxInput("annotation_reviewed_only", "Reviewed Swiss-Prot entries only", value = FALSE),
                actionButton("scan_uniprot_id_sources", "Scan ID sources", class = "btn-primary", style = "width:100%; margin-bottom: 8px;"),
                uiOutput("uniprot_id_source_ui"),
                textOutput("uniprot_status"),
                div(class = "muted", "Downloads UniProt annotations for the organism selected above, then uses the selected source as the new gene_id column. Requires internet access."),
                actionButton("build_uniprot_annotations", "Build description file", class = "btn-success", style = "width:100%;")
              )
            ),
            column(4,
              wellPanel(
                h4("Build from RefSeq GTF"),
                fileInput("refseq_gtf_file", "NCBI RefSeq GTF", accept = c(".gtf", ".gtf.gz")),
                uiOutput("refseq_gtf_id_source_ui"),
                checkboxInput("refseq_gtf_one_row_per_id", "One annotation row per selected ID", value = TRUE),
                textOutput("refseq_gtf_status"),
                div(class = "muted", "Parses GTF attributes and db_xref IDs, then uses the selected source as the new gene_id column."),
                actionButton("build_refseq_gtf_annotations", "Build description file", class = "btn-success", style = "width:100%;")
              )
            )
          ),
          fluidRow(
            column(12, uiOutput("arabidopsis_default_ui"))
          ),
          tags$hr(),
          fluidRow(
            column(12,
              wellPanel(
                h4("Current annotation source"),
                verbatimTextOutput("annotation_summary"),
                div(class = "download-row",
                  downloadButton("download_annotation_table", "Download annotation table"),
                  downloadButton("download_annotated_de", "Download annotated DE table")
                )
              )
            )
          ),
          fluidRow(
            column(12,
              uiOutput("annotation_preview_ui")
            )
          )
        ),

        tabPanel("DE results",
          fluidRow(
            column(6,
              h4("Volcano plot"),
              plotOutput("volcano_plot", width = "auto", height = "auto"),
              download_plot_ui("volcano", "Download volcano plot")
            ),
            column(6,
              h4("MA plot"),
              plotOutput("ma_plot", width = "auto", height = "auto"),
              textOutput("ma_message"),
              download_plot_ui("ma", "Download MA plot")
            )
          ),
          tags$hr(),
          fluidRow(
            column(12,
              h4("DE table & description"),
              div(class = "muted", "Search across all annotation columns! Use a space for 'AND' (e.g. 'kinase stress'), or use | for 'OR' (e.g. 'kinase|stress'). You can also filter specific columns using the boxes below the headers."),
              uiOutput("gene_counts_click_hint"),
              tags$hr(),
              DTOutput("search_annotations_table"),
              div(class = "download-row", downloadButton("download_search_annotations", "Download Table"))
            )
          )
        ),

        tabPanel("Genes groups",
          tabsetPanel(
            tabPanel("Gene families - enrichment",
              wellPanel(
                fluidRow(
                  column(4, uiOutput("gene_family_enrichment_status_ui")),
                  column(3, selectInput("gene_family_enrichment_direction", "Gene set", choices = c("up", "down", "all"), selected = "up")),
                  column(3, numericInput("gene_family_min_size", "Min genes per family", value = 3, min = 1, step = 1)),
                  column(2, br(), uiOutput("gene_family_enrichment_run_ui"))
                )
              ),
              fluidRow(
                column(12,
                  h4("Top enriched gene families"),
                  plotOutput("gene_family_enrichment_plot", width = "auto", height = "auto"),
                  download_plot_ui("gene_family_enrichment", "Download enrichment plot"),
                  tags$hr(),
                  DTOutput("gene_family_enrichment_table"),
                  div(class = "download-row", downloadButton("download_gene_family_enrichment_table", "Download enrichment table"))
                )
              )
            ),
            tabPanel("Gene families",
              wellPanel(
                fluidRow(
                  column(8, uiOutput("gene_family_selector_ui")),
                  column(4, br(), uiOutput("gene_family_load_ui"))
                ),
                uiOutput("gene_family_status_ui")
              ),
              fluidRow(
                column(12,
                  h4("Gene family volcano plot"),
                  plotOutput("gene_family_volcano", width = "auto", height = "auto"),
                  download_plot_ui("gene_family_volcano", "Download family volcano"),
                  tags$hr(),
                  h4("Genes in selected family"),
                  DTOutput("gene_family_table"),
                  div(class = "download-row", downloadButton("download_gene_family_table", "Download family table"))
                )
              )
            ),
            tabPanel("Custom groups - RA lab (At)",
              wellPanel(
                fluidRow(
                  column(8, uiOutput("gene_group_selector_ui")),
                  column(4, br(), uiOutput("gene_group_run_ui"))
                ),
                div(class = "muted", style = "margin-top: 8px;", "Downloads gene-set lists from GitHub. Requires internet access and GO.db package.")
              ),
              fluidRow(
                column(12,
                  h4("Volcano plot - highlighted group"),
                  plotOutput("gene_group_volcano", width = "auto", height = "auto"),
                  download_plot_ui("gene_group_volcano", "Download group volcano"),
                  tags$hr(),
                  h4("Genes in selected group"),
                  DTOutput("gene_group_table"),
                  div(class = "download-row", downloadButton("download_gene_group_table", "Download group table"))
                )
              )
            )
          )
        ),

        tabPanel("GO",
          tabsetPanel(
            tabPanel("Enrichment",
              wellPanel(
                fluidRow(
                  column(4, selectInput("go_direction", "Gene set", choices = c("up", "down", "all"), selected = "up")),
                  column(8, selectInput("go_orgdb", "GO OrgDb", choices = go_orgdb_choices, selected = "org.At.tair.db"))
                ),
                fluidRow(
                  column(4, selectInput("go_algorithm", "topGO algorithm", choices = c("weight01", "classic", "elim"), selected = "weight01")),
                  column(4, selectInput("go_statistic", "Statistic", choices = c("fisher", "ks"), selected = "fisher")),
                  column(4, br(), actionButton("run_go", "Run GO enrichment", class = "btn-primary", style = "width:100%;"))
                ),
                div(class = "muted", "Choose any listed OrgDb package. Entries marked [installed] are ready to use; other packages must be installed before GO analysis can run. P-value & top-N filters apply on display without re-running.")
              ),
              fluidRow(
                column(12,
                  h4("GO bubble plot"),
                  plotOutput("go_bubble", width = "auto", height = "auto"),
                  download_plot_ui("go_bubble", "Download GO bubble plot"),
                  tags$hr(),
                  DTOutput("go_table"),
                  div(class = "download-row", downloadButton("download_go_table", "Download GO table"))
                )
              )
            ),
            tabPanel("GO genes",
              wellPanel(
                fluidRow(
                  column(7, uiOutput("go_gene_codes_ui")),
                  column(2, br(), actionButton("run_go_gene_lookup", "Find genes", class = "btn-primary", style = "width:100%;")),
                  column(3, div(class = "muted", style = "margin-top: 8px;", "Uses the selected GO OrgDb, Gene ID type, ontology, and current padj/log2FC thresholds."))
                )
              ),
              fluidRow(
                column(12,
                  h4("GO genes volcano plot"),
                  plotOutput("go_genes_volcano", width = "auto", height = "auto"),
                  download_plot_ui("go_genes_volcano", "Download GO genes volcano"),
                  tags$hr(),
                  h4("GO genes"),
                  DTOutput("go_genes_table"),
                  div(class = "download-row", downloadButton("download_go_genes_table", "Download GO genes table"))
                )
              )
            ),
            tabPanel("REVIGO-like reduction",
              wellPanel(
                div(class = "muted", "Runs GO enrichment for the selected gene set if it is not already cached, then reduces those terms."),
                fluidRow(
                  column(3, selectInput("revigo_direction", "Gene set", choices = c("up", "down", "all"), selected = "up")),
                  column(3, numericInput("revigo_top_n", "Max GO terms", value = 80, min = 5, max = 300, step = 5)),
                  column(3, sliderInput("revigo_threshold", "Similarity reduction threshold", min = 0.4, max = 0.95, value = 0.7, step = 0.05)),
                  column(3, selectInput("revigo_algorithm", "Semantic layout", choices = c("UMAP" = "umap", "PCA" = "pca"), selected = "umap"))
                ),
                fluidRow(
                  column(9, div(class = "muted", "These options require a new REVIGO-like run.")),
                  column(3, actionButton("run_revigo", "Run REVIGO-like analysis", class = "btn-primary", style = "width:100%;"))
                )
              ),
              fluidRow(
                column(3, checkboxInput("revigo_show_labels", "Show parent labels", value = TRUE)),
                column(3, numericInput("revigo_parent_label_count", "Parent labels to show", value = 8, min = 0, max = 100, step = 1)),
                column(6, selectInput("revigo_color_palette", "REVIGO color set",
                  choices = c(
                    "Okabe-Ito" = "okabe_ito",
                    "Set 1" = "set1",
                    "Set 2" = "set2",
                    "Dark 2" = "dark2",
                    "Paired" = "paired",
                    "Tableau" = "tableau",
                    "Viridis" = "viridis",
                    "Plasma" = "plasma",
                    "Pastel" = "pastel"
                  ),
                  selected = "set1"
                ))
              ),
              fluidRow(
                column(6,
                  h4("Scatter plot"),
                  plotOutput("revigo_plot", width = "auto", height = "auto"),
                  download_plot_ui("revigo", "Download scatter plot")
                ),
                column(6,
                  h4("Treemap"),
                  plotOutput("revigo_treemap", width = "auto", height = "auto"),
                  download_plot_ui("revigo_treemap", "Download treemap plot")
                ),
                column(12,
                  div(class = "download-row", downloadButton("download_revigo_table", "Download REVIGO table"))
                )
              )
            ),
            tabPanel("GO offspring",
              wellPanel(
                fluidRow(
                  column(8, uiOutput("parent_go_ids_ui")),
                  column(4, br(), actionButton("run_offspring", "Create offspring summary", class = "btn-primary", style = "width:100%;"))
                )
              ),
              fluidRow(
                column(12,
                  h4("GO offspring summary table"),
                  DTOutput("offspring_table"),
                  div(class = "download-row", downloadButton("download_offspring", "Download offspring table"))
                )
              )
            ),
            tabPanel("Abiotic stress (plants)",
              wellPanel(
                fluidRow(
                  column(6, selectInput("stress_dataset", "Stress test gene set", choices = c("all", "up", "down"), selected = "all")),
                  column(6, br(), actionButton("run_stress", "Run abiotic stress enrichment", class = "btn-primary", style = "width:100%;"))
                )
              ),
              fluidRow(
                column(12,
                  h4("Abiotic stress enrichment"),
                  plotOutput("stress_plot", width = "auto", height = "auto"),
                  download_plot_ui("stress", "Download stress plot"),
                  tags$hr(),
                  DTOutput("stress_table"),
                  div(class = "download-row", downloadButton("download_stress_table", "Download stress table"))
                )
              )
            )
          )
        ),

        tabPanel("KEGG",
          tabsetPanel(
            tabPanel("Enrichment Analysis",
              wellPanel(
                fluidRow(
                  column(3, textInput("kegg_species", "KEGG code", value = "ath")),
                  column(4, br(), actionButton("run_kegg", "Run KEGG Enrichment", class = "btn-primary", style="width:100%")),
                  column(5, div(class = "muted", "Uses the organism selected in the Organism annotations tab and the current padj/log2FC thresholds. The KEGG code is auto-filled when available."))
                )
              ),
              fluidRow(
                column(12,
                  h4("KEGG Pathway Enrichment Plot"),
                  plotOutput("kegg_bubble_plot", width = "auto", height = "auto"),
                  download_plot_ui("kegg_bubble", "Download Bubble Plot"),
                  tags$hr(),
                  h4("Enriched Pathways Table"),
                  DTOutput("kegg_enrichment_table"),
                  div(class = "download-row", downloadButton("download_kegg_table", "Download KEGG table"))
                )
              )
            ),
            tabPanel("Pathview Visualization",
              wellPanel(
                fluidRow(
                  column(8, uiOutput("pathview_selector_ui")),
                  column(4, br(), actionButton("run_pathview", "Generate Pathview Map", class = "btn-primary", style="width:100%"))
                ),
                tags$hr(),
                fluidRow(
                  column(4,
                    numericInput("pathview_gene_limit", "Color key limit", value = 2, min = 0.1, step = 0.05)
                  ),
                  column(4,
                    selectInput("pathview_key_pos", "Color key position",
                      choices = c(
                        "Top right" = "topright",
                        "Top left" = "topleft",
                        "Bottom right" = "bottomright",
                        "Bottom left" = "bottomleft"
                      ),
                      selected = "topright"
                    )
                  ),
                  column(4,
                    tags$label("Near zero color"),
                    tags$div(style = "display: flex; gap: 8px;",
                      tags$input(type = "color", id = "pathview_color_mid", value = "#d4d4d4",
                        style = "width: 44px; height: 32px; cursor: pointer; padding: 0;",
                        oninput = "Shiny.setInputValue('pathview_color_mid', this.value, {priority: 'event'})")
                    )
                  )
                ),
                fluidRow(
                  column(4,
                    checkboxInput("pathview_plot_col_key", "Show color key", value = TRUE)
                  ),
                ),
                div(class = "muted", style = "margin-top: 8px;", "Downloads pathway mapping from KEGG and colors genes by their log2FoldChange. Requires internet access.")
              ),
              fluidRow(
                column(12,
                  h4("Pathview Map"),
                  uiOutput("pathview_image_ui"),
                  uiOutput("pathview_table_ui")
                )
              )
            )
          )
        ),

        tabPanel("Hallmark",
          tabsetPanel(
            tabPanel("Enrichment Analysis",
              wellPanel(
                fluidRow(
                  column(3, selectInput("msigdb_direction", "Gene set", choices = c("up", "down", "all"), selected = "up")),
                  column(4, selectizeInput("msigdb_species", "MSigDB species", choices = msigdb_species_select_choices, selected = default_msigdb_species,
                                           options = list(create = FALSE, placeholder = "Select an MSigDB species"))),
                  column(2, numericInput("msigdb_min_set_size", "Min set size", value = 5, min = 1, step = 1)),
                  column(3, br(), uiOutput("msigdb_run_ui"))
                ),
                div(class = "muted", "Uses MSigDB Hallmark gene sets from the msigdbr package and the Gene ID type selected in Organism annotations.")
              ),
              fluidRow(
                column(12,
                  h4("Hallmark enrichment plot"),
                  plotOutput("msigdb_plot", width = "auto", height = "auto"),
                  download_plot_ui("msigdb", "Download Hallmark plot"),
                  tags$hr(),
                  h4("Hallmark enrichment table"),
                  DTOutput("msigdb_table"),
                  div(class = "download-row", downloadButton("download_msigdb_table", "Download Hallmark table"))
                )
              )
            ),
            tabPanel("Hallmark genes",
              wellPanel(
                fluidRow(
                  column(7, uiOutput("hallmark_gene_codes_ui")),
                  column(2, br(), actionButton("run_hallmark_gene_lookup", "Find genes", class = "btn-primary", style = "width:100%;")),
                  column(3, div(class = "muted", style = "margin-top: 8px;", "Uses the selected MSigDB species, Gene ID type, and current padj/log2FC thresholds."))
                )
              ),
              fluidRow(
                column(12,
                  h4("Hallmark genes volcano plot"),
                  plotOutput("hallmark_genes_volcano", width = "auto", height = "auto"),
                  download_plot_ui("hallmark_genes_volcano", "Download Hallmark genes volcano"),
                  tags$hr(),
                  h4("Hallmark genes"),
                  DTOutput("hallmark_genes_table"),
                  div(class = "download-row", downloadButton("download_hallmark_genes_table", "Download Hallmark genes table"))
                )
              )
            )
          )
        ),

        tabPanel("PMN (plants)",
          tabsetPanel(
            tabPanel("Enrichment Analysis",
              wellPanel(
                fluidRow(
                  column(3, selectInput("pmn_direction", "Gene set", choices = c("up", "down", "all"), selected = "up")),
                  column(4, selectizeInput("pmn_cyc_db", "PMN Cyc DB", choices = pmn_database_choices, selected = "AraCyc",
                                         options = list(create = TRUE, placeholder = "Select or type a PMN Cyc database"))),
                  column(2, numericInput("pmn_min_set_size", "Min set size", value = 3, min = 1, step = 1)),
                  column(3, br(), uiOutput("pmn_run_ui"))
                ),
                div(class = "muted", "Uses Plant Metabolic Network Cyc pathway gene sets for plant organisms. The Cyc DB is auto-filled from the organism selected in Organism annotations when a PMN match is known.")
              ),
              fluidRow(
                column(12,
                  h4("PMN pathway enrichment plot"),
                  plotOutput("pmn_plot", width = "auto", height = "auto"),
                  download_plot_ui("pmn", "Download PMN plot"),
                  tags$hr(),
                  h4("PMN pathway enrichment table"),
                  DTOutput("pmn_table"),
                  div(class = "download-row", downloadButton("download_pmn_table", "Download PMN table"))
                )
              )
            ),
            tabPanel("Pathway genes",
              wellPanel(
                fluidRow(
                  column(7, uiOutput("pmn_pathway_codes_ui")),
                  column(2, br(), actionButton("run_pmn_pathway_lookup", "Find genes", class = "btn-primary", style = "width:100%;")),
                  column(3, div(class = "muted", style = "margin-top: 8px;", "Uses the selected PMN Cyc DB above and current padj/log2FC thresholds."))
                )
              ),
              fluidRow(
                column(12,
                  h4("PMN pathway genes volcano plot"),
                  plotOutput("pmn_pathway_volcano", width = "auto", height = "auto"),
                  download_plot_ui("pmn_pathway_volcano", "Download pathway volcano"),
                  tags$hr(),
                  h4("PMN pathway genes"),
                  DTOutput("pmn_pathway_table"),
                  div(class = "download-row", downloadButton("download_pmn_pathway_table", "Download pathway gene table"))
                )
              )
            )
          )
        ),

        tabPanel("TE analysis (At)",
          tabsetPanel(
            tabPanel("Enrichment Analysis",
              wellPanel(
                fluidRow(
                  column(4, numericInput("te_enrich_pvalue", "p-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
                  column(4, br(), actionButton("run_te_enrich", "Run TE Enrichment", class = "btn-primary", style = "width:100%;")),
                  column(4, div(class = "muted", "Uses Wilcoxon rank-sum test to find overrepresented TE superfamilies."))
                )
              ),
              fluidRow(
                column(12,
                  h4("TE Superfamily Enrichment Plot"),
                  plotOutput("te_enrich_bubble", width = "auto", height = "auto"),
                  download_plot_ui("te_enrich_bubble", "Download Enrichment Plot"),
                  tags$hr(),
                  h4("Enriched Superfamilies Table"),
                  DTOutput("te_enrich_table"),
                  div(class = "download-row", downloadButton("download_te_enrich_table", "Download Table"))
                )
              )
            ),
            tabPanel("Volcano Plot",
              wellPanel(
                fluidRow(
                  column(4,
                    selectizeInput("te_super_families", "TE super-families", choices = te_super_family_choices,
                                   selected = head(te_super_family_choices, 3), multiple = TRUE,
                                   options = list(placeholder = "Type to search super-families", create = FALSE))
                  ),
                  column(3, numericInput("te_padj_cutoff", "padj cutoff", value = 0.05, min = 0, max = 1, step = 0.01)),
                  column(3, numericInput("te_lfc_cutoff", "|log2FC| cutoff", value = 1, min = 0, step = 0.25)),
                  column(2, br(), actionButton("run_te_volcano", "Build TE volcano", class = "btn-primary", style = "width:100%;"))
                ),
                div(class = "muted", "Uses the default Arabidopsis description file and TAIR10 TE metadata from GitHub.")
              ),
              fluidRow(
                column(12,
                  h4("TEG volcano plot"),
                  plotOutput("te_volcano_plot", width = "auto", height = "auto"),
                  download_plot_ui("te_volcano", "Download TE volcano plot"),
                  tags$hr(),
                  DTOutput("te_volcano_table"),
                  div(class = "download-row", downloadButton("download_te_volcano_table", "Download TE volcano table"))
                )
              )
            )
          )
        ),

        tabPanel("Log / Help",
          tabsetPanel(
            tabPanel("Log",
              br(),
              h4("Run Log"),
              div(class = "muted", "App steps, user-facing messages, and errors are recorded here. Package console output and warnings are intentionally not shown."),
              verbatimTextOutput("run_log")
            ),
            tabPanel("Help",
              br(),
              h4("Author"),
              div("Yonatan Yerushalmy"),
              div(class = "muted", "Plant's Metabolism and Molecular Genetic laboratory"),
              div(class = "muted", "Rachel Amir's group"),
              tags$hr(),
              h4("Source Code & Repository"),
              div(class = "muted", "You can download, edit, and access the raw files and source code for this app:"),
              tags$a(href = "https://github.com/Yo-yerush/RNAseq_analysis_app/archive/refs/heads/main.zip", "Download (.zip)", target = "_blank"),
              div(class = "muted", "--"),
              div(class = "muted", "GitHub repository:"),
              tags$a(href = "https://github.com/Yo-yerush/RNAseq_analysis_app", "https://github.com/Yo-yerush/RNAseq_analysis_app", target = "_blank"),
              tags$hr(),
              h4("Notes & Usage"),
              tags$ul(
                tags$li(strong("Data input: "), "Upload DE CSV/TSV/TXT files or run DESeq2 directly from RSEM ", code("*.genes.results"), " files. CSV, TSV, and TXT uploads can use comma or tab delimiters. PCA is shown here only after running DESeq2 from RSEM count data."),
                tags$li(strong("DESeq2 from RSEM: "), "Use the editable colData table to set conditions. Optional extra colData columns can be used as an adjusted effect or as a condition:effect interaction. The Data tab prints the exact model formula and contrast."),
                tags$li(strong("DE results: "), "Volcano and MA plots use ", code("gene_id"), ", ", code("log2FoldChange"), ", ", code("padj"), " and optional ", code("baseMean"), ". The annotation search table is shown below the plots."),
                tags$li(strong("Organism annotations: "), "Choose organism and Gene ID type, load a manual annotation table, or build one from UniProt. Human Ensembl IDs can be bridged through the selected OrgDb when available."),
                tags$li(strong("DE preview: "), "The compact DE summary and preview table are shown in the Data input tab with DE and normalized-count downloads."),
                tags$li(strong("GO analysis: "), "Runs topGO enrichment, REVIGO-like semantic reduction, GO offspring summaries, and abiotic-stress GO summaries. Requires a compatible OrgDb package and Gene ID type."),
                tags$li(strong("KEGG analysis: "), "Downloads/caches KEGG pathways by KEGG organism code. If needed, the app maps the selected Gene ID type to KEGG-compatible Entrez IDs through the selected OrgDb. Pathview uses the same mapping to color pathway genes by log2FC."),
                tags$li(strong("PMN analysis: "), "Runs Plant Metabolic Network pathway enrichment for plant Cyc databases such as AraCyc, OryzaCyc, CornCyc, and TomatoCyc. The app auto-selects a PMN database when the selected organism is mapped; otherwise select or type a Cyc database manually."),
                tags$li(strong("MSigDB/Hallmark: "), "Runs Hallmark over-representation analysis with ", code("msigdbr"), ". The run button appears only for species available in MSigDB."),
                tags$li(strong("TE analysis: "), "Uses the bundled Arabidopsis TAIR/Methylome TE annotation files for TE superfamily enrichment and TE volcano plots. Other organisms require compatible TE-level annotations."),
                tags$li(strong("Genes groups: "), "Gene-family analysis is available for Arabidopsis and human only. Arabidopsis uses the RA lab family file; human uses HGNC family tables after mapping the current Gene ID type to ", code("hgnc_id"), ". The Gene families sub-tab builds selected-family volcano plots and tables, and Gene families - enrichment runs top family enrichment dotplots."),
                tags$li(strong("Downloads: "), "No plots or tables are saved automatically. Use the download buttons on each tab.")
              ),
              tags$hr(),
              h4("AI Developer Instructions"),
              div(class = "muted", "Copy/paste these instructions into an AI coding tool when adding features or debugging this app:"),
              tags$ul(
                tags$li(strong("Start here: "), "Read ", code("app_140526/app.R"), ", ", code("app_140526/R/helpers.R"), ", ", code("app_140526/R/build_uniprot_description_file.R"), " and any relevant file in ", code("app_140526/legacy_scripts/"), ". Preserve existing Shiny patterns instead of redesigning the app."),
                tags$li(strong("App structure: "), code("app.R"), " contains UI, server reactives, observers, downloads, and the shared ", code("rv <- reactiveValues(...)"), " state. Put large reusable analysis functions in ", code("R/helpers.R"), " or a focused script under ", code("legacy_scripts/"), " and source it from ", code("app.R"), "."),
                tags$li(strong("Standard data contract: "), "Most features expect a DE table with ", code("gene_id"), ", ", code("log2FoldChange"), ", ", code("padj"), ", optional ", code("pValue"), " and optional ", code("baseMean"), ". Use ", code("read_any_table()"), ", ", code("standardize_de_table()"), ", ", code("normalize_annotation_table()"), " and ", code("merge_with_description()"), " instead of ad hoc parsing."),
                tags$li(strong("Organism and ID mapping: "), "Use the selected ", code("input$go_keytype"), ", ", code("input$go_orgdb"), ", ", code("input$kegg_species"), " and ", code("input$pmn_cyc_db"), " for organism-aware analyses. For KEGG/Pathview, use ", code("map_de_ids_for_kegg()"), " so Ensembl/SYMBOL IDs are converted to KEGG-compatible IDs. For PMN, use the organism locus IDs expected by the selected Cyc database."),
                tags$li(strong("Adding a new analysis tab: "), "Add a ", code("tabPanel()"), " in the main ", code("tabsetPanel(id = 'tabs')"), ", add tab-specific sidebar controls with ", code("conditionalPanel(\"input.tabs == 'Tab name'\", ...)"), ", add result objects to ", code("rv"), ", clear them in ", code("clear_analysis_results()"), ", and add table/plot/download outputs."),
                tags$li(strong("Run workflow pattern: "), "Use ", code("observeEvent(input$run_x, { withProgress(... tryCatch(...)) })"), ". On success, store raw tables and plots in ", code("rv"), " and call ", code("append_log()"), ". On failure, use ", code("showNotification(e$message, type = 'error')"), " and also log the error."),
                tags$li(strong("Plot downloads: "), "For ggplot outputs, expose a reactive plot object and connect it with ", code("download_plot_ui()"), " in the UI and ", code("download_plot_server()"), " in the server. Keep plot width/height controls tab-specific where possible."),
                tags$li(strong("File downloads: "), "Use safe, informative filenames. For organism-specific files use ", code("safe_filename_part()"), " and include organism name/tax ID when relevant."),
                tags$li(strong("Testing after edits: "), "At minimum run ", code("Rscript -e \"parse('app_140526/app.R')\""), " and source any changed helper scripts. For mapping changes, test a small known ID example, such as human ", code("ENSG00000141510"), " mapping to Entrez ", code("7157"), " or UniProt ", code("P04637"), "."),
                tags$li(strong("Do not break existing behavior: "), "Keep Arabidopsis TAIR workflows working while adding human/other organism support. Avoid changing input IDs unless all server references are updated.")
              )
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  volumes <- shinyFiles::getVolumes()()
  shinyFiles::shinyDirChoose(input, "choose_rsem_dir", roots = volumes, session = session)

  rv <- reactiveValues(
    de_base = NULL,
    de = NULL,
    norm_counts = NULL,
    pca = NULL,
    coldata = NULL,
    de_summary = NULL,
    de_design_formula = NULL,
    de_contrast = NULL,
    deseq_all_comparisons = NULL,
    go_trigger = 0,
    offspring = NULL,
    stress = NULL,
    stress_plot = NULL,
    revigo = NULL,
    revigo_direction = NULL,
    revigo_up = NULL,
    revigo_down = NULL,
    go_gene_lookup = NULL,
    go_genes = NULL,
    te_volcano = NULL,
    te_volcano_plot = NULL,
    gene_groups = NULL,
    gene_group_plots = NULL,
    gg_context = NULL,
    gg_cache = list(),
    gene_family_context = NULL,
    gene_family_cache = list(),
    gene_family_enrichment = NULL,
    msigdb_enrichment = NULL,
    msigdb_plot = NULL,
    hallmark_gene_lookup = NULL,
    hallmark_genes = NULL,
    kegg_enrichment = NULL,
    kegg_bubble = NULL,
    pmn_enrichment = NULL,
    pmn_plot = NULL,
    pmn_pathway_lookup = NULL,
    pmn_pathway_genes = NULL,
    pathview_pathway = NULL,
    pathview_pwid = NULL,
    pathview_table = NULL,
    selected_gene_counts_gene = NULL,
    annotation_df = load_description_file(),
    annotation_label = "Default Arabidopsis description file (GitHub)",
    annotation_replace_gene_id = FALSE,
    annotation_preview_ready = FALSE,
    uniprot_id_choices = NULL,
    uniprot_status = "Scan UniProt ID sources for the selected organism before building.",
    refseq_gtf_id_choices = NULL,
    refseq_gtf_status = "Upload a RefSeq GTF file to scan available ID sources.",
    taxonomy_results = NULL,
    selected_tax_id = 3702,
    selected_organism_label = "Arabidopsis thaliana",
    detected_gene_id_type = NULL,
    detected_gene_id_type_confidence = NULL,
    detected_gene_id_type_signature = NULL,
    selected_uniprot_available = TRUE,
    selected_go_available = requireNamespace("org.At.tair.db", quietly = TRUE),
    selected_kegg_available = TRUE,
    selected_pmn_available = TRUE,
    taxonomy_status = "Choose a listed organism, or type any organism name and select the UniProt match.",
    go_cache = list(),
    log = character()
  )

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  ensure_orgdb_installed <- function(pkg) {
    pkg <- pkg %||% ""
    if (!nzchar(pkg)) return(FALSE)
    if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)
    if (!grepl("^org\\..*\\.db$", pkg)) return(FALSE)

    append_log("Installing OrgDb package:", pkg, level = "STEP")
    showNotification(paste("Installing", pkg, "from Bioconductor. This can take a few minutes."), type = "message", duration = 10)
    tryCatch({
      old_repos <- getOption("repos")
      on.exit(options(repos = old_repos), add = TRUE)
      if (is.null(old_repos[["CRAN"]]) || identical(unname(old_repos[["CRAN"]]), "@CRAN@")) {
        options(repos = c(CRAN = "https://cloud.r-project.org"))
      }
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
      ok <- requireNamespace(pkg, quietly = TRUE)
      if (!isTRUE(ok)) stop("Installation finished, but the package is still not available.")
      updateSelectInput(session, "go_orgdb", choices = make_go_orgdb_choices(), selected = pkg)
      append_log("OrgDb package installed:", pkg)
      showNotification(paste("Installed", pkg), type = "message", duration = 8)
      TRUE
    }, error = function(e) {
      append_log("OrgDb install error for", pkg, ":", e$message)
      showNotification(paste("Could not install", pkg, ":", e$message), type = "error", duration = 15)
      FALSE
    })
  }

  organism_config_for_tax <- function(tax_id) {
    tax_id <- suppressWarnings(as.integer(tax_id))
    hit <- organism_analysis_config[organism_analysis_config$tax_id == tax_id, , drop = FALSE]
    if (nrow(hit) == 0) return(NULL)
    hit[1, , drop = FALSE]
  }

  taxonomy_label_for_tax <- function(tax_id) {
    tax_id <- suppressWarnings(as.integer(tax_id))
    if (is.na(tax_id)) return(NULL)
    if (!is.null(rv$taxonomy_results) && nrow(rv$taxonomy_results) > 0) {
      hit <- rv$taxonomy_results[rv$taxonomy_results$tax_id == tax_id, , drop = FALSE]
      if (nrow(hit) > 0) return(hit$scientific_name[1] %||% hit$label[1])
    }
    cfg <- organism_config_for_tax(tax_id)
    if (!is.null(cfg)) return(cfg$label[1])
    NULL
  }

  apply_selected_taxon_to_analysis <- function(tax_id, organism_name = NULL) {
    tax_id <- suppressWarnings(as.integer(tax_id))
    if (is.na(tax_id)) return()
    rv$selected_tax_id <- tax_id
    rv$selected_organism_label <- organism_name %||% taxonomy_label_for_tax(tax_id) %||% paste0("tax_id ", tax_id)
    rv$uniprot_id_choices <- NULL
    rv$uniprot_status <- "Scan UniProt ID sources for the selected organism before building."
    if (rv$selected_organism_label %in% msigdb_species_choices) {
      updateSelectizeInput(session, "msigdb_species", selected = rv$selected_organism_label)
    } else {
      updateSelectizeInput(session, "msigdb_species", selected = "")
    }
    pmn_db <- pmn_database_for_tax(tax_id, rv$selected_organism_label)
    rv$selected_pmn_available <- nzchar(pmn_db)
    updateSelectizeInput(session, "pmn_cyc_db", selected = pmn_db)
    rv$selected_uniprot_available <- TRUE
    cfg <- organism_config_for_tax(tax_id)
    if (!is.null(cfg)) {
      updateSelectInput(session, "go_orgdb", selected = cfg$orgdb)
      selected_keytype <- rv$detected_gene_id_type %||% cfg$go_keytype
      keytype_choices <- gene_id_type_choices_for_orgdb(cfg$orgdb)
      if (!selected_keytype %in% unname(keytype_choices)) selected_keytype <- cfg$go_keytype
      if (!selected_keytype %in% unname(keytype_choices)) selected_keytype <- unname(keytype_choices)[1]
      updateSelectInput(session, "go_keytype", choices = keytype_choices, selected = selected_keytype)
      updateTextInput(session, "kegg_species", value = cfg$kegg_species)
      rv$selected_go_available <- requireNamespace(cfg$orgdb, quietly = TRUE)
      rv$selected_kegg_available <- !is.na(cfg$kegg_species) && nzchar(cfg$kegg_species)
    } else {
      rv$selected_go_available <- FALSE
      kegg_code <- if (!is.null(organism_name)) find_kegg_species_code(organism_name) else NA_character_
      if (!is.na(kegg_code) && nzchar(kegg_code)) {
        updateTextInput(session, "kegg_species", value = kegg_code)
        rv$selected_kegg_available <- TRUE
      } else {
        rv$selected_kegg_available <- FALSE
      }
    }
    rv$go_cache <- list()
    rv$kegg_enrichment <- NULL
    rv$kegg_bubble <- NULL
    rv$msigdb_enrichment <- NULL
    rv$msigdb_plot <- NULL
    rv$pmn_enrichment <- NULL
    rv$pmn_plot <- NULL
    rv$pmn_pathway_lookup <- NULL
    rv$pmn_pathway_genes <- NULL
    rv$pathview_pathway <- NULL
    rv$pathview_pwid <- NULL
    rv$pathview_table <- NULL
    rv$gene_family_context <- NULL
    rv$gene_family_cache <- list()
    rv$gene_family_enrichment <- NULL
  }

  topgo_id_from_keytype <- function(keytype) {
    keytype <- toupper(keytype %||% "ENTREZID")
    switch(keytype,
      ENTREZID = "entrez",
      SYMBOL = "symbol",
      ALIAS = "alias",
      GENENAME = "genename",
      ENSEMBL = "ensembl",
      REFSEQ = "refseq",
      TAIR = "entrez",
      "entrez"
    )
  }

  auto_update_gene_id_type <- function(de_df, source_label = "loaded table") {
    if (is.null(de_df) || !"gene_id" %in% names(de_df)) return(invisible(NULL))
    ids <- head(as.character(de_df$gene_id), 50)
    signature <- paste(nrow(de_df), paste(ids, collapse = "|"), sep = ":")
    if (identical(signature, rv$detected_gene_id_type_signature)) return(invisible(NULL))
    rv$detected_gene_id_type_signature <- signature

    detected <- detect_gene_id_type_from_values(de_df$gene_id)
    if (is.null(detected) || is.null(detected$keytype) || !detected$keytype %in% gene_id_type_choices) {
      rv$detected_gene_id_type <- NULL
      rv$detected_gene_id_type_confidence <- NULL
      append_log("Could not auto-detect Gene ID type from", source_label, "- keep or choose it manually.")
      return(invisible(NULL))
    }

    rv$detected_gene_id_type <- detected$keytype
    rv$detected_gene_id_type_confidence <- detected$confidence
    updateSelectInput(session, "go_keytype", selected = detected$keytype)
    rv$go_cache <- list()
    append_log(
      "Auto-detected Gene ID type:",
      detected$keytype
    )
    invisible(detected)
  }

  append_log <- function(..., level = NULL) {
    timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
    msg <- paste(..., collapse = " ")
    if (grepl("^warning\\b|\\bwarning:", msg, ignore.case = TRUE)) return(invisible(NULL))
    if (is.null(level)) {
      level <- if (grepl("error|failed|could not|not found|no .*found", msg, ignore.case = TRUE)) {
        "ERROR"
      } else if (grepl("^step:", msg, ignore.case = TRUE)) {
        "STEP"
      } else {
        "MESSAGE"
      }
    }
    msg <- sub("^step:\\s*", "", msg, ignore.case = TRUE)

    if (level == "ERROR") {
      rv$log <- c(isolate(rv$log), paste(timestamp, paste0("[", toupper(level), "]"), msg))   
    } else {
      rv$log <- c(isolate(rv$log), paste(timestamp, msg))
    }
    
    # Show a UI notification for successful actions so the user sees what happened
    # We skip error messages because they already have dedicated showNotification calls
    if (identical(toupper(level), "MESSAGE") && !grepl("no rsem gene.results files found", tolower(msg))) {
      showNotification(msg, type = "message", duration = 4)
    }
    invisible(NULL)
  }

  safe_filename_part <- function(x, fallback = "organism") {
    x <- trimws(as.character(x %||% fallback))
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    if (!nzchar(x)) fallback else x
  }

  msigdb_species_available <- function(species) {
    !is.null(species) && length(species) == 1 && nzchar(species) && species %in% msigdb_species_choices
  }

  default_annotation_label_for_tax <- function(tax_id) {
    tax_id <- suppressWarnings(as.integer(tax_id))
    if (identical(tax_id, 3702L)) return("Default Arabidopsis description file")
    if (identical(tax_id, 9606L)) return("Default human description file")
    if (identical(tax_id, 511145L)) return("Default E. coli K-12 MG1655 description file")
    "Default organism description file (GitHub)"
  }

  clear_analysis_results <- function() {
    rv$go_trigger <- 0
    rv$offspring <- NULL
    rv$stress <- NULL
    rv$stress_plot <- NULL
    rv$revigo <- NULL
    rv$revigo_direction <- NULL
    rv$revigo_up <- NULL
    rv$revigo_down <- NULL
    rv$go_gene_lookup <- NULL
    rv$go_genes <- NULL
    rv$te_volcano <- NULL
    rv$te_volcano_plot <- NULL
    rv$gene_groups <- NULL
    rv$gene_group_plots <- NULL
    rv$gg_context <- NULL
    rv$gg_cache <- list()
    rv$gene_family_context <- NULL
    rv$gene_family_cache <- list()
    rv$gene_family_enrichment <- NULL
    rv$msigdb_enrichment <- NULL
    rv$msigdb_plot <- NULL
    rv$hallmark_gene_lookup <- NULL
    rv$hallmark_genes <- NULL
    rv$kegg_enrichment <- NULL
    rv$kegg_bubble <- NULL
    rv$pmn_enrichment <- NULL
    rv$pmn_plot <- NULL
    rv$pmn_pathway_lookup <- NULL
    rv$pmn_pathway_genes <- NULL
    rv$pathview_pathway <- NULL
    rv$pathview_pwid <- NULL
    rv$pathview_table <- NULL
    rv$go_cache <- list()
  }

  annotation_match_count_for <- function(base_de) {
    if (is.null(base_de) || is.null(rv$annotation_df)) return(NA_integer_)
    if (isTRUE(rv$annotation_replace_gene_id)) {
      lookup <- make_annotation_lookup(rv$annotation_df)
      return(sum(gene_join_key(base_de$gene_id) %in% lookup$lookup_key, na.rm = TRUE))
    }
    sum(gene_join_key(base_de$gene_id) %in% gene_join_key(rv$annotation_df$gene_id), na.rm = TRUE)
  }

  apply_annotation_to_base <- function(base_de, reset_results = TRUE) {
    req(base_de)
    rv$de <- merge_with_description(
      base_de,
      rv$annotation_df,
      replace_gene_id = isTRUE(rv$annotation_replace_gene_id)
    )
    match_count <- annotation_match_count_for(base_de)
    if (isTRUE(rv$annotation_replace_gene_id)) {
      replaced_count <- 0L
      if ("original_gene_id" %in% names(rv$de)) {
        replaced_count <- sum(
          !is.na(rv$de$original_gene_id) &
            gene_join_key(rv$de$gene_id) != gene_join_key(rv$de$original_gene_id),
          na.rm = TRUE
        )
      }
      append_log(
        "Applied annotations:",
        ifelse(is.na(match_count), 0, match_count),
        "matched genes;",
        replaced_count,
        "gene_id values replaced."
      )
    } else if (!is.null(rv$annotation_df)) {
      append_log("Applied annotations:", ifelse(is.na(match_count), 0, match_count), "matched genes.")
    }
    if (isTRUE(reset_results)) clear_analysis_results()
  }

  annotation_match_count <- function() {
    annotation_match_count_for(rv$de_base)
  }

  apply_current_annotation <- function(reset_results = TRUE) {
    req(rv$de_base)
    apply_annotation_to_base(rv$de_base, reset_results = reset_results)
  }

  go_cache_key <- function(direction) {
    paste(direction, input$ontology, input$alpha, input$lfc_cutoff,
          input$go_algorithm, input$go_statistic, input$go_orgdb, input$go_keytype, sep = "|")
  }
  get_cached_go <- function(direction) {
    rv$go_cache[[go_cache_key(direction)]]
  }
  set_cached_go <- function(direction, value) {
    rv$go_cache[[go_cache_key(direction)]] <- value
    value
  }
  compute_go_cached <- function(direction) {
    cached <- get_cached_go(direction)
    if (!is.null(cached)) return(cached)
    # Always compute all terms (p_cutoff = 1); filtering happens at display time
    g <- run_topgo_enrichment(rv$de, direction = direction, ontology = input$ontology,
                              alpha = input$alpha, lfc_cutoff = input$lfc_cutoff,
                              p_cutoff = 1, algorithm = input$go_algorithm,
                              statistic = input$go_statistic,
                              orgdb = input$go_orgdb %||% "org.At.tair.db",
                              topgo_id = topgo_id_from_keytype(input$go_keytype),
                              keytype = input$go_keytype %||% "TAIR")
    set_cached_go(direction, g)
  }

  # Reactive: filtered GO data for current direction — updates without re-running
  go_display_data <- reactive({
    req(rv$go_trigger > 0)
    cached <- get_cached_go(input$go_direction)
    if (is.null(cached)) return(NULL)
    cached[!is.na(cached$pValue_num) & cached$pValue_num <= input$go_pcut, , drop = FALSE]
  })
  # Reactive: GO bubble plot built from filtered data
  go_display_plot <- reactive({
    d <- go_display_data()
    req(!is.null(d) && nrow(d) > 0)
    make_go_bubble_plot(d, paste("GO enrichment:", input$go_direction, input$ontology),
                        input$top_n, direction = input$go_direction, point_alpha = input$plot_alpha,
                        color_up = input$color_up %||% "#B2182B",
                        color_down = input$color_down %||% "#2166AC",
                        color_all = input$color_ns %||% "#B3B3B3",
                        plot_theme = input$plot_theme %||% "classic",
                        font_family = input$plot_font_family %||% "serif")
  })

  output$parent_go_ids_ui <- renderUI({
    d <- tryCatch(go_display_data(), error = function(e) NULL)
    if ((is.null(d) || nrow(d) == 0) && rv$go_trigger > 0) {
      d <- tryCatch(get_cached_go(input$go_direction), error = function(e) NULL)
    }
    ids <- if (!is.null(d) && nrow(d) > 0 && "GO.ID" %in% names(d)) {
      head(stats::na.omit(as.character(d$GO.ID)), 3)
    } else {
      "GO:0006950"
    }
    textAreaInput(
      "parent_go_ids",
      "Parent GO IDs for offspring summary",
      rows = 2,
      value = paste(ids, collapse = ", "),
      placeholder = "Example: GO:0006950. Separate multiple GO IDs with comma, semicolon, or space."
    )
  })

  observeEvent(input$choose_rsem_dir, {
    path <- shinyFiles::parseDirPath(volumes, input$choose_rsem_dir)
    if (length(path) && nzchar(path)) updateTextInput(session, "rsem_path", value = path)
  })

  observeEvent(input$deseq_input_type, {
    rv$coldata <- NULL
  }, ignoreInit = TRUE)

  observeEvent(input$scan_rsem, {
    req(input$rsem_path)
    append_log("Scanning RSEM folder:", input$rsem_path, level = "STEP")
    tbl <- scan_rsem_files(input$rsem_path)
    if (nrow(tbl) == 0) {
      showNotification("No .genes.results files found in this folder", type = "error")
      append_log("No RSEM gene.results files found in", input$rsem_path)
      return()
    }
    rv$coldata <- tbl[, c("sample_id", "condition", "sample_label")]
    append_log("Scanned", nrow(tbl), "RSEM gene.results files.")
  })

  observeEvent(input$scan_featurecounts, {
    req(input$featurecounts_file$datapath)
    append_log("Loading featureCounts samples:", input$featurecounts_file$name, level = "STEP")
    tryCatch({
      tbl <- scan_featurecounts_file(input$featurecounts_file$datapath)
      if (nrow(tbl) == 0) {
        showNotification("No sample count columns found in this featureCounts file", type = "error")
        append_log("No featureCounts sample columns found in", input$featurecounts_file$name)
        return()
      }
      rv$coldata <- tbl[, c("sample_id", "condition", "sample_label")]
      append_log("Loaded", nrow(tbl), "featureCounts samples.")
    }, error = function(e) {
      showNotification(paste("featureCounts load error:", e$message), type = "error", duration = 12)
      append_log("featureCounts load error:", e$message)
    })
  })

  observeEvent(input$coldata_file, {
    req(input$coldata_file$datapath)
    append_log("Loading colData file:", input$coldata_file$name, level = "STEP")
    cd <- read_any_table(input$coldata_file$datapath, source_name = input$coldata_file$name)
    cd <- normalize_coldata(cd)
    rv$coldata <- cd
    append_log("Loaded colData with", nrow(cd), "samples.")
  })

  observeEvent(input$show_coldata_example, {
    showModal(modalDialog(
      title = "Example colData table",
      tags$p("Use sample IDs that match the count-file sample names. Extra columns can be used as an adjusted effect or interaction."),
      tags$table(class = "table table-striped table-condensed",
        tags$thead(
          tags$tr(
            tags$th("sample_id"),
            tags$th("condition"),
            tags$th("batch")
          )
        ),
        tags$tbody(
          tags$tr(tags$td("WT_1"), tags$td("WT"), tags$td("B1")),
          tags$tr(tags$td("WT_2"), tags$td("WT"), tags$td("B2")),
          tags$tr(tags$td("Treatment_1"), tags$td("Treatment"), tags$td("B1")),
          tags$tr(tags$td("Treatment_2"), tags$td("Treatment"), tags$td("B2"))
        )
      ),
      tags$p(class = "muted", "Legacy column names are also accepted, for example x/sample/exp."),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  output$coldata_table <- renderDT({
    req(rv$coldata)
    datatable(rv$coldata, editable = TRUE, rownames = FALSE, options = list(pageLength = 12, scrollX = TRUE))
  })

  observeEvent(input$coldata_table_cell_edit, {
    info <- input$coldata_table_cell_edit
    rv$coldata[info$row, info$col + 1] <- DT::coerceValue(info$value, rv$coldata[info$row, info$col + 1])
  })

  observeEvent(input$add_coldata_effect_col, {
    req(rv$coldata)
    new_col <- make.unique(c(names(rv$coldata), "effect"), sep = "_")
    new_col <- new_col[length(new_col)]
    rv$coldata[[new_col]] <- ""
    append_log("Added editable colData column:", new_col)
  })

  output$contrast_ui <- renderUI({
    cd <- rv$coldata
    if (is.null(cd) || !"condition" %in% names(cd)) return(div(class = "muted", "Scan folder or load colData first."))
    conditions <- unique(as.character(cd$condition))
    conditions <- conditions[!is.na(conditions) & nzchar(conditions)]
    tagList(
      selectInput("treatment", "Treatment", choices = conditions, selected = conditions[min(2, length(conditions))]),
      selectInput("control", "Control", choices = conditions, selected = conditions[1])
    )
  })

  output$effect_ui <- renderUI({
    cd <- rv$coldata
    if (is.null(cd)) return(NULL)
    effect_cols <- setdiff(names(cd), c("sample_id", "file", "de_effect", "de_group"))
    effect_cols <- effect_cols[vapply(effect_cols, function(col) {
      vals <- unique(as.character(cd[[col]]))
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (identical(col, "condition")) {
        length(vals) > 2 && length(vals) < nrow(cd)
      } else {
        length(vals) > 1 && length(vals) < nrow(cd)
      }
    }, logical(1))]
    if (length(effect_cols) == 0) return(NULL)

    selected_col <- input$effect_col
    if (is.null(selected_col) || !selected_col %in% effect_cols) selected_col <- ""
    controls <- list(
      selectInput("effect_col", "Additional effect column", choices = c("None" = "", effect_cols), selected = selected_col)
    )
    if (nzchar(selected_col)) {
      effect_levels <- unique(as.character(cd[[selected_col]]))
      effect_levels <- effect_levels[!is.na(effect_levels) & nzchar(effect_levels)]
      selected_level <- input$effect_level
      if (is.null(selected_level) || !selected_level %in% effect_levels) selected_level <- effect_levels[1]
      controls <- c(controls, list(
        selectInput("effect_level", "Effect level/reference", choices = effect_levels, selected = selected_level)
      ))
      if (identical(selected_col, "condition")) {
        controls <- c(controls, list(
          div(class = "muted", "The condition column is already used in the DESeq2 design. Choose treatment/control above. Interaction requires a separate effect column."),
          tags$hr()
        ))
      } else {
        controls <- c(controls, list(
          checkboxInput("use_interaction", "Use `condition & effect` interaction", value = FALSE),
          div(class = "muted", " "),
          tags$hr()
        ))
      }
    }
    do.call(tagList, controls)
  })

  loaded_de <- reactive({
    if (input$data_mode == "csv") {
      req(input$de_file$datapath)
      standardize_de_table(read_any_table(input$de_file$datapath, source_name = input$de_file$name), merge_default_description = FALSE)
    } else {
      rv$de_base
    }
  })

  comparison_coldata <- function(coldata, treatment, control) {
    if (identical(treatment, control)) stop("Treatment and control must be different.")
    cd <- normalize_coldata(coldata)
    cond <- as.character(cd$condition)
    keep <- !is.na(cond) & cond %in% c(treatment, control)
    out <- cd[keep, , drop = FALSE]
    if (!any(as.character(out$condition) == treatment)) stop("Treatment is not present in colData condition column.")
    if (!any(as.character(out$condition) == control)) stop("Control is not present in colData condition column.")
    out
  }

  run_deseq_for_current_input <- function(coldata, treatment, control, all_vs_control = FALSE) {
    deseq_input_type <- input$deseq_input_type %||% "rsem"
    if (identical(deseq_input_type, "featurecounts")) {
      run_deseq2_from_featurecounts(
        counts_file = input$featurecounts_file$datapath,
        coldata = coldata,
        treatment = treatment,
        control = control,
        lfc_shrink = isTRUE(input$lfc_shrink),
        min_count = input$min_count,
        effect_col = if (!is.null(input$effect_col)) input$effect_col else "",
        effect_level = if (!is.null(input$effect_level)) input$effect_level else "",
        use_interaction = isTRUE(input$use_interaction) && !identical(input$effect_col, "condition"),
        all_vs_control = all_vs_control
      )
    } else {
      run_deseq2_from_rsem(
        folder = input$rsem_path,
        coldata = coldata,
        treatment = treatment,
        control = control,
        lfc_shrink = isTRUE(input$lfc_shrink),
        min_count = input$min_count,
        effect_col = if (!is.null(input$effect_col)) input$effect_col else "",
        effect_level = if (!is.null(input$effect_level)) input$effect_level else "",
        use_interaction = isTRUE(input$use_interaction) && !identical(input$effect_col, "condition"),
        all_vs_control = all_vs_control
      )
    }
  }

  observe({
    if (identical(input$data_mode, "csv")) {
      base_de <- tryCatch(loaded_de(), error = function(e) {
        showNotification(e$message, type = "error", duration = 8)
        NULL
      })
      rv$de_base <- base_de
      if (!is.null(base_de)) {
        auto_update_gene_id_type(base_de, input$data_mode)
        apply_annotation_to_base(base_de, reset_results = FALSE)
      } else {
        rv$de <- NULL
      }
      rv$norm_counts <- NULL
      rv$pca <- NULL
      rv$de_summary <- NULL
      rv$de_design_formula <- NULL
      rv$de_contrast <- NULL
      rv$deseq_all_comparisons <- NULL
      clear_analysis_results()
    }
  })

  observeEvent(input$run_deseq, {
    req(rv$coldata, input$treatment, input$control)
    deseq_input_type <- input$deseq_input_type %||% "rsem"
    if (identical(deseq_input_type, "featurecounts")) {
      req(input$featurecounts_file$datapath)
    } else {
      req(input$rsem_path)
    }
    append_log("Running DESeq2", paste0("(", deseq_input_type, "):"), input$treatment, "vs", input$control, level = "STEP")
    withProgress(message = "Running DESeq2", value = 0.1, {
      tryCatch({
        incProgress(0.2, detail = "Subsetting to comparison samples")
        de_coldata <- comparison_coldata(rv$coldata, input$treatment, input$control)
        incProgress(0.2, detail = if (identical(deseq_input_type, "featurecounts")) "Importing featureCounts matrix" else "Importing RSEM files")
        res <- run_deseq_for_current_input(
          coldata = de_coldata,
          treatment = input$treatment,
          control = input$control,
          all_vs_control = FALSE
        )
        incProgress(0.6, detail = "Preparing DE tables")
        rv$de_base <- res$de_table
        auto_update_gene_id_type(res$de_table, "DESeq2 results")
        apply_annotation_to_base(res$de_table, reset_results = FALSE)
        rv$norm_counts <- res$norm_counts
        rv$pca <- res$pca_table
        rv$de_summary <- res$summary
        rv$de_design_formula <- res$design_formula
        rv$de_contrast <- res$contrast
        rv$deseq_all_comparisons <- NULL
        clear_analysis_results()
        append_log("DESeq2 finished on comparison samples only:", input$treatment, "vs", input$control, "with", nrow(rv$de), "tested genes.")
        # updateTabsetPanel(session, "tabs", selected = "DE results")
      }, error = function(e) {
        showNotification(e$message, type = "error", duration = 12)
        append_log("DESeq2 error:", e$message)
      })
    })
  })

  observeEvent(input$run_deseq_pca_venn, {
    req(rv$coldata, input$treatment, input$control)
    deseq_input_type <- input$deseq_input_type %||% "rsem"
    if (identical(deseq_input_type, "featurecounts")) {
      req(input$featurecounts_file$datapath)
    } else {
      req(input$rsem_path)
    }
    append_log("Running all-sample PCA/Venn vs", input$control, level = "STEP")
    withProgress(message = "Running PCA/Venn", value = 0.1, {
      tryCatch({
        all_coldata <- normalize_coldata(rv$coldata)
        all_conditions <- unique(as.character(all_coldata$condition))
        all_conditions <- all_conditions[!is.na(all_conditions) & nzchar(all_conditions)]
        comparison_levels <- setdiff(all_conditions, input$control)
        if (length(comparison_levels) == 0) stop("No comparison conditions are available for Venn.")

        incProgress(0.2, detail = "Running DESeq2 on all samples")
        all_res <- run_deseq_for_current_input(
          coldata = all_coldata,
          treatment = input$treatment,
          control = input$control,
          all_vs_control = !isTRUE(input$pca_venn_pairwise)
        )
        rv$pca <- all_res$pca_table

        if (isTRUE(input$pca_venn_pairwise)) {
          comparison_tables <- list()
          for (i in seq_along(comparison_levels)) {
            level <- comparison_levels[[i]]
            incProgress(0.7 / length(comparison_levels), detail = paste("Running", level, "vs", input$control))
            pair_coldata <- comparison_coldata(all_coldata, level, input$control)
            pair_res <- run_deseq_for_current_input(
              coldata = pair_coldata,
              treatment = level,
              control = input$control,
              all_vs_control = FALSE
            )
            comparison_tables[[paste0(level, "_vs_", input$control)]] <- pair_res$de_table
          }
          rv$deseq_all_comparisons <- list(
            control = input$control,
            comparisons = comparison_levels,
            tables = comparison_tables
          )
          append_log("PCA/Venn finished with", length(comparison_tables), "separate pairwise comparisons. Main DE result was not changed.")
        } else {
          rv$deseq_all_comparisons <- all_res$all_comparisons
          append_log("PCA/Venn finished with one all-sample DESeq2 run. Main DE result was not changed.")
        }
      }, error = function(e) {
        showNotification(e$message, type = "error", duration = 12)
        append_log("PCA/Venn error:", e$message)
      })
    })
  })

  output$data_summary_box <- renderUI({
    df <- rv$de
    if (is.null(df)) return(NULL)
    dfc <- classify_de(df, input$alpha, input$lfc_cutoff)
    contrast_line <- NULL
    # if (!is.null(input$treatment) && nzchar(input$treatment %||% "") &&
    #     !is.null(input$control) && nzchar(input$control %||% "")) {
    #   contrast_line <- paste0(input$treatment, " vs ", input$control)
    # }
    formula_line <- if (!is.null(rv$de_design_formula) && nzchar(rv$de_design_formula)) {
      paste0("\nDESeq2 formula: ", rv$de_design_formula)
    } else {
      NULL
    }
    if (!is.null(rv$de_contrast) && nzchar(rv$de_contrast)) {
      contrast_line <- paste0("Contrast: ", sub("^condition:\\s*", "", rv$de_contrast))
    }
    wellPanel(
      tags$strong("DE summary"),
      tags$pre(style = "margin: 8px 0 0 0; background: transparent; border: none; padding: 0;",
        paste(
          paste0("Genes: ", nrow(dfc)),
          paste0("Up: ", sum(dfc$DE_class == "up", na.rm = TRUE)),
          paste0("Down: ", sum(dfc$DE_class == "down", na.rm = TRUE)),
          formula_line,
          contrast_line,
          sep = "\n"
        )
      )
    )
  })

  output$all_comparison_venn_ui <- renderUI({
    if (is.null(rv$deseq_all_comparisons)) return(NULL)
    comparison_names <- names(rv$deseq_all_comparisons$tables %||% list())
    current_comparison <- paste0(input$treatment %||% "", "_vs_", input$control %||% "")
    default_comparisons <- if (current_comparison %in% comparison_names) {
      c(current_comparison, setdiff(comparison_names, current_comparison))
    } else {
      comparison_names
    }
    default_comparisons <- head(default_comparisons, 4)
    tagList(
      h4("Shared significant genes across all comparisons vs control"),
      div(class = "muted", "Uses the current padj and |log2FC| thresholds. The selected comparison below remains the treatment/control chosen in the left panel."),
      uiOutput("venn_selection_note"),
      fluidRow(
        column(8, plotOutput("all_comparison_venn_plot", width = "auto", height = "auto")),
        column(4,
          selectizeInput("venn_comparisons", "Comparisons to show",
            choices = stats::setNames(comparison_names, comparison_display_label(comparison_names)),
            selected = default_comparisons,
            multiple = TRUE,
            options = list(maxItems = 5, plugins = list("remove_button"))
          ),
          selectInput("venn_color_palette", "Venn color set",
            choices = c(
              "Default" = "default",
              "Okabe-Ito" = "okabe_ito",
              "Set 1" = "set1",
              "Set 2" = "set2",
              "Dark 2" = "dark2",
              "Paired" = "paired",
              "Tableau" = "tableau",
              "Viridis" = "viridis",
              "Plasma" = "plasma",
              "Pastel" = "pastel"
            ),
            selected = "default"
          ),
          tags$label("Circle line color", style = "font-weight: 600; font-size: 0.95em;"),
          tags$div(style = "display: flex; align-items: center; gap: 8px; margin-top: 4px;",
            tags$input(type = "color", id = "venn_line_color", value = "#404040",
              style = "width: 42px; height: 18px; border: none; cursor: pointer; padding: 0;",
              oninput = "Shiny.setInputValue('venn_line_color', this.value, {priority: 'event'})"),
            tags$span("Outline", style = "font-size: 0.9em;")
          ),
          actionButton("show_venn_intersection", "Show shared genes", class = "btn-info")
        )
      ),
      div(class = "download-row", downloadButton("download_all_deseq_comparisons", "Download all comparisons table")),
      tags$hr()
    )
  })

  output$venn_selection_note <- renderUI({
    req(rv$deseq_all_comparisons)
    comparison_names <- names(rv$deseq_all_comparisons$tables %||% list())
    selected <- input$venn_comparisons
    if (is.null(selected) || length(selected) == 0) {
      selected <- head(comparison_names, 4)
    }
    selected <- intersect(selected, comparison_names)
    div(class = "muted",
      paste0("Showing ", length(selected), " of ", length(comparison_names), " comparison",
        ifelse(length(comparison_names) == 1, "", "s"), " vs control. Default starts with 4; select up to 5 comparisons.")
    )
  })

  selected_venn_comparisons <- reactive({
    req(rv$deseq_all_comparisons)
    comparison_names <- names(rv$deseq_all_comparisons$tables %||% list())
    selected <- input$venn_comparisons
    if (is.null(selected) || length(selected) == 0) {
      selected <- head(comparison_names, 4)
    }
    head(intersect(selected, comparison_names), 5)
  })

  venn_intersection_table <- reactive({
    req(rv$deseq_all_comparisons)
    selected <- selected_venn_comparisons()
    if (length(selected) == 0) return(data.frame(gene_id = character(), stringsAsFactors = FALSE))

    sig_sets <- deseq_comparison_significant_sets(
      rv$deseq_all_comparisons,
      alpha = input$alpha,
      lfc_cutoff = input$lfc_cutoff
    )
    sig_sets <- sig_sets[selected]
    shared <- Reduce(intersect, sig_sets)
    shared <- sort(unique(shared))
    out <- data.frame(gene_id = shared, stringsAsFactors = FALSE)
    if (length(shared) == 0) return(out)

    for (comparison in selected) {
      tbl <- rv$deseq_all_comparisons$tables[[comparison]]
      tbl <- classify_de(tbl, alpha = input$alpha, lfc_cutoff = input$lfc_cutoff)
      keep_cols <- intersect(c("gene_id", "log2FoldChange", "padj", "DE_class"), names(tbl))
      tbl <- tbl[tbl$gene_id %in% shared, keep_cols, drop = FALSE]
      names(tbl)[names(tbl) != "gene_id"] <- paste0(comparison, "_", names(tbl)[names(tbl) != "gene_id"])
      out <- merge(out, tbl, by = "gene_id", all.x = TRUE, sort = FALSE)
    }
    out[match(shared, out$gene_id), , drop = FALSE]
  })

  observeEvent(input$show_venn_intersection, {
    req(rv$deseq_all_comparisons)
    selected <- selected_venn_comparisons()
    title <- if (length(selected) == 0) {
      "Shared genes"
    } else {
      paste("Shared genes:", paste(comparison_display_label(selected), collapse = " + "))
    }
    showModal(modalDialog(
      title = title,
      div(class = "muted",
        "Genes shown are significant in all selected Venn comparisons using the current padj and |log2FC| thresholds."
      ),
      br(),
      DTOutput("venn_intersection_table"),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  output$venn_intersection_table <- renderDT({
    d <- venn_intersection_table()
    datatable(d, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE))
  })

  output$all_comparison_venn_plot <- renderPlot({
    req(rv$deseq_all_comparisons)
    sig_sets <- deseq_comparison_significant_sets(
      rv$deseq_all_comparisons,
      alpha = input$alpha,
      lfc_cutoff = input$lfc_cutoff
    )
    selected_comparisons <- input$venn_comparisons
    if (!is.null(selected_comparisons) && length(selected_comparisons) > 0) {
      selected_comparisons <- head(intersect(selected_comparisons, names(sig_sets)), 5)
      if (length(selected_comparisons) > 0) {
        sig_sets <- sig_sets[selected_comparisons]
      }
    }
    make_comparison_venn_plot(
      sig_sets,
      max_sets = 5,
      title = paste0("Significant genes vs ", rv$deseq_all_comparisons$control),
      fill_palette = input$venn_color_palette %||% "default",
      line_color = input$venn_line_color %||% "#404040",
      selected_comparison = paste0(input$treatment %||% "", "_vs_", input$control %||% ""),
      font_family = input$plot_font_family %||% "serif"
    )
  }, width = function() input$venn_plot_width, height = function() input$venn_plot_height)

  output$de_preview <- renderDT({
    req(rv$de)
    keep_cols <- c("gene_id", "Symbol", "log2FoldChange", "padj")
    d <- rv$de[, intersect(keep_cols, names(rv$de)), drop = FALSE]
    # d <- head(d, 5000)
    button_data <- add_gene_count_buttons(d)
    datatable(button_data$data, rownames = FALSE, filter = "top", escape = !button_data$has_buttons,
              options = list(
                pageLength = 10,
                scrollX = TRUE,
                autoWidth = FALSE,
                columnDefs = if (button_data$has_buttons) gene_count_button_defs else list()
              ),
              callback = if (button_data$has_buttons) gene_count_button_callback else JS(""))
  })
  
  output$search_annotations_table <- renderDT({
    req(rv$de)
    cols <- names(rv$de)
    first_cols <- intersect(c("gene_id", "Symbol", "log2FoldChange", "padj"), cols)
    other_cols <- setdiff(cols, first_cols)
    d <- rv$de[, c(first_cols, other_cols), drop = FALSE]
    rv$search_annotations_modal_data <- d
    button_data <- add_gene_count_buttons(d)
    d <- add_row_detail_buttons(button_data$data, "search_annotations_table")
    has_count_buttons <- button_data$has_buttons
    search_column_defs <- list(
      list(
        targets = "_all",
        render = JS(
          "function(data, type, row, meta) {",
          "  if (type === 'display' && data != null && typeof data === 'string' && (data.indexOf('gene-counts-btn') !== -1 || data.indexOf('row-detail-btn') !== -1)) {",
          "    return data;",
          "  }",
          "  if (type === 'display' && data != null && typeof data === 'string' && data.length > 50) {",
          "    return '<span title=\"' + data.replace(/\"/g, '&quot;') + '\">' + data.substr(0, 50) + '...</span>';",
          "  } else {",
          "    return data;",
          "  }",
          "}"
        )
      )
    )
    search_column_defs <- c(row_detail_button_defs(0), search_column_defs)
    if (isTRUE(has_count_buttons)) {
      search_column_defs <- c(search_column_defs, gene_count_button_defs)
    }
    
    # Use JavaScript to truncate long text cells and show full text on hover
    datatable(d, rownames = FALSE, filter = "top", selection = "none", escape = FALSE,
      options = list(
        search = list(regex = TRUE, smart = TRUE),
        pageLength = 15, 
        scrollX = TRUE,
        columnDefs = search_column_defs
      ),
      callback = JS(paste(c(
        as.character(row_detail_button_callback),
        if (isTRUE(has_count_buttons)) as.character(gene_count_button_callback) else ""
      ), collapse = "\n"))
    )
  })

  observeEvent(input$row_detail_clicked, {
    click <- input$row_detail_clicked
    table_id <- click$table %||% ""
    row_i <- suppressWarnings(as.integer(click$row %||% NA_integer_))
    if (is.na(row_i) || row_i < 1) return()

    if (identical(table_id, "search_annotations_table")) {
      d <- rv$search_annotations_modal_data
      title <- "Annotation search row details"
    } else if (identical(table_id, "annotation_preview")) {
      d <- rv$annotation_preview_modal_data
      title <- "Annotation preview row details"
    } else {
      return()
    }

    if (is.null(d) || row_i > nrow(d)) return()
    show_row_details_modal(d[row_i, , drop = FALSE], title)
  }, ignoreInit = TRUE)

  output$gene_counts_click_hint <- renderUI({
    if (!identical(input$data_mode, "rsem") || is.null(rv$norm_counts)) return(NULL)
    div(class = "muted", "Use the Plot button beside a gene to view normalized counts for the selected treatment/control comparison.")
  })

  add_gene_count_buttons <- function(d) {
    d <- as.data.frame(d, check.names = FALSE)
    has_count_buttons <- identical(input$data_mode, "rsem") && !is.null(rv$norm_counts) && "gene_id" %in% names(d)
    if (!isTRUE(has_count_buttons)) {
      return(list(data = d, has_buttons = FALSE))
    }
    count_buttons <- vapply(as.character(d$gene_id), function(gene_id) {
      gene_id_escaped <- htmltools::htmlEscape(gene_id, attribute = TRUE)
      paste0(
        "<button type=\"button\" class=\"btn btn-link btn-xs gene-counts-btn\" style=\"font-size: 18px; padding: 0; line-height: 1;\" ",
        "data-gene=\"", gene_id_escaped, "\">▣</button>"
      )
    }, character(1))
    list(data = data.frame(" " = count_buttons, d, check.names = FALSE), has_buttons = TRUE)
  }

  gene_count_button_callback <- JS(
    "table.on('click', 'button.gene-counts-btn', function(e) {",
    "  e.stopPropagation();",
    "  Shiny.setInputValue('gene_counts_gene_clicked', this.getAttribute('data-gene'), {priority: 'event'});",
    "});"
  )

  gene_count_button_defs <- list(
    list(targets = 0, orderable = FALSE, searchable = FALSE, width = "10px", className = "dt-left")
  )

  selected_gene_counts_plot_reactive <- reactive({
    req(rv$selected_gene_counts_gene, rv$norm_counts, rv$coldata)
    make_gene_norm_counts_boxplot(
      rv$norm_counts,
      rv$coldata,
      rv$selected_gene_counts_gene,
      treatment = input$treatment %||% NULL,
      control = input$control %||% NULL,
      plot_theme = input$plot_theme %||% "classic",
      font_family = input$plot_font_family %||% "serif",
      color_trnt = input$gene_counts_color_trnt %||% "#ac783e",
      color_ctrl = input$gene_counts_color_ctrl %||% "#505050",
      point_size = input$gene_counts_point_size %||% 2.4,
      point_alpha = input$gene_counts_point_alpha %||% 0.88,
      jitter_width = input$gene_counts_jitter_width %||% 0.12,
      box_width = input$gene_counts_box_width %||% 0.55
    )
  })

  output$selected_gene_counts_plot_ui <- renderUI({
    height <- input$gene_counts_plot_height %||% 420
    plotOutput("selected_gene_counts_plot", width = "100%", height = paste0(height, "px"))
  })

  output$selected_gene_counts_plot <- renderPlot({
    selected_gene_counts_plot_reactive()
  }, width = function() input$gene_counts_plot_width %||% 650,
     height = function() input$gene_counts_plot_height %||% 420)

  output$download_gene_counts <- download_plot_server(
    selected_gene_counts_plot_reactive,
    reactive(input$format_gene_counts %||% "png"),
    reactive(paste0("normalized_counts_", safe_filename_part(rv$selected_gene_counts_gene %||% "gene", "gene"))),
    reactive(input$gene_counts_plot_width %||% 650),
    reactive(input$gene_counts_plot_height %||% 420)
  )

  observeEvent(input$gene_counts_gene_clicked, {
    if (!identical(input$data_mode, "rsem") || is.null(rv$norm_counts) || is.null(rv$coldata)) return()
    gene_id <- input$gene_counts_gene_clicked
    if (is.na(gene_id) || !nzchar(gene_id)) return()
    rv$selected_gene_counts_gene <- gene_id
    showModal(modalDialog(
      title = paste("Normalized counts:", gene_id),
      fluidRow(
        column(3,
          numericInput("gene_counts_plot_width", "Width (px)", value = 250, min = 100, max = 1600, step = 25)
        ),
        column(3,
          numericInput("gene_counts_plot_height", "Height (px)", value = 200, min = 100, max = 1200, step = 25)
        ),
        column(3,
          sliderInput("gene_counts_point_size", "Point size", min = 0.1, max = 8, value = 2.4, step = 0.1)
        ),
        column(3,
          sliderInput("gene_counts_point_alpha", "Point transparency", min = 0, max = 1, value = 0.88, step = 0.05)
        )
      ),
      fluidRow(
        column(3,
          tags$label("Treatment color"),
          tags$div(style = "display: flex; align-items: center; gap: 8px;",
            tags$input(type = "color", id = "gene_counts_color_trnt", value = "#ac783e",
              style = "width: 42px; height: 26px; border: none; cursor: pointer; padding: 0;",
              oninput = "Shiny.setInputValue('gene_counts_color_trnt', this.value, {priority: 'event'})")
          )
        ),
        column(3,
          tags$label("Control color"),
          tags$div(style = "display: flex; align-items: center; gap: 8px;",
            tags$input(type = "color", id = "gene_counts_color_ctrl", value = "#505050",
              style = "width: 42px; height: 26px; border: none; cursor: pointer; padding: 0;",
              oninput = "Shiny.setInputValue('gene_counts_color_ctrl', this.value, {priority: 'event'})")
          )
        ),
        column(3,
          sliderInput("gene_counts_jitter_width", "Jitter width", min = 0, max = 0.5, value = 0.12, step = 0.01)
        ),
        column(3,
          sliderInput("gene_counts_box_width", "Box width", min = 0.2, max = 0.95, value = 0.55, step = 0.05)
        )
      ),
      uiOutput("selected_gene_counts_plot_ui"),
      download_plot_ui("gene_counts", "Download counts plot"),
      easyClose = TRUE,
      size = "l",
      footer = modalButton("Close")
    ))
  }, ignoreInit = TRUE)

  output$annotation_summary <- renderText({
    ann <- rv$annotation_df
    if (is.null(ann)) {
      return("No annotation table loaded. Search an organism and build from UniProt, or load a manual description file.")
    }
    match_count <- annotation_match_count()
    paste0(
      "Source: ", rv$annotation_label %||% "annotation table", "\n",
      "Annotation rows: ", nrow(ann), "\n",
      "Annotation columns: ", max(ncol(ann) - 2, 0), "\n",
      "Replace result gene_id with annotation gene_id: ", isTRUE(rv$annotation_replace_gene_id), "\n",
      "Loaded genes matched: ", ifelse(is.na(match_count), "load a DE table first", match_count), "\n",
      "Applied to DE table: ", !is.null(rv$de)
    )
  })

  output$annotation_preview <- renderDT({
    req(rv$annotation_preview_ready, rv$annotation_df)
    d <- rv$annotation_df
    d$annotation_key <- NULL
    d <- head(d, 5000)
    rv$annotation_preview_modal_data <- d
    display_d <- add_row_detail_buttons(d, "annotation_preview")
    datatable(display_d, rownames = FALSE, filter = "top", selection = "none", escape = FALSE,
              options = list(
                pageLength = 12,
                scrollX = TRUE,
                columnDefs = row_detail_button_defs(0)
              ),
              callback = row_detail_button_callback)
  })

  output$annotation_preview_ui <- renderUI({
    if (!isTRUE(rv$annotation_preview_ready) || is.null(rv$annotation_df)) return(NULL)
    tagList(
      h4("Annotation preview"),
      DTOutput("annotation_preview")
    )
  })

  output$uniprot_id_source_ui <- renderUI({
    choices <- rv$uniprot_id_choices
    if (is.null(choices) || length(choices) == 0) {
      return(selectInput("uniprot_id_source", "Use as gene_id", choices = c("Scan ID sources first" = ""), selected = ""))
    }

    selected <- if ("Gene.Names..ordered.locus." %in% unname(choices)) {
      "Gene.Names..ordered.locus."
    } else if ("GeneID" %in% unname(choices)) {
      "GeneID"
    } else if ("RefSeq" %in% unname(choices)) {
      "RefSeq"
    } else if ("Entry" %in% unname(choices)) {
      "Entry"
    } else {
      unname(choices)[1]
    }
    selectInput("uniprot_id_source", "Use as gene_id", choices = choices, selected = selected)
  })

  output$uniprot_status <- renderText({
    rv$uniprot_status %||% ""
  })

  observeEvent(input$annotation_reviewed_only, {
    rv$uniprot_id_choices <- NULL
    rv$uniprot_status <- "Reviewed-only setting changed. Scan UniProt ID sources again before building."
  }, ignoreInit = TRUE)

  observeEvent(input$scan_uniprot_id_sources, {
    tax_id <- suppressWarnings(as.integer(rv$selected_tax_id))
    validate(need(!is.na(tax_id) && tax_id > 0, "Enter a valid NCBI taxonomy ID."))
    append_log("Scanning UniProt ID sources for tax_id", tax_id, level = "STEP")
    rv$uniprot_status <- paste("Scanning UniProt ID sources for tax_id", tax_id, "...")
    withProgress(message = "Scanning UniProt ID sources", value = 0.1, {
      tryCatch({
        incProgress(0.4, detail = paste("Downloading UniProt IDs for tax_id", tax_id))
        choices <- uniprot_id_source_choices(
          tax_id = tax_id,
          reviewed_only = isTRUE(input$annotation_reviewed_only),
          max_records = 500
        )
        if (length(choices) == 0) stop("No usable UniProt ID sources were found for tax_id ", tax_id, ".")
        rv$uniprot_id_choices <- choices
        rv$uniprot_status <- paste("Found", length(choices), "UniProt ID sources. Choose one, then build the description file.")
        showNotification(paste("Found", length(choices), "UniProt ID sources."), type = "message", duration = 5)
        append_log("Found", length(choices), "UniProt ID sources for tax_id", tax_id)
      }, error = function(e) {
        rv$uniprot_id_choices <- NULL
        rv$uniprot_status <- paste("UniProt ID-source scan error:", e$message)
        showNotification(e$message, type = "error", duration = 15)
        append_log("UniProt ID-source scan error:", e$message)
      })
    })
  })

  output$refseq_gtf_id_source_ui <- renderUI({
    choices <- rv$refseq_gtf_id_choices
    if (is.null(choices) || length(choices) == 0) {
      return(selectInput("refseq_gtf_id_source", "Use as gene_id", choices = c("Upload a GTF file first" = ""), selected = ""))
    }

    selected <- if ("db_xref_GenBank" %in% unname(choices)) {
      "db_xref_GenBank"
    } else if ("protein_id" %in% unname(choices)) {
      "protein_id"
    } else {
      unname(choices)[1]
    }
    selectInput("refseq_gtf_id_source", "Use as gene_id", choices = choices, selected = selected)
  })

  output$refseq_gtf_status <- renderText({
    rv$refseq_gtf_status %||% ""
  })

  observeEvent(input$refseq_gtf_file, {
    req(input$refseq_gtf_file$datapath)
    append_log("Scanning RefSeq GTF ID sources:", input$refseq_gtf_file$name, level = "STEP")
    rv$refseq_gtf_status <- paste("Scanning ID sources in", input$refseq_gtf_file$name, "...")
    withProgress(message = "Scanning RefSeq GTF", value = 0.2, {
      tryCatch({
        incProgress(0.3, detail = "Reading GTF attributes")
        choices <- refseq_gtf_id_source_choices(input$refseq_gtf_file$datapath, max_rows = 50000)
        if (length(choices) == 0) stop("No usable GTF attributes or db_xref IDs were found.")
        rv$refseq_gtf_id_choices <- choices
        rv$refseq_gtf_status <- paste("Found", length(choices), "ID sources. Choose one, then build the description file.")
        showNotification(paste("Found", length(choices), "RefSeq GTF ID sources."), type = "message", duration = 5)
        append_log("Found", length(choices), "RefSeq GTF ID sources in", input$refseq_gtf_file$name)
      }, error = function(e) {
        rv$refseq_gtf_id_choices <- NULL
        rv$refseq_gtf_status <- paste("RefSeq GTF scan error:", e$message)
        showNotification(e$message, type = "error", duration = 12)
        append_log("RefSeq GTF scan error:", e$message)
      })
    })
  }, ignoreInit = TRUE)

  observeEvent(input$load_annotation_file, {
    req(input$annotation_file$datapath)
    append_log("Loading annotation table:", input$annotation_file$name, level = "STEP")
    tryCatch({
      ann <- normalize_annotation_table(read_any_table(input$annotation_file$datapath))
      rv$annotation_df <- ann
      rv$annotation_label <- paste("Uploaded file:", input$annotation_file$name)
      rv$annotation_replace_gene_id <- FALSE
      rv$annotation_preview_ready <- TRUE
      if (!is.null(rv$de_base)) apply_current_annotation(reset_results = TRUE)
      append_log("Loaded annotation table:", nrow(ann), "rows from", input$annotation_file$name)
    }, error = function(e) {
      showNotification(e$message, type = "error", duration = 12)
      append_log("Annotation load error:", e$message)
    })
  })

  observeEvent(input$build_refseq_gtf_annotations, {
    req(input$refseq_gtf_file$datapath)
    id_source <- input$refseq_gtf_id_source %||% ""
    validate(need(nzchar(id_source), "Upload a GTF file and select an ID source."))
    append_log("Building RefSeq GTF annotation table from", input$refseq_gtf_file$name, "using", id_source, level = "STEP")
    rv$refseq_gtf_status <- paste("Building description file using", id_source, "...")
    withProgress(message = "Building RefSeq GTF description file", value = 0.1, {
      tryCatch({
        incProgress(0.2, detail = "Parsing GTF attributes and db_xref IDs")
        ann <- build_refseq_gtf_description_file(
          input$refseq_gtf_file$datapath,
          gene_id_source = id_source,
          one_row_per_id = isTRUE(input$refseq_gtf_one_row_per_id)
        )
        incProgress(0.4, detail = "Normalizing annotation table")
        ann <- normalize_annotation_table(ann)
        rv$annotation_df <- ann
        rv$annotation_label <- paste0("RefSeq GTF ", input$refseq_gtf_file$name, " using ", id_source)
        rv$annotation_replace_gene_id <- TRUE
        rv$annotation_preview_ready <- TRUE

        detected <- detect_gene_id_type_from_values(ann$gene_id)
        if (!is.null(detected) && !is.null(detected$keytype) && detected$keytype %in% gene_id_type_choices_for_orgdb(input$go_orgdb %||% NULL)) {
          rv$detected_gene_id_type <- detected$keytype
          rv$detected_gene_id_type_confidence <- detected$confidence
          updateSelectInput(session, "go_keytype", selected = detected$keytype)
        }

        incProgress(0.3, detail = "Applying annotations to DE table")
        if (!is.null(rv$de_base)) apply_current_annotation(reset_results = TRUE)
        rv$refseq_gtf_status <- paste("Built", nrow(ann), "annotation rows from", input$refseq_gtf_file$name, "using", id_source)
        showNotification(paste("Built RefSeq GTF annotation table:", nrow(ann), "rows."), type = "message", duration = 6)
        append_log("Built RefSeq GTF annotation table:", nrow(ann), "rows from", input$refseq_gtf_file$name)
      }, error = function(e) {
        rv$refseq_gtf_status <- paste("RefSeq GTF annotation error:", e$message)
        showNotification(e$message, type = "error", duration = 15)
        append_log("RefSeq GTF annotation error:", e$message)
      })
    })
  })

  output$taxonomy_search_status <- renderText({
    rv$taxonomy_status
  })

  output$gene_id_type_detection <- renderText({
    if (is.null(rv$detected_gene_id_type)) {
      return("Auto-detects after loading a DE table; choose manually if needed.")
    }
    paste0(
      "Auto-detected from Data input: ",
      rv$detected_gene_id_type,
      if (!is.null(rv$detected_gene_id_type_confidence)) {
        paste0(" (", round(100 * rv$detected_gene_id_type_confidence), "% match)")
      } else {
        ""
      }
    )
  })

  output$organism_availability_ui <- renderUI({
    status_item <- function(ok, label, extra = NULL) {
      div(style = "display:inline-block; margin-right:18px;",
          strong(if (isTRUE(ok)) "✔️" else "×"),
          span(paste0(" ", label, if (isTRUE(ok) && !is.null(extra) && nzchar(extra)) paste0(": ", extra) else "")))
    }
    tagList(
      div(class = "muted", paste0("Selected: ", rv$selected_organism_label, " (", rv$selected_tax_id, ")")),
      div(style = "margin-top:6px;",
        status_item(rv$selected_uniprot_available, "UniProt"),
        status_item(rv$selected_go_available, "GO", input$go_orgdb %||% ""),
        status_item(rv$selected_kegg_available, "KEGG", input$kegg_species %||% ""),
        status_item(rv$selected_pmn_available, "PMN", input$pmn_cyc_db %||% "")
      )
    )
  })

  output$go_selected_organism <- renderText({
    paste0("Organism: ", rv$selected_organism_label, "\nTax ID: ", rv$selected_tax_id)
  })

  output$kegg_selected_organism <- renderText({
    paste0("Organism: ", rv$selected_organism_label, "\nTax ID: ", rv$selected_tax_id)
  })

  output$arabidopsis_default_ui <- renderUI({
    tax_id <- suppressWarnings(as.integer(rv$selected_tax_id))
    url <- default_description_file_url(tax_id)
    if (is.na(tax_id) || is.null(url)) return(NULL)
    label <- if (identical(tax_id, 3702L)) {
      "Use default Arabidopsis description file"
    } else if (identical(tax_id, 9606L)) {
      "Use default human description file"
    } else if (identical(tax_id, 511145L)) {
      "Use default E. coli K-12 MG1655 description file"
    } else {
      "Use default organism description file"
    }
    actionButton("use_default_organism_desc", label,
                 class = "btn-info", style = "width:100%; margin-bottom: 10px;")
  })

  search_and_fill_taxonomy <- function(query) {
    append_log("Searching UniProt taxonomy:", query, level = "STEP")
    tryCatch({
      hits <- search_uniprot_taxonomy(query, size = 50)
      rv$taxonomy_results <- hits
      choices <- stats::setNames(as.character(hits$tax_id), hits$label)
      updateSelectizeInput(session, "annotation_taxon_choice", choices = choices, selected = as.character(hits$tax_id[1]))
      apply_selected_taxon_to_analysis(hits$tax_id[1], hits$scientific_name[1])
      rv$taxonomy_status <- paste("Found", nrow(hits), "UniProt taxonomy matches. Selected:", hits$label[1])
      append_log("UniProt taxonomy search:", query, "-", nrow(hits), "matches.")
    }, error = function(e) {
      rv$taxonomy_results <- NULL
      rv$taxonomy_status <- e$message
      updateSelectizeInput(session, "annotation_taxon_choice", choices = default_taxon_choices, selected = "3702")
      apply_selected_taxon_to_analysis(3702, "Arabidopsis thaliana")
      showNotification(e$message, type = "error", duration = 12)
      append_log("UniProt taxonomy search error:", e$message)
    })
  }

  observeEvent(input$annotation_taxon_choice, {
    choice <- input$annotation_taxon_choice
    if (!is.null(choice) && nzchar(choice)) {
      tax_id <- suppressWarnings(as.integer(choice))
      if (is.na(tax_id)) {
        search_and_fill_taxonomy(choice)
      } else if (is.null(organism_config_for_tax(tax_id)) && is.null(taxonomy_label_for_tax(tax_id))) {
        search_and_fill_taxonomy(choice)
      } else {
        apply_selected_taxon_to_analysis(tax_id)
      }
    }
  }, ignoreInit = TRUE)

  observeEvent(input$go_orgdb, {
    cfg <- organism_analysis_config[organism_analysis_config$orgdb == input$go_orgdb, , drop = FALSE]
    go_available <- ensure_orgdb_installed(input$go_orgdb)
    keytype_choices <- gene_id_type_choices_for_orgdb(input$go_orgdb)
    selected_keytype <- input$go_keytype %||% rv$detected_gene_id_type %||% "ENTREZID"
    if (nrow(cfg) > 0) {
      selected_keytype <- rv$detected_gene_id_type %||% selected_keytype %||% cfg$go_keytype[1]
      if (!selected_keytype %in% unname(keytype_choices)) selected_keytype <- cfg$go_keytype[1]
    }
    if (!selected_keytype %in% unname(keytype_choices)) selected_keytype <- unname(keytype_choices)[1]
    updateSelectInput(session, "go_keytype", choices = keytype_choices, selected = selected_keytype)
    rv$selected_go_available <- go_available
    rv$go_cache <- list()
  }, ignoreInit = TRUE)

  observeEvent(input$kegg_species, {
    rv$selected_kegg_available <- !is.null(input$kegg_species) && nzchar(input$kegg_species)
  }, ignoreInit = TRUE)

  observeEvent(input$use_default_organism_desc, {
    tax_id <- suppressWarnings(as.integer(rv$selected_tax_id))
    url <- default_description_file_url(tax_id)
    validate(need(!is.null(url), "No bundled default description file is configured for this organism."))
    label <- default_annotation_label_for_tax(tax_id)
    append_log("Loading", label, level = "STEP")
    tryCatch({
      ann <- load_description_file(url)
      if (is.null(ann) || nrow(ann) == 0) stop("Could not load annotation table from: ", url)
      rv$annotation_df <- ann
      rv$annotation_label <- label
      rv$annotation_replace_gene_id <- FALSE
      rv$annotation_preview_ready <- TRUE
      if (!is.null(rv$de_base)) apply_current_annotation(reset_results = TRUE)
      append_log("Loaded", label, ":", nrow(ann), "rows.")
    }, error = function(e) {
      showNotification(e$message, type = "error", duration = 12)
      append_log("Default annotation error:", e$message)
    })
  })

  observeEvent(input$build_uniprot_annotations, {
    tax_id <- suppressWarnings(as.integer(rv$selected_tax_id))
    validate(need(!is.na(tax_id) && tax_id > 0, "Enter a valid NCBI taxonomy ID."))
    id_source <- input$uniprot_id_source %||% ""
    validate(need(nzchar(id_source), "Scan UniProt ID sources and select one before building."))
    append_log("Building UniProt annotation table for tax_id", tax_id, "using", id_source, level = "STEP")
    rv$uniprot_status <- paste("Building UniProt description file using", id_source, "...")
    withProgress(message = "Building UniProt description file", value = 0.1, {
      tryCatch({
        incProgress(0.2, detail = paste("Downloading UniProt annotations for tax_id", tax_id))
        ann <- build_uniprot_description_file(
          input_data = NULL,
          tax_id = tax_id,
          reviewed_only = isTRUE(input$annotation_reviewed_only),
          gene_id_source = id_source,
          gene_id_type = input$go_keytype %||% NULL,
          orgdb = input$go_orgdb %||% NULL
        )
        rv$annotation_label <- paste0("UniProt tax_id ", tax_id, " using ", id_source, if (isTRUE(input$annotation_reviewed_only)) " reviewed only" else "")
        ann <- normalize_annotation_table(ann)
        incProgress(0.6, detail = "Applying annotations to DE table")
        rv$annotation_df <- ann
        rv$annotation_replace_gene_id <- TRUE
        rv$annotation_preview_ready <- TRUE

        detected <- detect_gene_id_type_from_values(ann$gene_id)
        if (!is.null(detected) && !is.null(detected$keytype) && detected$keytype %in% gene_id_type_choices_for_orgdb(input$go_orgdb %||% NULL)) {
          rv$detected_gene_id_type <- detected$keytype
          rv$detected_gene_id_type_confidence <- detected$confidence
          updateSelectInput(session, "go_keytype", selected = detected$keytype)
        }

        if (!is.null(rv$de_base)) apply_current_annotation(reset_results = TRUE)
        match_count <- annotation_match_count()
        rv$uniprot_status <- paste("Built", nrow(ann), "UniProt annotation rows using", id_source)
        showNotification(paste("Built UniProt annotation table:", nrow(ann), "rows."), type = "message", duration = 6)
        append_log("Built UniProt annotation table for tax_id", tax_id, "using", id_source, ":", ifelse(is.na(match_count), 0, match_count), "matched genes.")
      }, error = function(e) {
        rv$uniprot_status <- paste("UniProt annotation error:", e$message)
        showNotification(e$message, type = "error", duration = 15)
        append_log("UniProt annotation error:", e$message)
      })
    })
  })

  volcano_reactive <- reactive({
    req(rv$de)
    make_volcano_plot(rv$de, input$alpha, input$lfc_cutoff, "Volcano plot",
                      point_size = input$plot_point_size, point_alpha = input$plot_alpha,
                      color_up = input$color_up %||% "#B2182B", color_down = input$color_down %||% "#2166AC", color_ns = input$color_ns %||% "#B3B3B3",
                      plot_theme = input$plot_theme %||% "classic", font_family = input$plot_font_family %||% "serif")
  })
  output$volcano_plot <- renderPlot({ volcano_reactive() }, width = function() input$de_plot_width, height = function() input$de_plot_height)

  ma_reactive <- reactive({
    req(rv$de)
    make_ma_plot(rv$de, input$alpha, input$lfc_cutoff, "MA plot",
                 point_size = input$plot_point_size, point_alpha = input$plot_alpha,
                 color_up = input$color_up %||% "#B2182B", color_down = input$color_down %||% "#2166AC", color_ns = input$color_ns %||% "#B3B3B3",
                 plot_theme = input$plot_theme %||% "classic", font_family = input$plot_font_family %||% "serif")
  })
  output$ma_plot <- renderPlot({
    tryCatch(ma_reactive(), error = function(e) {
      plot.new(); text(0.5, 0.5, e$message, cex = 1.1)
    })
  }, width = function() input$de_plot_width, height = function() input$de_plot_height)
  output$ma_message <- renderText({
    if (is.null(rv$de)) return("")
    if (!"baseMean" %in% names(rv$de)) "MA plot requires baseMean. Run DESeq2 from RSEM files or upload a DE table containing baseMean." else ""
  })

  # PCA conditions selector (default = treatment + control from DESeq2 run)
  output$pca_conditions_ui <- renderUI({
    req(rv$pca)
    all_conditions <- unique(rv$pca$condition)
    # Default to the treatment and control used for the comparison
    default <- all_conditions
    if (!is.null(input$treatment) && !is.null(input$control)) {
      used <- intersect(c(input$treatment, input$control), all_conditions)
      if (length(used) > 0) default <- used
    }
    selectizeInput("pca_conditions", "PCA conditions", choices = all_conditions,
                   selected = default, multiple = TRUE,
                   options = list(placeholder = "Select conditions to display"))
  })

  pca_reactive <- reactive({
    req(rv$pca)
    make_pca_plot(rv$pca, "PCA", point_size = input$plot_point_size * 3.2, point_alpha = input$plot_alpha,
                  show_labels = isTRUE(input$pca_show_labels),
                  conditions = input$pca_conditions,
                  plot_theme = input$plot_theme %||% "classic", font_family = input$plot_font_family %||% "serif",
                  color_palette = input$pca_color_palette %||% "default")
  })
  output$pca_plot <- renderPlot({
    if (is.null(rv$pca)) {
      plot.new(); text(0.5, 0.5, "PCA requires expression/count data. Run DESeq2 from RSEM files.", cex = 1.1)
    } else pca_reactive()
  }, width = function() input$de_plot_width, height = function() input$de_plot_height)
  output$pca_message <- renderText({
    if (is.null(rv$pca)) "PCA is not available from a DE-only CSV." else ""
  })

  observeEvent(input$run_go, {
    req(rv$de)
    append_log("Running GO enrichment for", input$go_direction, "genes", level = "STEP")
    withProgress(message = "Running GO enrichment", value = 0.2, {
      tryCatch({
        incProgress(0.2, detail = "GO enrichment")
        g <- compute_go_cached(input$go_direction)  # stores all terms in cache
        rv$go_trigger <- rv$go_trigger + 1
        append_log("GO enrichment done for", input$go_direction, "genes:", nrow(g), "total terms (filtering by p-value on display).")
      }, error = function(e) {
        showNotification(e$message, type = "error", duration = 12)
        append_log("GO error:", e$message)
      })
    })
  })

  output$go_bubble <- renderPlot({
    d <- tryCatch(go_display_data(), error = function(e) NULL)
    if (is.null(d) || nrow(d) == 0) {
      plot.new(); text(0.5, 0.5, if (rv$go_trigger == 0) "Run GO enrichment first." else "No terms pass the current p-value cutoff.", cex = 1.1)
    } else tryCatch(go_display_plot(), error = function(e) { plot.new(); text(0.5, 0.5, e$message, cex = 1.1) })
  }, width = function() input$go_plot_width, height = function() input$go_plot_height)
  output$go_table <- renderDT({
    d <- tryCatch(go_display_data(), error = function(e) NULL)
    req(!is.null(d) && nrow(d) > 0)
    datatable(d, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE))
  })

  output$go_gene_codes_ui <- renderUI({
    d <- tryCatch(go_display_data(), error = function(e) NULL)
    if ((is.null(d) || nrow(d) == 0) && rv$go_trigger > 0) {
      d <- tryCatch(get_cached_go(input$go_direction), error = function(e) NULL)
    }
    if (!is.null(d) && nrow(d) > 0 && "GO.ID" %in% names(d)) {
      ids <- stats::na.omit(as.character(d$GO.ID))
      ids <- unique(ids[nzchar(ids)])
      labels <- ids
      if ("Term" %in% names(d)) {
        terms_by_id <- stats::setNames(as.character(d$Term), as.character(d$GO.ID))
        labels <- paste0(ids, " - ", terms_by_id[ids])
        labels[is.na(labels) | !nzchar(labels)] <- ids[is.na(labels) | !nzchar(labels)]
      }
      return(selectizeInput(
        "go_gene_codes",
        "GO ID(s)",
        choices = stats::setNames(ids, labels),
        selected = head(ids, 3),
        multiple = TRUE,
        options = list(placeholder = "Type to search GO IDs", plugins = list("remove_button"))
      ))
    }
    selectizeInput(
      "go_gene_codes",
      "GO ID(s)",
      choices = c("GO:0008150"),
      selected = "GO:0008150",
      multiple = TRUE,
      options = list(create = TRUE, placeholder = "Example: GO:0008150", plugins = list("remove_button"))
    )
  })

  observeEvent(input$run_go_gene_lookup, {
    req(rv$de)
    append_log("Looking up GO genes for", input$go_gene_codes, level = "STEP")
    withProgress(message = "Looking up GO genes", value = 0.2, {
      tryCatch({
        genes <- make_go_gene_table(
          rv$de,
          go_ids = input$go_gene_codes,
          ontology = input$ontology %||% "BP",
          orgdb = input$go_orgdb %||% "org.At.tair.db",
          keytype = input$go_keytype %||% "TAIR",
          alpha = input$alpha,
          lfc_cutoff = input$lfc_cutoff
        )
        rv$go_genes <- genes
        if (nrow(genes) > 0) {
          rv$go_gene_lookup <- stats::aggregate(
            gene_id ~ GO_ID + GO_term,
            data = genes,
            FUN = function(x) length(unique(x))
          )
          names(rv$go_gene_lookup)[names(rv$go_gene_lookup) == "gene_id"] <- "Matched_loaded_genes"
        } else {
          rv$go_gene_lookup <- data.frame()
        }
        append_log("GO gene lookup completed:", nrow(rv$go_genes), "matched gene rows.")
      }, error = function(e) {
        rv$go_gene_lookup <- NULL
        rv$go_genes <- NULL
        showNotification(paste("GO gene lookup error:", e$message), type = "error", duration = 15)
        append_log("GO gene lookup error:", e$message)
      })
    })
  })

  output$go_genes_table <- renderDT({
    req(rv$go_genes)
    d <- rv$go_genes
    req(nrow(d) > 0)
    keep_cols <- c("gene_id", "original_gene_id", "GO_ID", "GO_term", "Symbol", "log2FoldChange", "pValue", "padj", "DE_class", "short_description", "Protein_name")
    d <- d[, intersect(keep_cols, names(d)), drop = FALSE]
    d <- d[order(d$GO_ID, d$padj, d$pValue, na.last = TRUE), , drop = FALSE]
    d <- d[!duplicated(d[, intersect(c("gene_id", "GO_ID"), names(d)), drop = FALSE]), , drop = FALSE]
    button_data <- add_gene_count_buttons(d)
    datatable(button_data$data, rownames = FALSE, filter = "top", escape = !button_data$has_buttons,
              options = list(
                pageLength = 15,
                scrollX = TRUE,
                autoWidth = FALSE,
                columnDefs = if (button_data$has_buttons) gene_count_button_defs else list()
              ),
              callback = if (button_data$has_buttons) gene_count_button_callback else JS(""))
  })

  go_genes_volcano_reactive <- reactive({
    req(rv$go_genes)
    d <- rv$go_genes
    req(nrow(d) > 0)
    plot_df <- d[order(d$padj, d$pValue, na.last = TRUE), , drop = FALSE]
    plot_df <- plot_df[!duplicated(plot_df$gene_id), , drop = FALSE]
    go_label <- paste(unique(d$GO_ID), collapse = ", ")
    if (nchar(go_label) > 80) go_label <- paste0(substr(go_label, 1, 77), "...")
    make_gene_group_volcano_plot(
      plot_df,
      paste0("GO genes: ", go_label),
      alpha = input$alpha,
      lfc_cutoff = input$lfc_cutoff,
      color_up = input$color_up %||% "#B2182B",
      color_down = input$color_down %||% "#2166AC",
      color_ns = input$color_ns %||% "#B3B3B3",
      plot_theme = input$plot_theme %||% "classic",
      font_family = input$plot_font_family %||% "serif"
    )
  })

  output$go_genes_volcano <- renderPlot({
    if (is.null(rv$go_gene_lookup)) {
      plot.new(); text(0.5, 0.5, "Enter GO ID(s) and click 'Find genes'.", cex = 1.1)
    } else if (is.null(rv$go_genes) || nrow(rv$go_genes) == 0) {
      plot.new(); text(0.5, 0.5, "No loaded genes matched the selected GO term(s).", cex = 1.1)
    } else {
      tryCatch(go_genes_volcano_reactive(), error = function(e) {
        plot.new(); text(0.5, 0.5, e$message, cex = 1.1)
      })
    }
  }, width = function() input$go_plot_width, height = function() input$go_plot_height)

  output$msigdb_run_ui <- renderUI({
    if (length(msigdb_species_choices) == 0) {
      return(div(class = "muted", "Install msigdbr to enable Hallmark analysis."))
    }
    if (!msigdb_species_available(input$msigdb_species)) {
      return(div(class = "muted", paste0(rv$selected_organism_label, " is not available in MSigDB. Select another available MSigDB species to run.")))
    }
    actionButton("run_msigdb", "Run Hallmark analysis", class = "btn-primary", style = "width:100%;")
  })

  observeEvent(input$run_msigdb, {
    req(rv$de)
    validate(need(msigdb_species_available(input$msigdb_species), "Select an available MSigDB species."))
    append_log("Running MSigDB Hallmark analysis for", input$msigdb_direction, "genes", level = "STEP")
    withProgress(message = "Running MSigDB Hallmark analysis", value = 0.2, {
      tryCatch({
        rv$msigdb_enrichment <- run_msigdb_hallmark_enrichment(
          rv$de,
          direction = input$msigdb_direction %||% "up",
          species = input$msigdb_species %||% rv$selected_organism_label,
          keytype = input$go_keytype %||% "SYMBOL",
          alpha = input$alpha,
          lfc_cutoff = input$lfc_cutoff,
          min_set_size = input$msigdb_min_set_size %||% 5
        )
        rv$msigdb_plot <- NULL
        append_log("MSigDB Hallmark analysis completed:", nrow(rv$msigdb_enrichment), "Hallmark sets tested.")
      }, error = function(e) {
        rv$msigdb_plot <- NULL
        showNotification(paste("MSigDB/Hallmark Error:", e$message), type = "error", duration = 15)
        append_log("MSigDB/Hallmark error:", e$message)
      })
    })
  })

  msigdb_display_data <- reactive({
    req(rv$msigdb_enrichment)
    d <- rv$msigdb_enrichment
    d[!is.na(d$pAdjusted) & d$pAdjusted <= (input$msigdb_pcut %||% 0.05), , drop = FALSE]
  })

  msigdb_plot_reactive <- reactive({
    d <- msigdb_display_data()
    req(nrow(d) > 0)
    make_msigdb_hallmark_plot(
      d,
      title = paste("MSigDB Hallmark:", input$msigdb_direction %||% "up"),
      top_n = input$msigdb_top_n %||% 20,
      point_alpha = input$plot_alpha,
      color_up = input$color_up %||% "#B2182B",
      color_down = input$color_down %||% "#2166AC",
      color_all = input$color_ns %||% "#B3B3B3",
      plot_theme = input$plot_theme %||% "classic",
      font_family = input$plot_font_family %||% "serif"
    )
  })

  output$msigdb_plot <- renderPlot({
    if (is.null(rv$msigdb_enrichment)) {
      plot.new(); text(0.5, 0.5, "Run Hallmark analysis first.", cex = 1.1)
    } else {
      d <- tryCatch(msigdb_display_data(), error = function(e) NULL)
      if (is.null(d) || nrow(d) == 0) {
        plot.new(); text(0.5, 0.5, "No Hallmark terms pass the current FDR cutoff.", cex = 1.1)
      } else {
        tryCatch(msigdb_plot_reactive(), error = function(e) { plot.new(); text(0.5, 0.5, e$message, cex = 1.1) })
      }
    }
  }, width = function() input$msigdb_plot_width, height = function() input$msigdb_plot_height)

  output$msigdb_table <- renderDT({
    d <- tryCatch(msigdb_display_data(), error = function(e) NULL)
    req(!is.null(d) && nrow(d) > 0)
    datatable(d, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE))
  })

  output$hallmark_gene_codes_ui <- renderUI({
    d <- tryCatch(msigdb_display_data(), error = function(e) NULL)
    if ((is.null(d) || nrow(d) == 0) && !is.null(rv$msigdb_enrichment)) {
      d <- rv$msigdb_enrichment
    }
    if (!is.null(d) && nrow(d) > 0 && "Hallmark" %in% names(d)) {
      ids <- stats::na.omit(as.character(d$Hallmark))
      ids <- unique(ids[nzchar(ids)])
      labels <- ids
      if ("Term" %in% names(d)) {
        terms_by_id <- stats::setNames(as.character(d$Term), as.character(d$Hallmark))
        labels <- paste0(ids, " - ", terms_by_id[ids])
        labels[is.na(labels) | !nzchar(labels)] <- ids[is.na(labels) | !nzchar(labels)]
      }
      return(selectizeInput(
        "hallmark_gene_codes",
        "Hallmark code(s)",
        choices = stats::setNames(ids, labels),
        selected = head(ids, 3),
        multiple = TRUE,
        options = list(placeholder = "Type to search Hallmark codes", plugins = list("remove_button"))
      ))
    }
    selectizeInput(
      "hallmark_gene_codes",
      "Hallmark code(s)",
      choices = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
      selected = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
      multiple = TRUE,
      options = list(create = TRUE, placeholder = "Example: HALLMARK_OXIDATIVE_PHOSPHORYLATION", plugins = list("remove_button"))
    )
  })

  observeEvent(input$run_hallmark_gene_lookup, {
    req(rv$de)
    validate(need(msigdb_species_available(input$msigdb_species), "Select an available MSigDB species."))
    append_log("Looking up Hallmark genes for", input$hallmark_gene_codes, level = "STEP")
    withProgress(message = "Looking up Hallmark genes", value = 0.2, {
      tryCatch({
        genes <- make_msigdb_hallmark_gene_table(
          rv$de,
          hallmark_ids = input$hallmark_gene_codes,
          species = input$msigdb_species %||% rv$selected_organism_label,
          keytype = input$go_keytype %||% "SYMBOL",
          alpha = input$alpha,
          lfc_cutoff = input$lfc_cutoff
        )
        rv$hallmark_genes <- genes
        if (nrow(genes) > 0) {
          rv$hallmark_gene_lookup <- stats::aggregate(
            gene_id ~ Hallmark + Hallmark_term,
            data = genes,
            FUN = function(x) length(unique(x))
          )
          names(rv$hallmark_gene_lookup)[names(rv$hallmark_gene_lookup) == "gene_id"] <- "Matched_loaded_genes"
        } else {
          rv$hallmark_gene_lookup <- data.frame()
        }
        append_log("Hallmark gene lookup completed:", nrow(rv$hallmark_genes), "matched gene rows.")
      }, error = function(e) {
        rv$hallmark_gene_lookup <- NULL
        rv$hallmark_genes <- NULL
        showNotification(paste("Hallmark gene lookup error:", e$message), type = "error", duration = 15)
        append_log("Hallmark gene lookup error:", e$message)
      })
    })
  })

  output$hallmark_genes_table <- renderDT({
    req(rv$hallmark_genes)
    d <- rv$hallmark_genes
    req(nrow(d) > 0)
    keep_cols <- c("gene_id", "original_gene_id", "Hallmark", "Hallmark_term", "Symbol", "log2FoldChange", "pValue", "padj", "DE_class", "short_description", "Protein_name")
    d <- d[, intersect(keep_cols, names(d)), drop = FALSE]
    d <- d[order(d$Hallmark, d$padj, d$pValue, na.last = TRUE), , drop = FALSE]
    d <- d[!duplicated(d[, intersect(c("gene_id", "Hallmark"), names(d)), drop = FALSE]), , drop = FALSE]
    button_data <- add_gene_count_buttons(d)
    datatable(button_data$data, rownames = FALSE, filter = "top", escape = !button_data$has_buttons,
              options = list(
                pageLength = 15,
                scrollX = TRUE,
                autoWidth = FALSE,
                columnDefs = if (button_data$has_buttons) gene_count_button_defs else list()
              ),
              callback = if (button_data$has_buttons) gene_count_button_callback else JS(""))
  })

  hallmark_genes_volcano_reactive <- reactive({
    req(rv$hallmark_genes)
    d <- rv$hallmark_genes
    req(nrow(d) > 0)
    plot_df <- d[order(d$padj, d$pValue, na.last = TRUE), , drop = FALSE]
    plot_df <- plot_df[!duplicated(plot_df$gene_id), , drop = FALSE]
    hallmark_label <- paste(unique(d$Hallmark), collapse = ", ")
    if (nchar(hallmark_label) > 80) hallmark_label <- paste0(substr(hallmark_label, 1, 77), "...")
    make_gene_group_volcano_plot(
      plot_df,
      paste0("Hallmark genes: ", hallmark_label),
      alpha = input$alpha,
      lfc_cutoff = input$lfc_cutoff,
      color_up = input$color_up %||% "#B2182B",
      color_down = input$color_down %||% "#2166AC",
      color_ns = input$color_ns %||% "#B3B3B3",
      plot_theme = input$plot_theme %||% "classic",
      font_family = input$plot_font_family %||% "serif"
    )
  })

  output$hallmark_genes_volcano <- renderPlot({
    if (is.null(rv$hallmark_gene_lookup)) {
      plot.new(); text(0.5, 0.5, "Enter Hallmark code(s) and click 'Find genes'.", cex = 1.1)
    } else if (is.null(rv$hallmark_genes) || nrow(rv$hallmark_genes) == 0) {
      plot.new(); text(0.5, 0.5, "No loaded genes matched the selected Hallmark set(s).", cex = 1.1)
    } else {
      tryCatch(hallmark_genes_volcano_reactive(), error = function(e) {
        plot.new(); text(0.5, 0.5, e$message, cex = 1.1)
      })
    }
  }, width = function() input$msigdb_plot_width, height = function() input$msigdb_plot_height)

  output$pmn_run_ui <- renderUI({
    selected_db <- trimws(input$pmn_cyc_db %||% "")
    if (!nzchar(selected_db)) {
      return(div(class = "muted", paste0(rv$selected_organism_label, " is not mapped to an available PMN Cyc database. Select or type a PMN Cyc DB to run.")))
    }
    actionButton("run_pmn", "Run PMN analysis", class = "btn-primary", style = "width:100%;")
  })

  observeEvent(input$pmn_cyc_db, {
    rv$selected_pmn_available <- nzchar(trimws(input$pmn_cyc_db %||% ""))
  }, ignoreInit = TRUE)

  observeEvent(input$run_pmn, {
    req(rv$de)
    validate(need(nzchar(trimws(input$pmn_cyc_db %||% "")), "Select a PMN Cyc database."))
    append_log("Running PMN enrichment for", input$pmn_direction, "genes in", input$pmn_cyc_db, level = "STEP")
    withProgress(message = "Running PMN pathway enrichment", value = 0.2, {
      tryCatch({
        rv$pmn_enrichment <- run_pmn_enrichment(
          rv$de,
          direction = input$pmn_direction %||% "up",
          cyc_db = input$pmn_cyc_db,
          alpha = input$alpha,
          lfc_cutoff = input$lfc_cutoff,
          min_set_size = input$pmn_min_set_size %||% 3
        )
        rv$pmn_plot <- NULL
        append_log("PMN enrichment completed:", nrow(rv$pmn_enrichment), "pathways tested.")
      }, error = function(e) {
        rv$pmn_plot <- NULL
        showNotification(paste("PMN Error:", e$message), type = "error", duration = 15)
        append_log("PMN error:", e$message)
      })
    })
  })

  pmn_display_data <- reactive({
    req(rv$pmn_enrichment)
    d <- rv$pmn_enrichment
    d[!is.na(d$p.adjusted) & d$p.adjusted <= (input$pmn_pcut %||% 0.05), , drop = FALSE]
  })

  output$pmn_pathway_codes_ui <- renderUI({
    d <- tryCatch(pmn_display_data(), error = function(e) NULL)
    if ((is.null(d) || nrow(d) == 0) && !is.null(rv$pmn_enrichment)) {
      d <- rv$pmn_enrichment
    }
    if (!is.null(d) && nrow(d) > 0 && "pathway.code" %in% names(d)) {
      ids <- stats::na.omit(as.character(d$pathway.code))
      ids <- ids[nzchar(ids)]
      ids <- unique(ids)
      labels <- ids
      if ("pathway.name" %in% names(d)) {
        names_by_id <- stats::setNames(as.character(d$pathway.name), as.character(d$pathway.code))
        labels <- paste0(ids, " - ", names_by_id[ids])
        labels[is.na(labels) | !nzchar(labels)] <- ids[is.na(labels) | !nzchar(labels)]
      }
      return(selectizeInput(
        "pmn_pathway_codes",
        "PMN pathway code(s)",
        choices = stats::setNames(ids, labels),
        selected = head(ids, 3),
        multiple = TRUE,
        options = list(placeholder = "Type to search PMN pathway codes", plugins = list("remove_button"))
      ))
    }
    selectizeInput(
      "pmn_pathway_codes",
      "PMN pathway code(s)",
      choices = c("PWY-6477"),
      selected = "PWY-6477",
      multiple = TRUE,
      options = list(create = TRUE, placeholder = "Example: PWY-6477", plugins = list("remove_button"))
    )
  })

  pmn_plot_reactive <- reactive({
    d <- pmn_display_data()
    req(nrow(d) > 0)
    plot_pmn_bubble(
      d,
      p_value_threshold = input$pmn_pcut %||% 0.05,
      top_n = input$pmn_top_n %||% 20,
      point_alpha = input$plot_alpha,
      color_up = input$color_up %||% "#B2182B",
      color_down = input$color_down %||% "#2166AC",
      color_all = input$color_ns %||% "#B3B3B3",
      plot_theme = input$plot_theme %||% "classic",
      font_family = input$plot_font_family %||% "serif"
    )
  })

  output$pmn_plot <- renderPlot({
    if (is.null(rv$pmn_enrichment)) {
      plot.new(); text(0.5, 0.5, "Run PMN analysis first.", cex = 1.1)
    } else {
      d <- tryCatch(pmn_display_data(), error = function(e) NULL)
      if (is.null(d) || nrow(d) == 0) {
        plot.new(); text(0.5, 0.5, "No PMN pathways pass the current FDR cutoff.", cex = 1.1)
      } else {
        tryCatch(pmn_plot_reactive(), error = function(e) { plot.new(); text(0.5, 0.5, e$message, cex = 1.1) })
      }
    }
  }, width = function() input$pmn_plot_width, height = function() input$pmn_plot_height)

  output$pmn_table <- renderDT({
    d <- tryCatch(pmn_display_data(), error = function(e) NULL)
    req(!is.null(d) && nrow(d) > 0)
    datatable(d, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE))
  })

  observeEvent(input$run_pmn_pathway_lookup, {
    req(rv$de)
    validate(need(nzchar(trimws(input$pmn_cyc_db %||% "")), "Select a PMN Cyc database."))
    append_log("Looking up PMN pathway genes for", input$pmn_pathway_codes, "in", input$pmn_cyc_db, level = "STEP")
    withProgress(message = "Looking up PMN pathway genes", value = 0.2, {
      tryCatch({
        pmn_data <- get_pmn_genes(cyc_db = input$pmn_cyc_db)
        rv$pmn_pathway_lookup <- make_pmn_pathway_summary(
          rv$de,
          pathway_codes = input$pmn_pathway_codes,
          cyc_db = input$pmn_cyc_db,
          alpha = input$alpha,
          lfc_cutoff = input$lfc_cutoff,
          pmn_data = pmn_data
        )
        rv$pmn_pathway_genes <- make_pmn_pathway_gene_table(
          rv$de,
          pathway_codes = input$pmn_pathway_codes,
          cyc_db = input$pmn_cyc_db,
          alpha = input$alpha,
          lfc_cutoff = input$lfc_cutoff,
          pmn_data = pmn_data
        )
        append_log("PMN pathway lookup completed:", nrow(rv$pmn_pathway_lookup), "pathway rows and", nrow(rv$pmn_pathway_genes), "matched genes.")
      }, error = function(e) {
        rv$pmn_pathway_lookup <- NULL
        rv$pmn_pathway_genes <- NULL
        showNotification(paste("PMN pathway lookup error:", e$message), type = "error", duration = 15)
        append_log("PMN pathway lookup error:", e$message)
      })
    })
  })

  output$pmn_pathway_table <- renderDT({
    req(rv$pmn_pathway_genes)
    d <- rv$pmn_pathway_genes
    req(nrow(d) > 0)
    keep_cols <- c("gene_id", "pathway.code", "pathway.name", "Symbol", "log2FoldChange", "pValue", "padj", "DE_class")
    d <- d[, intersect(keep_cols, names(d)), drop = FALSE]
    d <- d[order(d$pathway.code, d$padj, d$pValue, na.last = TRUE), , drop = FALSE]
    d <- d[!duplicated(d[, intersect(c("gene_id", "pathway.code"), names(d)), drop = FALSE]), , drop = FALSE]
    button_data <- add_gene_count_buttons(d)
    datatable(button_data$data, rownames = FALSE, filter = "top", escape = !button_data$has_buttons,
              options = list(
                pageLength = 15,
                scrollX = TRUE,
                autoWidth = FALSE,
                columnDefs = if (button_data$has_buttons) gene_count_button_defs else list()
              ),
              callback = if (button_data$has_buttons) gene_count_button_callback else JS(""))
  })

  pmn_pathway_volcano_reactive <- reactive({
    req(rv$pmn_pathway_genes)
    d <- rv$pmn_pathway_genes
    req(nrow(d) > 0)
    plot_df <- d[order(d$padj, d$pValue, na.last = TRUE), , drop = FALSE]
    plot_df <- plot_df[!duplicated(plot_df$gene_id), , drop = FALSE]
    pathway_label <- paste(unique(d$pathway.code), collapse = ", ")
    if (nchar(pathway_label) > 80) pathway_label <- paste0(substr(pathway_label, 1, 77), "...")
    make_gene_group_volcano_plot(
      plot_df,
      paste0("PMN pathway genes: ", pathway_label),
      alpha = input$alpha,
      lfc_cutoff = input$lfc_cutoff,
      color_up = input$color_up %||% "#B2182B",
      color_down = input$color_down %||% "#2166AC",
      color_ns = input$color_ns %||% "#B3B3B3",
      plot_theme = input$plot_theme %||% "classic",
      font_family = input$plot_font_family %||% "serif"
    )
  })

  output$pmn_pathway_volcano <- renderPlot({
    if (is.null(rv$pmn_pathway_lookup)) {
      plot.new(); text(0.5, 0.5, "Enter PMN pathway code(s) and click 'Find genes'.", cex = 1.1)
    } else if (is.null(rv$pmn_pathway_genes) || nrow(rv$pmn_pathway_genes) == 0) {
      plot.new(); text(0.5, 0.5, "No loaded genes matched the selected PMN pathway(s).", cex = 1.1)
    } else {
      tryCatch(pmn_pathway_volcano_reactive(), error = function(e) {
        plot.new(); text(0.5, 0.5, e$message, cex = 1.1)
      })
    }
  }, width = function() input$pmn_plot_width, height = function() input$pmn_plot_height)

  observeEvent(input$run_offspring, {
    req(rv$de)
    append_log("Creating GO offspring summary", level = "STEP")
    tryCatch({
      rv$offspring <- make_go_offspring_summary(
        rv$de, input$parent_go_ids, input$alpha, input$lfc_cutoff,
        orgdb = input$go_orgdb %||% "org.At.tair.db",
        keytype = input$go_keytype %||% "TAIR"
      )
      append_log("GO offspring summary created:", nrow(rv$offspring), "parent terms.")
    }, error = function(e) {
      showNotification(e$message, type = "error", duration = 12)
      append_log("GO offspring error:", e$message)
    })
  })
  output$offspring_table <- renderDT({
    req(rv$offspring)
    datatable(rv$offspring, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  observeEvent(input$run_stress, {
    req(rv$de)
    append_log("Running abiotic stress enrichment for", input$stress_dataset, "genes", level = "STEP")
    tryCatch({
      rv$stress <- make_abiotic_stress_table(
        rv$de, dataset = input$stress_dataset, alpha = input$alpha, lfc_cutoff = input$lfc_cutoff,
        orgdb = input$go_orgdb %||% "org.At.tair.db",
        keytype = input$go_keytype %||% "TAIR"
      )
      rv$stress_plot <- NULL
      append_log("Abiotic stress enrichment finished for", input$stress_dataset, "gene set.")
    }, error = function(e) {
      showNotification(e$message, type = "error", duration = 12)
      append_log("Stress enrichment error:", e$message)
    })
  })
  stress_plot_reactive <- reactive({
    req(rv$stress)
    make_abiotic_stress_plot(
      rv$stress,
      paste("Abiotic stress enrichment:", input$stress_dataset),
      plot_theme = input$plot_theme %||% "classic",
      font_family = input$plot_font_family %||% "serif"
    )
  })

  output$stress_plot <- renderPlot({
    if (is.null(rv$stress)) {
      plot.new(); text(0.5, 0.5, "Run abiotic stress enrichment first.", cex = 1.1)
    } else {
      stress_plot_reactive()
    }
  }, width = function() input$go_plot_width, height = function() input$go_plot_height)
  output$stress_table <- renderDT({
    req(rv$stress)
    datatable(rv$stress, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  observeEvent(input$run_revigo, {
    req(rv$de)
    append_log("Running REVIGO-like GO reduction", level = "STEP")
    withProgress(message = "Running REVIGO-like semantic reduction", value = 0.1, {
      tryCatch({
        direction <- input$revigo_direction %||% "up"
        direction_label <- switch(direction,
          up = "Upregulated genes",
          down = "Downregulated genes",
          all = "All significant genes",
          direction
        )
        incProgress(0.3, detail = paste("Computing GO enrichment for", direction))
        go_terms <- compute_go_cached(direction)
        if (is.null(go_terms) || nrow(go_terms) < 2) {
          rv$revigo <- NULL
          rv$revigo_direction <- direction
          stop(paste("Not enough", direction, "GO terms to reduce."))
        }
        incProgress(0.4, detail = "Reducing GO terms")
        rv$revigo <- run_rrvgo_reduce(go_terms, ontology = input$ontology, top_n = input$revigo_top_n,
                                      threshold = input$revigo_threshold, title = direction_label,
                                      orgdb = input$go_orgdb %||% "org.At.tair.db",
                                      plot_theme = input$plot_theme %||% "classic",
                                      font_family = input$plot_font_family %||% "serif",
                                      algorithm = input$revigo_algorithm %||% "umap",
                                      color_palette = input$revigo_color_palette %||% "set1",
                                      show_labels = isTRUE(input$revigo_show_labels),
                                      max_parent_labels = input$revigo_parent_label_count %||% 8)
        rv$revigo_direction <- direction
        append_log("REVIGO-like analysis finished for", direction, "genes.")
      }, error = function(e) {
        showNotification(e$message, type = "error", duration = 12)
        append_log("REVIGO-like error:", e$message)
      })
    })
  })

  output$revigo_title <- renderText({
    if (is.null(rv$revigo_direction)) {
      "REVIGO-like GO terms"
    } else {
      paste("REVIGO-like GO terms:", rv$revigo_direction)
    }
  })

  revigo_plot_reactive <- reactive({
    req(rv$revigo)
    if (is.null(rv$revigo$plot_data)) return(rv$revigo$plot)
    direction <- rv$revigo_direction %||% "selected"
    direction_label <- switch(direction,
      up = "Upregulated genes",
      down = "Downregulated genes",
      all = "All significant genes",
      direction
    )
    plot_revigo_data(
      rv$revigo$plot_data,
      title = direction_label,
      color_palette = input$revigo_color_palette %||% "set1",
      show_labels = isTRUE(input$revigo_show_labels),
      max_parent_labels = input$revigo_parent_label_count %||% 8,
      plot_theme = input$plot_theme %||% "classic",
      font_family = input$plot_font_family %||% "serif"
    )$plot
  })

  revigo_treemap_reactive <- reactive({
    req(rv$revigo, rv$revigo$table)
    direction <- rv$revigo_direction %||% "selected"
    direction_label <- switch(direction,
      up = "Upregulated genes",
      down = "Downregulated genes",
      all = "All significant genes",
      direction
    )
    make_revigo_treemap_plot(rv$revigo$table, title = paste("REVIGO-like treemap:", direction_label))
  })

  output$revigo_plot <- renderPlot({
    if (is.null(rv$revigo)) {
      plot.new(); text(0.5, 0.5, "Select a gene set and run REVIGO-like analysis.", cex = 1.1)
    } else {
      revigo_plot_reactive()
    }
  }, width = function() input$go_plot_width, height = function() input$go_plot_height)

  output$revigo_treemap <- renderPlot({
    if (is.null(rv$revigo)) {
      plot.new(); text(0.5, 0.5, "Select a gene set and run REVIGO-like analysis.", cex = 1.1)
    } else {
      revigo_treemap_reactive()
    }
  }, width = function() input$go_plot_width, height = function() input$go_plot_height)
  # ---- TE Analysis -------------------------------------------------------
  observeEvent(input$run_te_enrich, {
    req(rv$de)
    append_log("Running TE superfamily enrichment", level = "STEP")
    withProgress(message = "Running TE Enrichment...", value = 0.3, {
      tryCatch({
        res <- run_te_enrichment(rv$de, pvalue_cutoff = input$te_enrich_pvalue)
        rv$te_enrichment <- res
        
        incProgress(0.5, detail = "Generating Plot")
        rv$te_enrich_bubble <- NULL
        
        append_log("TE enrichment completed.")
      }, error = function(e) {
        showNotification(paste("TE Enrichment Error:", e$message), type = "error", duration = 15)
        append_log("TE Enrichment error:", e$message)
      })
    })
  })
  
  output$te_enrich_table <- renderDT({
    req(rv$te_enrichment)
    datatable(rv$te_enrichment, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE))
  })
  
  te_enrich_bubble_reactive <- reactive({
    req(rv$te_enrichment)
    plot_te_enrichment(
      rv$te_enrichment,
      p_value_threshold = input$te_enrich_pvalue,
      color_up = input$color_up %||% "#B2182B",
      color_down = input$color_down %||% "#2166AC",
      plot_theme = input$plot_theme %||% "classic",
      font_family = input$plot_font_family %||% "serif"
    )
  })
  
  output$te_enrich_bubble <- renderPlot({
    p <- tryCatch(te_enrich_bubble_reactive(), error = function(e) NULL)
    if (is.null(p)) {
      plot.new(); text(0.5, 0.5, "No significantly enriched TE superfamilies.", cex = 1.1)
    } else {
      p
    }
  }, width = function() input$teg_plot_width, height = function() input$teg_plot_height)
  
  observeEvent(input$run_te_volcano, {
    req(rv$de)
    req(length(input$te_super_families) > 0)
    append_log("Building TE volcano for selected superfamilies", level = "STEP")
    withProgress(message = "Building TE volcano", value = 0.2, {
      tryCatch({
        out <- make_retro_te_volcano(
          rv$de,
          super_families = input$te_super_families,
          padj_cutoff = input$te_padj_cutoff,
          lfc_cutoff = input$te_lfc_cutoff,
          point_size = input$plot_point_size,
          point_alpha = input$plot_alpha,
          plot_theme = input$plot_theme %||% "classic",
          font_family = input$plot_font_family %||% "serif"
        )
        rv$te_volcano <- out$data
        rv$te_volcano_plot <- out$plot
        append_log("TEG volcano built:", nrow(out$data), "rows.")
      }, error = function(e) {
        showNotification(e$message, type = "error", duration = 12)
        append_log("TEG volcano error:", e$message)
      })
    })
  })
  te_volcano_plot_reactive <- reactive({
    req(rv$te_volcano)
    plot_retro_te_volcano_data(
      rv$te_volcano,
      color_col = "Transposon_Super_Family",
      label_name = "Transposon\nSuper-Family",
      padj_cutoff = input$te_padj_cutoff,
      lfc_cutoff = input$te_lfc_cutoff,
      point_size = input$plot_point_size,
      point_alpha = input$plot_alpha,
      plot_theme = input$plot_theme %||% "classic",
      font_family = input$plot_font_family %||% "serif"
    )
  })

  output$te_volcano_plot <- renderPlot({
    if (is.null(rv$te_volcano_plot)) {
      plot.new(); text(0.5, 0.5, "Run TE volcano first.", cex = 1.1)
    } else {
      te_volcano_plot_reactive()
    }
  }, width = function() input$teg_plot_width, height = function() input$teg_plot_height)
  output$te_volcano_table <- renderDT({
    req(rv$te_volcano)
    button_data <- add_gene_count_buttons(rv$te_volcano)
    datatable(button_data$data, rownames = FALSE, filter = "top", escape = !button_data$has_buttons,
              options = list(
                pageLength = 12,
                scrollX = TRUE,
                autoWidth = FALSE,
                columnDefs = if (button_data$has_buttons) gene_count_button_defs else list()
              ),
              callback = if (button_data$has_buttons) gene_count_button_callback else JS(""))
  })

  output$run_log <- renderText({ paste(rv$log, collapse = "\n") })

  output$download_de <- downloadHandler(
    filename = function() paste0("DE_results_", Sys.Date(), ".csv"),
    content = function(file) write.csv(rv$de, file, row.names = FALSE)
  )

  output$download_all_deseq_comparisons <- downloadHandler(
    filename = function() paste0("DESeq2_all_vs_", safe_filename_part(rv$deseq_all_comparisons$control %||% "control"), "_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$deseq_all_comparisons, rv$deseq_all_comparisons$tables)
      tables <- Map(function(tbl, comparison) {
        tbl$Comparison <- comparison
        tbl$Control <- rv$deseq_all_comparisons$control
        tbl[, c("Comparison", "Control", setdiff(names(tbl), c("Comparison", "Control"))), drop = FALSE]
      }, rv$deseq_all_comparisons$tables, names(rv$deseq_all_comparisons$tables))
      write.csv(do.call(rbind, tables), file, row.names = FALSE)
    }
  )

  output$download_annotation_table <- downloadHandler(
    filename = function() {
      org <- safe_filename_part(rv$selected_organism_label)
      paste0("annotations_", org, "_taxid_", rv$selected_tax_id, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(rv$annotation_df)
      d <- rv$annotation_df
      d$annotation_key <- NULL
      write.csv(d, file, row.names = FALSE)
    }
  )

  output$download_annotated_de <- downloadHandler(
    filename = function() paste0("annotated_DE_results_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$de)
      write.csv(rv$de, file, row.names = FALSE)
    }
  )
  
  output$download_search_annotations <- downloadHandler(
    filename = function() paste0("Search_Annotations_", Sys.Date(), ".csv"),
    content = function(file) { req(rv$de); write.csv(rv$de, file, row.names = FALSE) }
  )
  
  output$download_norm_counts <- downloadHandler(
    filename = function() paste0("normalized_counts_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$norm_counts)
      write.csv(rv$norm_counts, file, row.names = FALSE)
    }
  )
  output$download_volcano <- download_plot_server(volcano_reactive, reactive(input$format_volcano), "volcano", reactive(input$de_plot_width), reactive(input$de_plot_height))
  output$download_ma <- download_plot_server(ma_reactive, reactive(input$format_ma), "MA_plot", reactive(input$de_plot_width), reactive(input$de_plot_height))
  output$download_pca <- download_plot_server(pca_reactive, reactive(input$format_pca), "PCA", reactive(input$de_plot_width), reactive(input$de_plot_height))
  output$download_pca_table <- downloadHandler(
    filename = function() paste0("PCA_table_", Sys.Date(), ".csv"),
    content = function(file) { req(rv$pca); write.csv(rv$pca, file, row.names = FALSE) }
  )
  output$download_go_bubble <- download_plot_server(go_display_plot, reactive(input$format_go_bubble), "GO_bubble", reactive(input$go_plot_width), reactive(input$go_plot_height))
  output$download_go_table <- downloadHandler(
    filename = function() paste0("GO_table_", input$go_direction, "_", input$ontology, "_", Sys.Date(), ".csv"),
    content = function(file) {
      d <- tryCatch(go_display_data(), error = function(e) NULL)
      req(!is.null(d)); write.csv(d, file, row.names = FALSE)
    }
  )
  output$download_go_genes_volcano <- download_plot_server(
    go_genes_volcano_reactive,
    reactive(input$format_go_genes_volcano),
    reactive(paste0("GO_genes_volcano_", gsub("[^A-Za-z0-9]+", "_", paste(input$go_gene_codes %||% "GO", collapse = "_")))),
    reactive(input$go_plot_width),
    reactive(input$go_plot_height)
  )
  output$download_go_genes_table <- downloadHandler(
    filename = function() paste0("GO_genes_", input$ontology %||% "BP", "_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$go_genes)
      write.csv(rv$go_genes, file, row.names = FALSE)
    }
  )
  output$download_msigdb <- download_plot_server(msigdb_plot_reactive, reactive(input$format_msigdb), "MSigDB_Hallmark", reactive(input$msigdb_plot_width), reactive(input$msigdb_plot_height))
  output$download_msigdb_table <- downloadHandler(
    filename = function() paste0("MSigDB_Hallmark_", input$msigdb_direction, "_", Sys.Date(), ".csv"),
    content = function(file) {
      d <- tryCatch(msigdb_display_data(), error = function(e) NULL)
      req(!is.null(d)); write.csv(d, file, row.names = FALSE)
    }
  )
  output$download_hallmark_genes_volcano <- download_plot_server(
    hallmark_genes_volcano_reactive,
    reactive(input$format_hallmark_genes_volcano),
    reactive(paste0("Hallmark_genes_volcano_", gsub("[^A-Za-z0-9]+", "_", paste(input$hallmark_gene_codes %||% "Hallmark", collapse = "_")))),
    reactive(input$msigdb_plot_width),
    reactive(input$msigdb_plot_height)
  )
  output$download_hallmark_genes_table <- downloadHandler(
    filename = function() paste0("Hallmark_genes_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$hallmark_genes)
      write.csv(rv$hallmark_genes, file, row.names = FALSE)
    }
  )
  output$download_pmn <- download_plot_server(
    pmn_plot_reactive,
    reactive(input$format_pmn),
    reactive(paste0("PMN_", gsub("[^A-Za-z0-9]+", "_", input$pmn_cyc_db %||% "Cyc"))),
    reactive(input$pmn_plot_width),
    reactive(input$pmn_plot_height)
  )
  output$download_pmn_table <- downloadHandler(
    filename = function() paste0("PMN_", gsub("[^A-Za-z0-9]+", "_", input$pmn_cyc_db %||% "Cyc"), "_", input$pmn_direction, "_", Sys.Date(), ".csv"),
    content = function(file) {
      d <- tryCatch(pmn_display_data(), error = function(e) NULL)
      req(!is.null(d)); write.csv(d, file, row.names = FALSE)
    }
  )
  output$download_pmn_pathway_volcano <- download_plot_server(
    pmn_pathway_volcano_reactive,
    reactive(input$format_pmn_pathway_volcano),
    reactive(paste0("PMN_pathway_volcano_", gsub("[^A-Za-z0-9]+", "_", input$pmn_cyc_db %||% "Cyc"))),
    reactive(input$pmn_plot_width),
    reactive(input$pmn_plot_height)
  )
  output$download_pmn_pathway_table <- downloadHandler(
    filename = function() paste0("PMN_pathway_genes_", gsub("[^A-Za-z0-9]+", "_", input$pmn_cyc_db %||% "Cyc"), "_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$pmn_pathway_genes)
      write.csv(rv$pmn_pathway_genes, file, row.names = FALSE)
    }
  )
  output$download_offspring <- downloadHandler(
    filename = function() paste0("GO_offspring_summary_", Sys.Date(), ".csv"),
    content = function(file) { req(rv$offspring); write.csv(rv$offspring, file, row.names = FALSE) }
  )
  output$download_stress <- download_plot_server(stress_plot_reactive, reactive(input$format_stress), "abiotic_stress", reactive(input$go_plot_width), reactive(input$go_plot_height))
  output$download_stress_table <- downloadHandler(
    filename = function() paste0("abiotic_stress_enrichment_", input$stress_dataset, "_", Sys.Date(), ".csv"),
    content = function(file) { req(rv$stress); write.csv(rv$stress, file, row.names = FALSE) }
  )
  output$download_revigo <- download_plot_server(
    revigo_plot_reactive,
    reactive(input$format_revigo),
    reactive(paste0("REVIGO_", rv$revigo_direction %||% "selected")),
    reactive(input$go_plot_width),
    reactive(input$go_plot_height)
  )
  output$download_revigo_treemap <- download_plot_server(
    revigo_treemap_reactive,
    reactive(input$format_revigo_treemap),
    reactive(paste0("REVIGO_treemap_", rv$revigo_direction %||% "selected")),
    reactive(input$go_plot_width),
    reactive(input$go_plot_height)
  )
  output$download_revigo_table <- downloadHandler(
    filename = function() paste0("REVIGO_", rv$revigo_direction %||% "selected", "_table_", input$ontology, "_", Sys.Date(), ".csv"),
    content = function(file) { req(rv$revigo); write.csv(rv$revigo$table, file, row.names = FALSE) }
  )
  output$download_te_enrich_bubble <- download_plot_server(te_enrich_bubble_reactive, reactive(input$format_te_enrich_bubble), "TE_Enrichment_bubble", reactive(input$teg_plot_width), reactive(input$teg_plot_height))
  output$download_te_enrich_table <- downloadHandler(
    filename = function() paste0("TE_Enrichment_table_", Sys.Date(), ".csv"),
    content = function(file) { req(rv$te_enrichment); write.csv(rv$te_enrichment, file, row.names = FALSE) }
  )
  output$download_te_volcano <- download_plot_server(te_volcano_plot_reactive, reactive(input$format_te_volcano), "TEG_volcano", reactive(input$teg_plot_width), reactive(input$teg_plot_height))
  output$download_te_volcano_table <- downloadHandler(
    filename = function() paste0("TEG_volcano_table_", Sys.Date(), ".csv"),
    content = function(file) { req(rv$te_volcano); write.csv(rv$te_volcano, file, row.names = FALSE) }
  )

  # ---- Gene Groups ---------------------------------------------------------
  output$gene_group_run_ui <- renderUI({
    tax_id <- suppressWarnings(as.integer(rv$selected_tax_id))
    if (!identical(tax_id, 3702L)) {
      return(div(class = "muted", "Available only for Arabidopsis thaliana."))
    }
    actionButton("run_gene_groups", "Build gene groups", class = "btn-primary", style = "width:100%;")
  })

  # Phase 1: "Build" button → download reference data + GO terms (one-time, slow)
  observeEvent(input$run_gene_groups, {
    req(rv$de)
    tax_id <- suppressWarnings(as.integer(rv$selected_tax_id))
    req(identical(tax_id, 3702L))
    append_log("Loading gene group reference data", level = "STEP")
    withProgress(message = "Loading reference data from GitHub + GO terms...", value = 0.1, {
      tryCatch({
        incProgress(0.8, detail = "Downloading gene-set lists and GO offspring terms")
        ctx <- setup_gene_groups(rv$de, alpha = input$alpha, lfc_cutoff = input$lfc_cutoff)
        rv$gg_context   <- ctx
        rv$gg_cache     <- list()   # clear per-group cache when re-built
        updateSelectInput(session, "gene_group_select",
                          choices  = stats::setNames(GENE_GROUP_NAMES, gsub("_", " ", GENE_GROUP_NAMES, fixed = TRUE)),
                          selected = GENE_GROUP_NAMES[1])
        append_log("Gene group reference data loaded. Select a group to compute it.")
      }, error = function(e) {
        showNotification(paste("Gene groups setup error:", e$message), type = "error", duration = 15)
        append_log("Gene groups setup error:", e$message)
      })
    })
  })

  # Phase 2: group selection → compute only the chosen group (fast, cached)
  observeEvent(input$gene_group_select, {
    req(rv$gg_context, input$gene_group_select)
    gn <- input$gene_group_select
    if (!is.null(rv$gg_cache[[gn]])) return()  # already computed
    append_log("Computing gene group:", gn, level = "STEP")
    withProgress(message = paste("Computing:", gsub("_", " ", gn)), value = 0.5, {
      tryCatch({
        result <- compute_single_group(gn, rv$gg_context,
                                       color_up = input$color_up %||% "#B2182B", color_down = input$color_down %||% "#2166AC",
                                       color_ns = input$color_ns %||% "#B3B3B3",
                                       plot_theme = input$plot_theme %||% "classic",
                                       font_family = input$plot_font_family %||% "serif")
        rv$gg_cache[[gn]] <- result
        append_log("Gene group computed:", gn,
                   if (!is.null(result$data)) paste0("(", nrow(result$data), " genes)") else "(0 genes)")
      }, error = function(e) {
        rv$gg_cache[[gn]] <- list(data = NULL, plot = NULL)
        showNotification(paste("Group error:", e$message), type = "warning", duration = 10)
        append_log("Gene group error [", gn, "]:", e$message)
      })
    })
  })

  output$gene_group_selector_ui <- renderUI({
    req(rv$gg_context)
    group_labels <- gsub("_", " ", GENE_GROUP_NAMES, fixed = TRUE)
    selectInput("gene_group_select", "Select gene group",
                choices  = stats::setNames(GENE_GROUP_NAMES, group_labels),
                selected = GENE_GROUP_NAMES[1])
  })

  gene_group_plot_reactive <- reactive({
    req(rv$gg_cache, input$gene_group_select)
    cached <- rv$gg_cache[[input$gene_group_select]]
    req(!is.null(cached), !is.null(cached$data))
    make_gene_group_volcano_plot(
      cached$data,
      input$gene_group_select,
      alpha = input$alpha,
      lfc_cutoff = input$lfc_cutoff,
      color_up = input$color_up %||% "#B2182B",
      color_down = input$color_down %||% "#2166AC",
      color_ns = input$color_ns %||% "#B3B3B3",
      plot_theme = input$plot_theme %||% "classic",
      font_family = input$plot_font_family %||% "serif"
    )
  })

  output$gene_group_volcano <- renderPlot({
    if (is.null(rv$gg_context)) {
      plot.new(); text(0.5, 0.5, "Click 'Build gene groups' first.", cex = 1.1)
    } else if (is.null(rv$gg_cache[[input$gene_group_select]])) {
      plot.new(); text(0.5, 0.5, "Computing...", cex = 1.1)
    } else {
      p <- tryCatch(gene_group_plot_reactive(), error = function(e) NULL)
      if (is.null(p)) { plot.new(); text(0.5, 0.5, "No genes in this group.", cex = 1.1) }
      else p
    }
  }, width = function() input$grp_plot_width, height = function() input$grp_plot_height)

  output$gene_group_table <- renderDT({
    req(rv$gg_cache, input$gene_group_select)
    d <- rv$gg_cache[[input$gene_group_select]]$data
    req(!is.null(d) && nrow(d) > 0)
    keep_cols <- c("gene_id", "Symbol", "log2FoldChange", "padj")
    d <- d[, intersect(keep_cols, names(d)), drop = FALSE]
    button_data <- add_gene_count_buttons(d)
    datatable(button_data$data, rownames = FALSE, filter = "top", escape = !button_data$has_buttons,
              options = list(
                pageLength = 15,
                scrollX = TRUE,
                autoWidth = FALSE,
                columnDefs = if (button_data$has_buttons) gene_count_button_defs else list()
              ),
              callback = if (button_data$has_buttons) gene_count_button_callback else JS(""))
  })

  output$download_gene_group_volcano <- download_plot_server(
    gene_group_plot_reactive,
    reactive(input$format_gene_group_volcano),
    reactive(paste0("gene_group_", gsub(" ", "_", input$gene_group_select %||% "group"))),
    reactive(input$grp_plot_width),
    reactive(input$grp_plot_height)
  )

  output$download_gene_group_table <- downloadHandler(
    filename = function() paste0("gene_group_", input$gene_group_select, "_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(rv$gg_cache, input$gene_group_select)
      d <- rv$gg_cache[[input$gene_group_select]]$data
      req(!is.null(d))
      write.csv(d, file, row.names = FALSE)
    }
  )

  gene_family_available <- reactive({
    suppressWarnings(as.integer(rv$selected_tax_id)) %in% c(3702L, 9606L)
  })

  gene_family_backend_label <- reactive({
    tax_id <- suppressWarnings(as.integer(rv$selected_tax_id))
    if (identical(tax_id, 3702L)) return("Arabidopsis TAIR gene-family file")
    if (identical(tax_id, 9606L)) return("human HGNC gene-family database")
    paste0(rv$selected_organism_label %||% "Selected organism", " is not available for gene-family analysis.")
  })

  output$gene_family_status_ui <- renderUI({
    if (isTRUE(gene_family_available())) {
      div(class = "muted", style = "margin-top: 8px;",
          paste0("Uses the ", gene_family_backend_label(), ". The selected organism is ",
                 rv$selected_organism_label, " (", rv$selected_tax_id, ")."))
    } else {
      div(class = "muted", style = "margin-top: 8px;",
          "Gene-family analysis is available only for Arabidopsis thaliana and Homo sapiens. Select one in Organism annotations.")
    }
  })

  output$gene_family_enrichment_status_ui <- renderUI({
    if (isTRUE(gene_family_available())) {
      div(class = "muted", paste0("Uses the ", gene_family_backend_label(), ". If families are not loaded yet, the app loads them when enrichment starts."))
    } else {
      div(class = "muted", "Gene-family enrichment is available only for Arabidopsis thaliana and Homo sapiens.")
    }
  })

  output$gene_family_load_ui <- renderUI({
    if (!isTRUE(gene_family_available())) return(div(class = "muted", "Not available"))
    actionButton("run_gene_families", "Load gene families", class = "btn-primary", style = "width:100%;")
  })

  output$gene_family_enrichment_run_ui <- renderUI({
    if (!isTRUE(gene_family_available())) return(div(class = "muted", "Not available for selected organism."))
    actionButton("run_gene_family_enrichment", "Run family enrichment", class = "btn-primary", style = "width:100%;")
  })

  observeEvent(input$run_gene_families, {
    req(rv$de, gene_family_available())
    append_log("Loading gene family reference data for", rv$selected_organism_label, level = "STEP")
    withProgress(message = "Loading gene family reference data from GitHub...", value = 0.4, {
      tryCatch({
        ctx <- setup_gene_family_analysis(
          rv$de,
          tax_id = rv$selected_tax_id,
          alpha = input$alpha,
          lfc_cutoff = input$lfc_cutoff,
          gene_id_type = input$go_keytype %||% "TAIR",
          orgdb = input$go_orgdb %||% "org.At.tair.db"
        )
        rv$gene_family_context <- ctx
        rv$gene_family_cache <- list()
        rv$gene_family_enrichment <- NULL
        append_log("Gene family reference data loaded:", length(ctx$choices), "families.")
      }, error = function(e) {
        showNotification(paste("Gene family setup error:", e$message), type = "error", duration = 15)
        append_log("Gene family setup error:", e$message)
      })
    })
  })

  output$gene_family_selector_ui <- renderUI({
    req(rv$gene_family_context)
    choices <- rv$gene_family_context$choices
    selected <- head(choices, 1)
    if (!is.null(rv$gene_family_enrichment) && nrow(rv$gene_family_enrichment) > 0 &&
        "Gene_Family" %in% names(rv$gene_family_enrichment)) {
      top_families <- head(unique(as.character(rv$gene_family_enrichment$Gene_Family)), 3)
      top_families <- intersect(top_families, choices)
      if (length(top_families) > 0) {
        choices <- c(top_families, setdiff(choices, top_families))
        selected <- top_families
      }
    }
    selectizeInput(
      "gene_family_select", "Gene families",
      choices = choices,
      selected = selected,
      multiple = TRUE,
      options = list(placeholder = "Type to search gene families", plugins = list("remove_button"))
    )
  })

  gene_family_data_reactive <- reactive({
    req(rv$gene_family_context, input$gene_family_select)
    compute_at_gene_family(input$gene_family_select, rv$gene_family_context)
  })

  gene_family_plot_reactive <- reactive({
    d <- gene_family_data_reactive()
    family_title <- paste(input$gene_family_select %||% "gene family", collapse = ", ")
    make_at_gene_family_volcano_plot(
      d,
      family_title,
      alpha = input$alpha,
      lfc_cutoff = input$lfc_cutoff,
      color_up = input$color_up %||% "#B2182B",
      color_down = input$color_down %||% "#2166AC",
      color_ns = input$color_ns %||% "#B3B3B3",
      plot_theme = input$plot_theme %||% "classic",
      font_family = input$plot_font_family %||% "serif"
    )
  })

  output$gene_family_volcano <- renderPlot({
    if (is.null(rv$gene_family_context)) {
      plot.new(); text(0.5, 0.5, "Click 'Load gene families' first.", cex = 1.1)
    } else {
      p <- tryCatch(gene_family_plot_reactive(), error = function(e) NULL)
      if (is.null(p)) { plot.new(); text(0.5, 0.5, "No genes in the selected family.", cex = 1.1) }
      else p
    }
  }, width = function() input$grp_plot_width, height = function() input$grp_plot_height)

  output$gene_family_table <- renderDT({
    d <- gene_family_data_reactive()
    req(!is.null(d) && nrow(d) > 0)
    keep_cols <- c("gene_id", "hgnc_id", "Gene_Family", "Sub_Family", "Symbol", "log2FoldChange", "pValue", "padj", "DE_class", "short_description")
    d <- d[, intersect(keep_cols, names(d)), drop = FALSE]
    if ("gene_id" %in% names(d)) d <- d[!duplicated(d$gene_id), , drop = FALSE]
    button_data <- add_gene_count_buttons(d)
    datatable(button_data$data, rownames = FALSE, filter = "top", escape = !button_data$has_buttons,
              options = list(
                pageLength = 15,
                scrollX = TRUE,
                autoWidth = FALSE,
                columnDefs = if (button_data$has_buttons) gene_count_button_defs else list()
              ),
              callback = if (button_data$has_buttons) gene_count_button_callback else JS(""))
  })

  output$download_gene_family_volcano <- download_plot_server(
    gene_family_plot_reactive,
    reactive(input$format_gene_family_volcano),
    reactive(paste0("gene_family_", gsub("[^A-Za-z0-9]+", "_", paste(input$gene_family_select %||% "family", collapse = "_")))),
    reactive(input$grp_plot_width),
    reactive(input$grp_plot_height)
  )

  output$download_gene_family_table <- downloadHandler(
    filename = function() paste0("gene_family_", gsub("[^A-Za-z0-9]+", "_", paste(input$gene_family_select %||% "family", collapse = "_")), "_", Sys.Date(), ".csv"),
    content = function(file) {
      d <- gene_family_data_reactive()
      req(!is.null(d))
      keep_cols <- c("gene_id", "hgnc_id", "Gene_Family", "Sub_Family", "Symbol", "log2FoldChange", "pValue", "padj", "DE_class", "short_description")
      d <- d[, intersect(keep_cols, names(d)), drop = FALSE]
      if ("gene_id" %in% names(d)) d <- d[!duplicated(d$gene_id), , drop = FALSE]
      write.csv(d, file, row.names = FALSE)
    }
  )

  observeEvent(input$run_gene_family_enrichment, {
    req(rv$de, gene_family_available())
    append_log("Running gene family enrichment for", rv$selected_organism_label, level = "STEP")
    withProgress(message = "Running gene family enrichment...", value = 0.5, {
      tryCatch({
        if (is.null(rv$gene_family_context)) {
          incProgress(0.2, detail = "Loading gene family reference data")
          rv$gene_family_context <- setup_gene_family_analysis(
            rv$de,
            tax_id = rv$selected_tax_id,
            alpha = input$alpha,
            lfc_cutoff = input$lfc_cutoff,
            gene_id_type = input$go_keytype %||% "TAIR",
            orgdb = input$go_orgdb %||% "org.At.tair.db"
          )
        }
        res <- run_gene_family_enrichment(
          rv$de,
          tax_id = rv$selected_tax_id,
          direction = input$gene_family_enrichment_direction %||% "up",
          alpha = input$alpha,
          lfc_cutoff = input$lfc_cutoff,
          min_set_size = input$gene_family_min_size %||% 3,
          gene_id_type = input$go_keytype %||% "TAIR",
          orgdb = input$go_orgdb %||% "org.At.tair.db",
          ctx = rv$gene_family_context
        )
        rv$gene_family_enrichment <- res
        append_log("Gene family enrichment completed:", nrow(res), "families tested.")
      }, error = function(e) {
        showNotification(paste("Gene family enrichment error:", e$message), type = "error", duration = 15)
        append_log("Gene family enrichment error:", e$message)
      })
    })
  })

  gene_family_enrichment_plot_reactive <- reactive({
    req(rv$gene_family_enrichment)
    plot_at_gene_family_enrichment(
      rv$gene_family_enrichment,
      p_value_threshold = input$alpha,
      top_n = 10,
      color_up = input$color_up %||% "#B2182B",
      color_down = input$color_down %||% "#2166AC",
      plot_theme = input$plot_theme %||% "classic",
      font_family = input$plot_font_family %||% "serif"
    )
  })

  output$gene_family_enrichment_plot <- renderPlot({
    if (is.null(rv$gene_family_enrichment)) {
      plot.new(); text(0.5, 0.5, "Run family enrichment first.", cex = 1.1)
    } else {
      tryCatch({
        p <- gene_family_enrichment_plot_reactive()
        if (is.null(p)) {
          plot.new(); text(0.5, 0.5, "No gene-family enrichment rows are available to plot.", cex = 1.1)
        } else {
          p
        }
      }, error = function(e) {
        plot.new(); text(0.5, 0.5, paste("Gene-family plot error:", e$message), cex = 1.1)
      })
    }
  }, width = function() input$grp_plot_width, height = function() input$grp_plot_height)

  output$gene_family_enrichment_table <- renderDT({
    req(rv$gene_family_enrichment)
    datatable(rv$gene_family_enrichment, rownames = FALSE, filter = "top",
              options = list(pageLength = 15, scrollX = TRUE))
  })

  output$download_gene_family_enrichment <- download_plot_server(
    gene_family_enrichment_plot_reactive,
    reactive(input$format_gene_family_enrichment),
    "gene_family_enrichment",
    reactive(input$grp_plot_width),
    reactive(input$grp_plot_height)
  )

  output$download_gene_family_enrichment_table <- downloadHandler(
    filename = function() paste0("gene_family_enrichment_", input$gene_family_enrichment_direction %||% "up", "_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$gene_family_enrichment)
      write.csv(rv$gene_family_enrichment, file, row.names = FALSE)
    }
  )

  # ---- KEGG Analysis -------------------------------------------------------
  observeEvent(input$run_kegg, {
    req(rv$de)
    append_log("Running KEGG enrichment for", input$kegg_species %||% "ath", level = "STEP")
    withProgress(message = "Running KEGG Enrichment...", value = 0.3, {
      tryCatch({
        res <- run_kegg_enrichment(
          rv$de,
          padj_cutoff = input$alpha,
          lfc_cutoff = input$lfc_cutoff,
          kegg_species = input$kegg_species %||% "ath",
          gene_id_type = input$go_keytype %||% NULL,
          orgdb = input$go_orgdb %||% NULL
        )
        rv$kegg_enrichment <- res
        
        incProgress(0.5, detail = "Preparing plot")
        rv$kegg_bubble <- NULL
        
        append_log("KEGG enrichment completed:", if (!is.null(res)) nrow(res) else 0, "pathways tested.")
      }, error = function(e) {
        showNotification(paste("KEGG Error:", e$message), type = "error", duration = 15)
        append_log("KEGG error:", e$message)
      })
    })
  })
  
  output$kegg_enrichment_table <- renderDT({
    req(rv$kegg_enrichment)
    datatable(rv$kegg_enrichment, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE))
  })
  
  kegg_bubble_reactive <- reactive({
    req(rv$kegg_enrichment)
    plot_kegg_bubble(
      rv$kegg_enrichment,
      p_value_threshold = input$alpha,
      color_up = input$color_up %||% "#B2182B",
      color_down = input$color_down %||% "#2166AC",
      plot_theme = input$plot_theme %||% "classic",
      font_family = input$plot_font_family %||% "serif"
    )
  })
  
  output$kegg_bubble_plot <- renderPlot({
    p <- tryCatch(kegg_bubble_reactive(), error = function(e) NULL)
    if (is.null(p)) {
      plot.new(); text(0.5, 0.5, "No significantly enriched pathways.", cex = 1.1)
    } else {
      p
    }
  }, width = function() input$kegg_plot_width, height = function() input$kegg_plot_height)
  
  output$download_kegg_table <- downloadHandler(
    filename = function() paste0("KEGG_enrichment_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$kegg_enrichment)
      write.csv(rv$kegg_enrichment, file, row.names = FALSE)
    }
  )
  
  output$download_kegg_bubble <- download_plot_server(
    kegg_bubble_reactive,
    reactive(input$format_kegg_bubble),
    "KEGG_enrichment_bubble",
    reactive(input$kegg_plot_width),
    reactive(input$kegg_plot_height)
  )

  # ---- Pathview Visualization ----------------------------------------------
  output$pathview_selector_ui <- renderUI({
    req(rv$kegg_enrichment)
    sig_paths <- rv$kegg_enrichment[rv$kegg_enrichment$p.value <= input$alpha, ]
    req(nrow(sig_paths) > 0)
    choices <- setNames(sig_paths$pathway.code, paste0(sig_paths$pathway.code, " - ", sig_paths$pathway.name))
    selectInput("pathview_select", "Select Enriched Pathway", choices = choices)
  })

  pathview_result_table <- function(mapped_de, pv_res = NULL) {
    if (is.null(mapped_de) || nrow(mapped_de) == 0) return(data.frame())
    de <- as.data.frame(mapped_de, check.names = FALSE)
    if (!"original_gene_id" %in% names(de)) de$original_gene_id <- de$gene_id

    map_ids <- character()
    plot_gene <- tryCatch(pv_res$plot.data.gene, error = function(e) NULL)
    if (!is.null(plot_gene) && nrow(plot_gene) > 0) {
      id_col <- first_existing_col(plot_gene, c("kegg.names", "kegg.name", "all.mapped", "labels", "label"))
      if (!is.null(id_col)) {
        map_ids <- unlist(strsplit(as.character(plot_gene[[id_col]]), "[,;[:space:]]+"))
        map_ids <- unique(map_ids[!is.na(map_ids) & nzchar(map_ids)])
      }
    }
    if (length(map_ids) > 0) {
      map_keys <- unique(sub("^.*:", "", map_ids))
      gene_keys <- sub("^.*:", "", as.character(de$gene_id))
      original_keys <- sub("^.*:", "", as.character(de$original_gene_id))
      de <- de[gene_keys %in% map_keys | original_keys %in% map_keys, , drop = FALSE]
    }
    if (nrow(de) == 0) return(data.frame())

    symbol_col <- first_existing_col(de, c("Symbol", "symbol", "gene_symbol", "Gene.symbol"))
    p_col <- first_existing_col(de, c("pValue", "pvalue", "p.value", "P.Value", "padj"))
    ec_col <- first_existing_col(de, c("EC", "ec", "EC_number", "EC.number", "EC_numbers", "Enzyme", "enzyme", "KEGG_EC"))

    out <- data.frame(
      gene_id = as.character(de$original_gene_id),
      Symbol = if (!is.null(symbol_col)) as.character(de[[symbol_col]]) else "",
      logFC = suppressWarnings(as.numeric(de$log2FoldChange)),
      pvalue = if (!is.null(p_col)) suppressWarnings(as.numeric(de[[p_col]])) else NA_real_,
      EC = if (!is.null(ec_col)) as.character(de[[ec_col]]) else "",
      stringsAsFactors = FALSE
    )
    out <- out[!duplicated(out$gene_id), , drop = FALSE]
    out[order(out$pvalue, na.last = TRUE), , drop = FALSE]
  }

  generate_pathview_map <- function(pwid, progress_message = "Generating Pathview map (downloading from KEGG)...") {
    withProgress(message = progress_message, value = 0.5, {
      tryCatch({
        rv$pathview_table <- NULL
        # Filter for significant genes only so non-sig genes remain uncolored
        sig_df <- rv$de[
          !is.na(rv$de$padj) & !is.na(rv$de$log2FoldChange) &
            rv$de$padj < input$alpha & abs(rv$de$log2FoldChange) >= input$lfc_cutoff,
        ]
        sig_df <- map_de_ids_for_kegg(
          sig_df,
          gene_id_type = input$go_keytype %||% NULL,
          orgdb = input$go_orgdb %||% NULL
        )
        lfc <- sig_df$log2FoldChange
        names(lfc) <- sig_df$gene_id
        lfc <- lfc[!is.na(lfc)]
        lfc_agg <- stats::aggregate(lfc, by = list(gene_id = names(lfc)), FUN = function(x) x[which.max(abs(x))])
        lfc <- stats::setNames(lfc_agg$x, lfc_agg$gene_id)
        if (length(lfc) == 0) {
          stop("No significant genes could be mapped to KEGG IDs for Pathview. Check Gene ID type, OrgDb, and KEGG code.")
        }
        append_log("Pathview mapped", length(lfc), "significant genes to KEGG IDs.")

        run_id <- gsub("[^0-9A-Za-z]", "", paste0(format(Sys.time(), "%Y%m%d%H%M%OS3"), sample.int(1000000, 1)))
        pv_dir <- file.path(tempdir(), paste0("pathview_", pwid, "_", run_id))
        dir.create(pv_dir, showWarnings = FALSE)
        
        old_wd <- setwd(pv_dir)
        on.exit(setwd(old_wd), add = TRUE)
        
        # Load pathview explicitly to make internal datasets (like 'bods') available
        library(pathview)
        data("bods", package = "pathview")
        gene_limit <- suppressWarnings(as.numeric(input$pathview_gene_limit %||% 1))
        if (!is.finite(gene_limit) || is.na(gene_limit) || gene_limit <= 0) gene_limit <- 1
        
        pv_res <- pathview(
          gene.data = lfc,
          pathway.id = pwid,
          species = input$kegg_species %||% "ath",
          limit = list(gene = gene_limit, cpds = 1),
          bins = list(gene = 10, cpds = 10),
          low = list(gene = input$color_down %||% "#2166AC", cpds = input$color_up %||% "#2166AC"),
          mid = list(gene = input$pathview_color_mid %||% "#d4d4d4", cpds = input$pathview_color_mid %||% "#d4d4d4"),
          high = list(gene = input$color_up %||% "#B2182B", cpds = input$pathview_color_up %||% "#B2182B"),
          plot.col.key = isTRUE(input$pathview_plot_col_key),
          key.pos = input$pathview_key_pos %||% "topright",
          new.signature = FALSE,
          # Genes are mapped to KEGG IDs before calling Pathview; avoid forcing
          # Arabidopsis-only TAIR IDs for featureCounts or non-ATH inputs.
          gene.idtype = "KEGG"
        )
        rv$pathview_table <- pathview_result_table(sig_df, pv_res)
        
        img_path <- file.path(pv_dir, paste0(pwid, ".pathview.png"))
        if (file.exists(img_path)) {
          rv$pathview_pathway <- img_path
          rv$pathview_pwid <- pwid
          append_log("Pathview generated:", pwid)
        } else {
          showNotification("Pathview did not generate an image.", type = "warning")
        }
      }, error = function(e) {
        showNotification(paste("Pathview Error:", e$message), type = "error", duration = 15)
        append_log("Pathview error:", e$message)
      })
    })
  }

  observeEvent(input$pathview_select, {
    rv$pathview_pathway <- NULL
    rv$pathview_pwid <- NULL
    rv$pathview_table <- NULL
  }, ignoreInit = TRUE)
  
  observeEvent(input$run_pathview, {
    req(rv$de, input$pathview_select)
    append_log("Generating Pathview map:", input$pathview_select, level = "STEP")
    generate_pathview_map(input$pathview_select)
  })

  output$pathview_image_ui <- renderUI({
    if (is.null(rv$pathview_pathway)) {
      div(class = "text-muted", "Select a pathway and click 'Generate Pathview Map'")
    } else {
      tagList(
        imageOutput("pathview_render", width = "100%", height = "auto"),
        div(class = "download-row", downloadButton("download_pathview", "Download Map PNG"))
      )
    }
  })
  
  output$pathview_render <- renderImage({
    req(rv$pathview_pathway)
    list(src = rv$pathview_pathway, contentType = "image/png", alt = "Pathview Map", width = "100%")
  }, deleteFile = FALSE)

  output$pathview_table_ui <- renderUI({
    req(rv$pathview_pathway)
    tagList(
      tags$hr(),
      h4("Pathview mapped genes"),
      DTOutput("pathview_table")
    )
  })

  output$pathview_table <- renderDT({
    req(rv$pathview_table)
    button_data <- add_gene_count_buttons(rv$pathview_table)
    datatable(button_data$data, rownames = FALSE, filter = "top", escape = !button_data$has_buttons,
              options = list(
                pageLength = 15,
                scrollX = TRUE,
                autoWidth = FALSE,
                columnDefs = if (button_data$has_buttons) gene_count_button_defs else list()
              ),
              callback = if (button_data$has_buttons) gene_count_button_callback else JS(""))
  })
  
  output$download_pathview <- downloadHandler(
    filename = function() paste0("pathview_", basename(rv$pathview_pathway)),
    content = function(file) {
      req(rv$pathview_pathway)
      file.copy(rv$pathview_pathway, file)
    }
  )
}



shinyApp(ui, server)
