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

# Allow file uploads up to 30MB
options(shiny.maxRequestSize = 30 * 1024^2)

source(file.path("R", "helpers.R"), local = TRUE)
source(file.path("legacy_scripts", "volcano_TEG_overlap_with_TE_families_RNAseq.R"), local = TRUE)
source(file.path("legacy_scripts", "genes_into_groups.R"), local = TRUE)
source(file.path("legacy_scripts", "kegg_analysis.R"), local = TRUE)

wrap_html <- function(x, width = 45) {
  HTML(paste(strwrap(x, width = width), collapse = "<br>"))
}

example_csv <- file.path("example_data", "all_genes_results_mto1_vs_wt.csv")
te_file_path <- file.path("description_files", "TAIR10_Transposable_Elements.txt")
te_super_family_choices <- tryCatch({
  te_tbl <- load_te_file(te_file_path)
  sort(unique(te_tbl$Transposon_Super_Family))
}, error = function(e) character())

ui <- fluidPage(
  theme = shinythemes::shinytheme("united"),
  tags$head(tags$style(HTML("\n    .app-title { margin-top: 10px; margin-bottom: 4px; font-weight: 700; }\n    .muted { color: #666; font-size: 0.92em; }\n    .tab-content { padding: 16px; border: 1px solid #ddd; border-top: none; }\n    .download-row .btn { margin-right: 8px; margin-top: 6px; }\n    pre { white-space: pre-wrap; }\n  "))),

  # Add the theme selector if shinythemes is installed (for easy testing of themes)
  # if (requireNamespace("shinythemes", quietly = TRUE)) shinythemes::themeSelector(),

  titlePanel(div(class = "app-title", "RNA-seq Analysis Dashboard")),
  div(class = "muted", "Load DE results directly, or run DESeq2 from RSEM .genes.results files. Outputs are shown in the app and downloaded only when requested."),
  br(),

  sidebarLayout(
    sidebarPanel(width = 3,
      h4("1. Data input"),
      radioButtons("data_mode", NULL,
        choices = c(
          "Upload DE results CSV/TSV" = "csv",
          "Run DESeq2 from RSEM folder" = "rsem",
          "Example dataset (mto1 vs. wt)" = "example"
        ), selected = "csv"),

      conditionalPanel("input.data_mode == 'csv'",
        fileInput("de_file", "DE results table", accept = c(".csv", ".tsv", ".txt")),
        div(class = "muted", "Required columns: gene_id, log2FoldChange, padj. baseMean enables MA plot.")
      ),

      conditionalPanel("input.data_mode == 'rsem'",
        shinyDirButton("choose_rsem_dir", "Choose RSEM folder", "Select a folder"),
        textInput("rsem_path", "RSEM folder path", value = ""),
        actionButton("scan_rsem", "Scan folder", class = "btn-primary"),
        tags$hr(),
        fileInput("coldata_file", "Optional colData CSV/TSV", accept = c(".csv", ".tsv", ".txt")),
        div(class = "muted", "colData can use columns sample_id/condition, or the older x/sample/exp format. You can also edit the scanned table below."),
        uiOutput("contrast_ui"),
        checkboxInput("lfc_shrink", "Use ashr LFC shrinkage", value = FALSE),
        numericInput("min_count", "Filter genes: min total count", value = 10, min = 0, step = 1),
        actionButton("run_deseq", "Run DESeq2", class = "btn-success")
      ),

      tags$hr(),
      h4("2. Thresholds"),
      numericInput("alpha", "padj cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
      numericInput("lfc_cutoff", "|log2FC| cutoff", value = 1, min = 0, step = 0.25),
      numericInput("go_pcut", "GO p-value display cutoff", value = 0.01, min = 0, max = 1, step = 0.001),
      selectInput("ontology", "GO ontology", choices = c("BP", "MF", "CC"), selected = "BP"),
      numericInput("top_n", "Top GO terms to display", value = 20, min = 5, max = 200, step = 5),
      tags$hr(),
      h4("3. Plot Settings"),

      # Per-tab width/height sliders with sensible defaults per plot type
      conditionalPanel("input.tabs == 'DE plots'",
        sliderInput("de_plot_width",  "Plot width (px)",  min = 200, max = 1600, value = 350, step = 50),
        sliderInput("de_plot_height", "Plot height (px)", min = 200, max = 1200, value = 200, step = 50),
        tags$hr(),
        h5("PCA options"),
        checkboxInput("pca_show_labels", "Show sample labels", value = TRUE),
        uiOutput("pca_conditions_ui")
      ),
      conditionalPanel("input.tabs == 'GO analysis'",
        sliderInput("go_plot_width",  "Plot width (px)",  min = 200, max = 1600, value = 450, step = 50),
        sliderInput("go_plot_height", "Plot height (px)", min = 200, max = 1200, value = 400, step = 50)
      ),
      conditionalPanel("input.tabs == 'TE analysis'",
        sliderInput("teg_plot_width",  "Plot width (px)",  min = 200, max = 1600, value = 550, step = 50),
        sliderInput("teg_plot_height", "Plot height (px)", min = 100, max = 1200, value = 200, step = 50)
      ),
      conditionalPanel("input.tabs == 'Gene groups'",
        sliderInput("grp_plot_width",  "Plot width (px)",  min = 200, max = 1600, value = 350, step = 50),
        sliderInput("grp_plot_height", "Plot height (px)", min = 200, max = 1200, value = 200, step = 50)
      ),
      conditionalPanel("input.tabs == 'KEGG analysis'",
        sliderInput("kegg_plot_width",  "Plot width (px)",  min = 200, max = 1600, value = 750, step = 50),
        sliderInput("kegg_plot_height", "Plot height (px)", min = 200, max = 1200, value = 300, step = 50)
      ),
      # Point aesthetics (global)
      sliderInput("plot_point_size", "Point size", min = 0.1, max = 5, value = 1, step = 0.1),
      sliderInput("plot_alpha", "Point alpha", min = 0.1, max = 1, value = 0.65, step = 0.05),
      tags$hr(),
      tags$label("Point colors", style = "font-weight: 600; font-size: 0.95em;"),
      tags$div(style = "display: flex; align-items: center; gap: 8px; margin-top: 4px;",
        tags$input(type = "color", id = "color_up", value = "#B2182B",
          style = "width: 42px; height: 30px; border: none; cursor: pointer; padding: 0;",
          oninput = "Shiny.setInputValue('color_up', this.value, {priority: 'event'})"),
        tags$span("Upregulated", style = "font-size: 0.9em;")
      ),
      tags$div(style = "display: flex; align-items: center; gap: 8px; margin-top: 4px;",
        tags$input(type = "color", id = "color_down", value = "#2166AC",
          style = "width: 42px; height: 30px; border: none; cursor: pointer; padding: 0;",
          oninput = "Shiny.setInputValue('color_down', this.value, {priority: 'event'})"),
        tags$span("Downregulated", style = "font-size: 0.9em;")
      ),
      tags$div(style = "display: flex; align-items: center; gap: 8px; margin-top: 4px;",
        tags$input(type = "color", id = "color_ns", value = "#B3B3B3",
          style = "width: 42px; height: 30px; border: none; cursor: pointer; padding: 0;",
          oninput = "Shiny.setInputValue('color_ns', this.value, {priority: 'event'})"),
        tags$span("Not significant", style = "font-size: 0.9em;")
      )
    ),


    mainPanel(width = 9,
      tabsetPanel(id = "tabs",
        tabPanel("Data",
          fluidRow(
            column(4, verbatimTextOutput("data_summary")),
            column(8,
              conditionalPanel("input.data_mode == 'rsem'",
                h4("Editable colData"),
                DTOutput("coldata_table"),
                div(class = "muted", "Edit the condition column, then select treatment/control in the sidebar and run DESeq2.")
              ),
              h4("DE table preview"),
              DTOutput("de_preview"),
              div(class = "download-row", downloadButton("download_de", "Download DE table"), downloadButton("download_norm_counts", "Download normalized counts"))
            )
          )
        ),

        tabPanel("DE plots",
          fluidRow(
            column(12,
              h4("Volcano plot"),
              plotOutput("volcano_plot", width = "auto", height = "auto"),
              download_plot_ui("volcano", "Download volcano plot")
            )
          ),
          tags$hr(),
          fluidRow(
            column(12,
              h4("MA plot"),
              plotOutput("ma_plot", width = "auto", height = "auto"),
              textOutput("ma_message"),
              download_plot_ui("ma", "Download MA plot")
            )
          ),
          tags$hr(),
          fluidRow(
            column(12,
              h4("PCA"),
              plotOutput("pca_plot", width = "auto", height = "auto"),
              textOutput("pca_message"),
              download_plot_ui("pca", "Download PCA plot"),
              div(class = "download-row", downloadButton("download_pca_table", "Download PCA table"))
            )
          )
        ),

        tabPanel("GO analysis",
          tabsetPanel(
            tabPanel("Enrichment",
              wellPanel(
                fluidRow(
                  column(3, selectInput("go_direction", "Gene set", choices = c("up", "down", "all"), selected = "up")),
                  column(3, selectInput("go_algorithm", "topGO algorithm", choices = c("weight01", "classic", "elim"), selected = "weight01")),
                  column(3, selectInput("go_statistic", "topGO statistic", choices = c("fisher", "ks"), selected = "fisher")),
                  column(3, br(), actionButton("run_go", "Run GO enrichment", class = "btn-primary", style = "width:100%;"))
                ),
                div(class = "muted", "Runs all terms; p-value & top-N filters apply on display without re-running.")
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
            tabPanel("REVIGO-like reduction",
              wellPanel(
                div(class = "muted", "Uses cached GO results. Run GO enrichment for both up and down genes first."),
                fluidRow(
                  column(4, numericInput("revigo_top_n", "Max GO terms", value = 80, min = 5, max = 300, step = 5)),
                  column(4, sliderInput("revigo_threshold", "Similarity reduction threshold", min = 0.4, max = 0.95, value = 0.7, step = 0.05)),
                  column(4, br(), actionButton("run_revigo", "Run REVIGO-like analysis", class = "btn-primary", style = "width:100%;"))
                )
              ),
              fluidRow(
                column(12,
                  h4("Upregulated GO terms"),
                  plotOutput("revigo_up_plot", width = "auto", height = "auto"),
                  download_plot_ui("revigo_up", "Download up plot"),
                  div(class = "download-row", downloadButton("download_revigo_up_table", "Download up table")),
                  tags$hr(),
                  h4("Downregulated GO terms"),
                  plotOutput("revigo_down_plot", width = "auto", height = "auto"),
                  download_plot_ui("revigo_down", "Download down plot"),
                  div(class = "download-row", downloadButton("download_revigo_down_table", "Download down table"))
                )
              )
            ),
            tabPanel("GO offspring",
              wellPanel(
                fluidRow(
                  column(8, textAreaInput("parent_go_ids", "Parent GO IDs for offspring summary", rows = 2,
                                value = "GO:0006950, GO:0009628, GO:0009408, GO:0006979, GO:0009414, GO:0009611")),
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
            tabPanel("Abiotic stress",
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

        tabPanel("KEGG analysis",
          tabsetPanel(
            tabPanel("Enrichment Analysis",
              wellPanel(
                fluidRow(
                  column(4, actionButton("run_kegg", "Run KEGG Enrichment", class = "btn-primary", style="width:100%")),
                  column(8, div(class = "muted", "Uses Wilcoxon rank-sum test on dataset p-values against Arabidopsis KEGG pathways. May take a moment the first time to fetch pathway genes from KEGG servers."))
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
                div(class = "muted", style = "margin-top: 8px;", "Downloads pathway mapping from KEGG and colors genes by their log2FoldChange. Requires internet access.")
              ),
              fluidRow(
                column(12,
                  h4("Pathview Map"),
                  uiOutput("pathview_image_ui")
                )
              )
            )
          )
        ),

        tabPanel("TE analysis",
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
                div(class = "muted", "Uses description_files/Methylome.At_description_file.csv.gz and description_files/TAIR10_Transposable_Elements.txt")
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

        tabPanel("Gene groups",
          wellPanel(
            fluidRow(
              column(8, uiOutput("gene_group_selector_ui")),
              column(4, br(), actionButton("run_gene_groups", "Build gene groups", class = "btn-primary", style = "width:100%;"))
            ),
            div(class = "muted", style = "margin-top: 8px;", "Downloads gene-set lists from GitHub. Requires internet access and GO.db package.")
          ),
          fluidRow(
            column(12,
              h4("Volcano plot – highlighted group"),
              plotOutput("gene_group_volcano", width = "auto", height = "auto"),
              download_plot_ui("gene_group_volcano", "Download group volcano"),
              tags$hr(),
              h4("Genes in selected group"),
              DTOutput("gene_group_table"),
              div(class = "download-row", downloadButton("download_gene_group_table", "Download group table"))
            )
          )
        ),

        tabPanel("Search Annotations",
          fluidRow(
            column(12,
              br(),
              h4("Search Description File"),
              div(class = "muted", "Search across all annotation columns! Use a space for 'AND' (e.g. 'kinase stress'), or use | for 'OR' (e.g. 'kinase|stress'). You can also filter specific columns using the boxes below the headers."),
              tags$hr(),
              DTOutput("search_annotations_table"),
              div(class = "download-row", downloadButton("download_search_annotations", "Download Table"))
            )
          )
        ),

        tabPanel("Log / Help",
          tabsetPanel(
            tabPanel("Log",
              br(),
              h4("Run Log"),
              div(class = "muted", "Execution events, errors, and intermediate steps are recorded here for debugging."),
              verbatimTextOutput("run_log")
            ),
            tabPanel("Help",
              br(),
              h4("Notes & Usage"),
              tags$ul(
                tags$li("Volcano works from any DE table with gene_id, log2FoldChange and padj. If columns are not found by name, it will use the first 3 columns automatically."),
                tags$li("MA requires baseMean; this is produced by DESeq2/RSEM or can come in an uploaded CSV."),
                tags$li("PCA requires expression/count data and is therefore produced only after running DESeq2 from RSEM files."),
                tags$li("GO offspring/stress summaries use GO biological-process annotations from the CSV if present; otherwise the app tries org.At.tair.db."),
                tags$li("No plots or tables are saved automatically. Use the download buttons.")
              ),
              tags$hr(),
              h4("AI Developer Instructions"),
              div(class = "muted", "Guidelines for future feature loading & debugging by AI assistants:"),
              tags$ul(
                tags$li(strong("Architecture: "), "The app uses a standard Shiny UI/Server layout in ", code("app.R"), ". Heavy logic, data processing, and visualization functions reside in ", code("R/helpers.R"), " (or legacy scripts). Always add new large processing functions to ", code("helpers.R"), " to keep the app file clean."),
                tags$li(strong("Adding Features: "), "To add a new tab, define the UI under ", code("mainPanel > tabsetPanel"), " and the corresponding reactive logic in the ", code("server"), " function. Use the established pattern of caching results in ", code("rv"), " (reactiveValues) to prevent unnecessary recalculations."),
                tags$li(strong("Debugging: "), "Check the 'Log' tab first. Use ", code("append_log()"), " in the server code to expose intermediate steps, variable states, and errors directly to the user interface."),
                tags$li(strong("Data Standard: "), "Functions in ", code("helpers.R"), " should rely on the standardized DE table format. When introducing new inputs, ensure they pass through ", code("standardize_de_table()"), " or similar validation functions.")
              ),
              tags$hr(),
              h4("Source Code & Repository"),
              div(class = "muted", "You can download, edit, and access the raw files and source code for this app from its GitHub repository:"),
              tags$a(href = "https://github.com/Yo-yerush/RNAseq_analysis_app.git", "https://github.com/Yo-yerush/RNAseq_analysis_app.git", target = "_blank"),
              tags$hr(),
              h4("Author"),
              div("Yonatan Yerushalmy"),
              div(class = "muted", "Plant's Metabolism and Molecular Genetic laboratory"),
              div(class = "muted", "Prof' Rachel Amir's group")
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
    de = NULL,
    norm_counts = NULL,
    pca = NULL,
    coldata = NULL,
    de_summary = NULL,
    go_trigger = 0,
    offspring = NULL,
    stress = NULL,
    stress_plot = NULL,
    revigo_up = NULL,
    revigo_down = NULL,
    te_volcano = NULL,
    te_volcano_plot = NULL,
    gene_groups = NULL,
    gene_group_plots = NULL,
    gg_context = NULL,
    gg_cache = list(),
    kegg_enrichment = NULL,
    kegg_bubble = NULL,
    pathview_pathway = NULL,
    go_cache = list(),
    log = character()
  )

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  append_log <- function(...) {
    timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
    msg <- paste(..., collapse = " ")
    rv$log <- c(rv$log, paste(timestamp, "-", msg))
    
    # Show a UI notification for successful actions so the user sees what happened
    # We skip error messages because they already have dedicated showNotification calls
    if (!grepl("error|no rsem gene.results files found", tolower(msg))) {
      showNotification(msg, type = "message", duration = 4)
    }
  }

  go_cache_key <- function(direction) {
    paste(direction, input$ontology, input$alpha, input$lfc_cutoff,
          input$go_algorithm, input$go_statistic, sep = "|")
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
                              statistic = input$go_statistic)
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
                        input$top_n, direction = input$go_direction, point_alpha = input$plot_alpha)
  })

  observeEvent(input$choose_rsem_dir, {
    path <- shinyFiles::parseDirPath(volumes, input$choose_rsem_dir)
    if (length(path) && nzchar(path)) updateTextInput(session, "rsem_path", value = path)
  })

  observeEvent(input$scan_rsem, {
    req(input$rsem_path)
    tbl <- scan_rsem_files(input$rsem_path)
    if (nrow(tbl) == 0) {
      showNotification("No .genes.results files found in this folder", type = "error")
      append_log("No RSEM gene.results files found in", input$rsem_path)
      return()
    }
    rv$coldata <- tbl[, c("sample_id", "condition", "sample_label")]
    append_log("Scanned", nrow(tbl), "RSEM gene.results files.")
  })

  observeEvent(input$coldata_file, {
    req(input$coldata_file$datapath)
    cd <- read_any_table(input$coldata_file$datapath)
    cd <- normalize_coldata(cd)
    rv$coldata <- cd
    append_log("Loaded colData with", nrow(cd), "samples.")
  })

  output$coldata_table <- renderDT({
    req(rv$coldata)
    datatable(rv$coldata, editable = TRUE, rownames = FALSE, options = list(pageLength = 12, scrollX = TRUE))
  })

  observeEvent(input$coldata_table_cell_edit, {
    info <- input$coldata_table_cell_edit
    rv$coldata[info$row, info$col + 1] <- DT::coerceValue(info$value, rv$coldata[info$row, info$col + 1])
  })

  output$contrast_ui <- renderUI({
    cd <- rv$coldata
    if (is.null(cd) || !"condition" %in% names(cd)) return(div(class = "muted", "Scan folder or load colData first."))
    conditions <- unique(cd$condition)
    tagList(
      selectInput("treatment", "Treatment", choices = conditions, selected = conditions[min(2, length(conditions))]),
      selectInput("control", "Control", choices = conditions, selected = conditions[1])
    )
  })

  loaded_de <- reactive({
    if (input$data_mode == "example") {
      validate(need(file.exists(example_csv), "Example CSV was not found."))
      standardize_de_table(read_any_table(example_csv))
    } else if (input$data_mode == "csv") {
      req(input$de_file$datapath)
      standardize_de_table(read_any_table(input$de_file$datapath))
    } else {
      rv$de
    }
  })

  observe({
    if (input$data_mode %in% c("example", "csv")) {
      rv$de <- tryCatch(loaded_de(), error = function(e) {
        showNotification(e$message, type = "error", duration = 8)
        NULL
      })
      rv$norm_counts <- NULL
      rv$pca <- NULL
      rv$de_summary <- NULL
      rv$go_trigger <- 0
      rv$offspring <- NULL
      rv$stress <- NULL
      rv$stress_plot <- NULL
      rv$revigo_up <- NULL
      rv$revigo_down <- NULL
      rv$te_volcano <- NULL
      rv$te_volcano_plot <- NULL
      rv$gene_groups <- NULL
      rv$gene_group_plots <- NULL
      rv$gg_context <- NULL
      rv$gg_cache <- list()
      rv$kegg_enrichment <- NULL
      rv$kegg_bubble <- NULL
      rv$pathview_pathway <- NULL
      rv$go_cache <- list()
    }
  })

  observeEvent(input$run_deseq, {
    req(input$rsem_path, rv$coldata, input$treatment, input$control)
    withProgress(message = "Running DESeq2", value = 0.1, {
      tryCatch({
        incProgress(0.2, detail = "Importing RSEM files")
        res <- run_deseq2_from_rsem(
          folder = input$rsem_path,
          coldata = rv$coldata,
          treatment = input$treatment,
          control = input$control,
          lfc_shrink = isTRUE(input$lfc_shrink),
          min_count = input$min_count
        )
        incProgress(0.7, detail = "Preparing tables and PCA")
        rv$de <- res$de_table
        rv$norm_counts <- res$norm_counts
        rv$pca <- res$pca_table
        rv$de_summary <- res$summary
        rv$go_trigger <- 0
        rv$offspring <- NULL
        rv$stress <- NULL
        rv$stress_plot <- NULL
        rv$revigo_up <- NULL
        rv$revigo_down <- NULL
        rv$te_volcano <- NULL
        rv$te_volcano_plot <- NULL
        rv$go_cache <- list()
        append_log("DESeq2 finished:", input$treatment, "vs", input$control, "with", nrow(rv$de), "tested genes.")
        updateTabsetPanel(session, "tabs", selected = "DE plots")
      }, error = function(e) {
        showNotification(e$message, type = "error", duration = 12)
        append_log("DESeq2 error:", e$message)
      })
    })
  })

  output$data_summary <- renderText({
    df <- rv$de
    if (is.null(df)) return("No DE table loaded yet.")
    dfc <- classify_de(df, input$alpha, input$lfc_cutoff)
    paste0(
      "Genes: ", nrow(dfc), "\n",
      "Up: ", sum(dfc$DE_class == "up", na.rm = TRUE), "\n",
      "Down: ", sum(dfc$DE_class == "down", na.rm = TRUE), "\n",
      "Not significant: ", sum(dfc$DE_class == "not_significant", na.rm = TRUE), "\n",
      "Has baseMean: ", "baseMean" %in% names(dfc), "\n",
      "Has PCA: ", !is.null(rv$pca), "\n\n",
      if (!is.null(rv$de_summary)) paste0("DESeq2 summary:\n", rv$de_summary) else ""
    )
  })

  output$de_preview <- renderDT({
    req(rv$de)
    keep_cols <- c("gene_id", "Symbol", "log2FoldChange", "padj")
    d <- rv$de[, intersect(keep_cols, names(rv$de)), drop = FALSE]
    datatable(head(d, 5000), rownames = FALSE, filter = "top", options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$search_annotations_table <- renderDT({
    req(rv$de)
    cols <- names(rv$de)
    first_cols <- intersect(c("gene_id", "Symbol", "log2FoldChange", "padj"), cols)
    other_cols <- setdiff(cols, first_cols)
    d <- rv$de[, c(first_cols, other_cols), drop = FALSE]
    
    # Use JavaScript to truncate long text cells and show full text on hover
    datatable(d, rownames = FALSE, filter = "top", 
      options = list(
        search = list(regex = TRUE, smart = TRUE),
        pageLength = 15, 
        scrollX = TRUE,
        columnDefs = list(
          list(
            targets = "_all",
            render = JS(
              "function(data, type, row, meta) {",
              "  if (type === 'display' && data != null && typeof data === 'string' && data.length > 50) {",
              "    return '<span title=\"' + data.replace(/\"/g, '&quot;') + '\">' + data.substr(0, 50) + '...</span>';",
              "  } else {",
              "    return data;",
              "  }",
              "}"
            )
          )
        )
      )
    )
  })

  volcano_reactive <- reactive({
    req(rv$de)
    make_volcano_plot(rv$de, input$alpha, input$lfc_cutoff, "Volcano plot",
                      point_size = input$plot_point_size, point_alpha = input$plot_alpha,
                      color_up = input$color_up %||% "#B2182B", color_down = input$color_down %||% "#2166AC", color_ns = input$color_ns %||% "#B3B3B3")
  })
  output$volcano_plot <- renderPlot({ volcano_reactive() }, width = function() input$de_plot_width, height = function() input$de_plot_height)

  ma_reactive <- reactive({
    req(rv$de)
    make_ma_plot(rv$de, input$alpha, input$lfc_cutoff, "MA plot",
                 point_size = input$plot_point_size, point_alpha = input$plot_alpha,
                 color_up = input$color_up %||% "#B2182B", color_down = input$color_down %||% "#2166AC", color_ns = input$color_ns %||% "#B3B3B3")
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
                  conditions = input$pca_conditions)
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

  observeEvent(input$run_offspring, {
    req(rv$de)
    tryCatch({
      rv$offspring <- make_go_offspring_summary(rv$de, input$parent_go_ids, input$alpha, input$lfc_cutoff)
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
    tryCatch({
      rv$stress <- make_abiotic_stress_table(rv$de, dataset = input$stress_dataset, alpha = input$alpha, lfc_cutoff = input$lfc_cutoff)
      rv$stress_plot <- make_abiotic_stress_plot(rv$stress, paste("Abiotic stress enrichment:", input$stress_dataset))
      append_log("Abiotic stress enrichment finished for", input$stress_dataset, "gene set.")
    }, error = function(e) {
      showNotification(e$message, type = "error", duration = 12)
      append_log("Stress enrichment error:", e$message)
    })
  })
  output$stress_plot <- renderPlot({
    if (is.null(rv$stress_plot)) { plot.new(); text(0.5, 0.5, "Run abiotic stress enrichment first.", cex = 1.1) } else rv$stress_plot
  }, width = function() input$go_plot_width, height = function() input$go_plot_height)
  output$stress_table <- renderDT({
    req(rv$stress)
    datatable(rv$stress, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  observeEvent(input$run_revigo, {
    req(rv$de)
    withProgress(message = "Running REVIGO-like semantic reduction", value = 0.1, {
      tryCatch({
        incProgress(0.2, detail = "Loading cached GO results")
        go_up <- get_cached_go("up")
        go_down <- get_cached_go("down")
        if (is.null(go_up) || is.null(go_down)) {
          stop("Run GO enrichment for both up and down genes first (same parameters) so REVIGO can reuse cached results.")
        }
        incProgress(0.2, detail = "Reducing up GO terms")
        rv$revigo_up <- if (!is.null(go_up) && nrow(go_up) >= 2) {
          run_rrvgo_reduce(go_up, ontology = input$ontology, top_n = input$revigo_top_n,
                           threshold = input$revigo_threshold, title = "Upregulated genes")
        } else {
          append_log("REVIGO-like note: not enough upregulated GO terms to reduce.")
          NULL
        }
        incProgress(0.2, detail = "Reducing down GO terms")
        rv$revigo_down <- if (!is.null(go_down) && nrow(go_down) >= 2) {
          run_rrvgo_reduce(go_down, ontology = input$ontology, top_n = input$revigo_top_n,
                           threshold = input$revigo_threshold, title = "Downregulated genes")
        } else {
          append_log("REVIGO-like note: not enough downregulated GO terms to reduce.")
          NULL
        }
        append_log("REVIGO-like analysis finished.")
      }, error = function(e) {
        showNotification(e$message, type = "error", duration = 12)
        append_log("REVIGO-like error:", e$message)
      })
    })
  })
  output$revigo_up_plot <- renderPlot({
    if (is.null(rv$revigo_up)) { plot.new(); text(0.5, 0.5, "Run REVIGO-like analysis first.", cex = 1.1) } else rv$revigo_up$plot
  }, width = function() input$go_plot_width, height = function() input$go_plot_height)
  output$revigo_down_plot <- renderPlot({
    if (is.null(rv$revigo_down)) { plot.new(); text(0.5, 0.5, "Run REVIGO-like analysis first.", cex = 1.1) } else rv$revigo_down$plot
  }, width = function() input$go_plot_width, height = function() input$go_plot_height)
  # ---- TE Analysis -------------------------------------------------------
  observeEvent(input$run_te_enrich, {
    req(rv$de)
    withProgress(message = "Running TE Enrichment...", value = 0.3, {
      tryCatch({
        res <- run_te_enrichment(rv$de, pvalue_cutoff = input$te_enrich_pvalue)
        rv$te_enrichment <- res
        
        incProgress(0.5, detail = "Generating Plot")
        rv$te_enrich_bubble <- plot_te_enrichment(res, p_value_threshold = input$te_enrich_pvalue)
        
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
    req(rv$te_enrich_bubble)
    rv$te_enrich_bubble
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
    withProgress(message = "Building TE volcano", value = 0.2, {
      tryCatch({
        out <- make_retro_te_volcano(
          rv$de,
          super_families = input$te_super_families,
          padj_cutoff = input$te_padj_cutoff,
          lfc_cutoff = input$te_lfc_cutoff,
          point_size = input$plot_point_size,
          point_alpha = input$plot_alpha
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
  output$te_volcano_plot <- renderPlot({
    if (is.null(rv$te_volcano_plot)) {
      plot.new(); text(0.5, 0.5, "Run TE volcano first.", cex = 1.1)
    } else {
      rv$te_volcano_plot
    }
  }, width = function() input$teg_plot_width, height = function() input$teg_plot_height)
  output$te_volcano_table <- renderDT({
    req(rv$te_volcano)
    datatable(rv$te_volcano, rownames = FALSE, filter = "top", options = list(pageLength = 12, scrollX = TRUE))
  })

  output$run_log <- renderText({ paste(rv$log, collapse = "\n") })

  output$download_de <- downloadHandler(
    filename = function() paste0("DE_results_", Sys.Date(), ".csv"),
    content = function(file) write.csv(rv$de, file, row.names = FALSE)
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
  output$download_offspring <- downloadHandler(
    filename = function() paste0("GO_offspring_summary_", Sys.Date(), ".csv"),
    content = function(file) { req(rv$offspring); write.csv(rv$offspring, file, row.names = FALSE) }
  )
  output$download_stress <- download_plot_server(reactive(rv$stress_plot), reactive(input$format_stress), "abiotic_stress", reactive(input$go_plot_width), reactive(input$go_plot_height))
  output$download_stress_table <- downloadHandler(
    filename = function() paste0("abiotic_stress_enrichment_", input$stress_dataset, "_", Sys.Date(), ".csv"),
    content = function(file) { req(rv$stress); write.csv(rv$stress, file, row.names = FALSE) }
  )
  output$download_revigo_up <- download_plot_server(reactive(rv$revigo_up$plot), reactive(input$format_revigo_up), "REVIGO_up", reactive(input$go_plot_width), reactive(input$go_plot_height))
  output$download_revigo_down <- download_plot_server(reactive(rv$revigo_down$plot), reactive(input$format_revigo_down), "REVIGO_down", reactive(input$go_plot_width), reactive(input$go_plot_height))
  output$download_revigo_up_table <- downloadHandler(
    filename = function() paste0("REVIGO_up_table_", input$ontology, "_", Sys.Date(), ".csv"),
    content = function(file) { req(rv$revigo_up); write.csv(rv$revigo_up$table, file, row.names = FALSE) }
  )
  output$download_revigo_down_table <- downloadHandler(
    filename = function() paste0("REVIGO_down_table_", input$ontology, "_", Sys.Date(), ".csv"),
    content = function(file) { req(rv$revigo_down); write.csv(rv$revigo_down$table, file, row.names = FALSE) }
  )
  output$download_te_enrich_bubble <- download_plot_server(te_enrich_bubble_reactive, reactive(input$format_te_enrich_bubble), "TE_Enrichment_bubble", reactive(input$teg_plot_width), reactive(input$teg_plot_height))
  output$download_te_enrich_table <- downloadHandler(
    filename = function() paste0("TE_Enrichment_table_", Sys.Date(), ".csv"),
    content = function(file) { req(rv$te_enrichment); write.csv(rv$te_enrichment, file, row.names = FALSE) }
  )
  output$download_te_volcano <- download_plot_server(reactive(rv$te_volcano_plot), reactive(input$format_te_volcano), "TEG_volcano", reactive(input$teg_plot_width), reactive(input$teg_plot_height))
  output$download_te_volcano_table <- downloadHandler(
    filename = function() paste0("TEG_volcano_table_", Sys.Date(), ".csv"),
    content = function(file) { req(rv$te_volcano); write.csv(rv$te_volcano, file, row.names = FALSE) }
  )

  # ---- Gene Groups ---------------------------------------------------------
  # Phase 1: "Build" button → download reference data + GO terms (one-time, slow)
  observeEvent(input$run_gene_groups, {
    req(rv$de)
    withProgress(message = "Loading reference data from GitHub + GO terms...", value = 0.1, {
      tryCatch({
        incProgress(0.8, detail = "Downloading gene-set lists and GO offspring terms")
        ctx <- setup_gene_groups(rv$de, alpha = input$alpha, lfc_cutoff = input$lfc_cutoff)
        rv$gg_context   <- ctx
        rv$gg_cache     <- list()   # clear per-group cache when re-built
        updateSelectInput(session, "gene_group_select",
                          choices  = GENE_GROUP_NAMES,
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
    withProgress(message = paste("Computing:", gsub("_", " ", gn)), value = 0.5, {
      tryCatch({
        result <- compute_single_group(gn, rv$gg_context,
                                       color_up = input$color_up %||% "#B2182B", color_down = input$color_down %||% "#2166AC",
                                       color_ns = input$color_ns %||% "#B3B3B3")
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
    selectInput("gene_group_select", "Select gene group",
                choices  = GENE_GROUP_NAMES,
                selected = GENE_GROUP_NAMES[1])
  })

  gene_group_plot_reactive <- reactive({
    req(rv$gg_cache, input$gene_group_select)
    rv$gg_cache[[input$gene_group_select]]$plot
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
    datatable(d, rownames = FALSE, filter = "top",
              options = list(pageLength = 15, scrollX = TRUE))
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

  # ---- KEGG Analysis -------------------------------------------------------
  observeEvent(input$run_kegg, {
    req(rv$de)
    withProgress(message = "Running KEGG Enrichment...", value = 0.3, {
      tryCatch({
        res <- run_kegg_enrichment(rv$de, pvalue_cutoff = input$alpha)
        rv$kegg_enrichment <- res
        
        incProgress(0.5, detail = "Generating Plot")
        rv$kegg_bubble <- plot_kegg_bubble(res, p_value_threshold = input$alpha)
        
        append_log("KEGG enrichment completed.")
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
    req(rv$kegg_bubble)
    rv$kegg_bubble
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
  
  observeEvent(input$run_pathview, {
    req(rv$de, input$pathview_select)
    pwid <- input$pathview_select
    
    withProgress(message = "Generating Pathview map (downloading from KEGG)...", value = 0.5, {
      tryCatch({
        # Filter for significant genes only so non-sig genes remain uncolored
        sig_df <- rv$de[!is.na(rv$de$pValue) & rv$de$pValue < input$alpha, ]
        lfc <- sig_df$log2FoldChange
        names(lfc) <- sig_df$gene_id
        lfc <- lfc[!is.na(lfc)]
        
        pv_dir <- file.path(tempdir(), paste0("pathview_", pwid))
        dir.create(pv_dir, showWarnings = FALSE)
        
        old_wd <- setwd(pv_dir)
        on.exit(setwd(old_wd))
        
        # Load pathview explicitly to make internal datasets (like 'bods') available
        library(pathview)
        data("bods", package = "pathview")
        
        pv_res <- pathview(
          gene.data = lfc,
          pathway.id = pwid,
          species = "ath",
          limit = list(gene = input$lfc_cutoff, cpds = 1),
          low = list(gene = "blue", cpds = "blue"),
          mid = list(gene = "gray", cpds = "gray"),
          high = list(gene = "red", cpds = "red"),
          kegg.native = TRUE,
          gene.idtype = "KEGG"
        )
        
        img_path <- file.path(pv_dir, paste0(pwid, ".pathview.png"))
        if (file.exists(img_path)) {
          rv$pathview_pathway <- img_path
          append_log("Pathview generated:", pwid)
        } else {
          showNotification("Pathview did not generate an image.", type = "warning")
        }
      }, error = function(e) {
        showNotification(paste("Pathview Error:", e$message), type = "error", duration = 15)
        append_log("Pathview error:", e$message)
      })
    })
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
  
  output$download_pathview <- downloadHandler(
    filename = function() paste0("pathview_", basename(rv$pathview_pathway)),
    content = function(file) {
      req(rv$pathview_pathway)
      file.copy(rv$pathview_pathway, file)
    }
  )
}



shinyApp(ui, server)
