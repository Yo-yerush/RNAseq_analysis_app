if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

run_all_safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!nzchar(x)) "run_all"
  x
}

run_all_timestamp <- function() {
  format(Sys.time(), "%d%m%y_%H%M")
}

run_all_comparison_label <- function(comparison_name = NULL) {
  comparison_name <- trimws(as.character(comparison_name %||% ""))
  if (!nzchar(comparison_name)) return("")
  run_all_safe_filename(comparison_name)
}

run_all_write_csv <- function(x, file) {
  if (is.null(x)) stop("No table is available.")
  x <- as.data.frame(x, check.names = FALSE)
  utils::write.csv(x, file, row.names = FALSE)
  normalizePath(file, winslash = "/", mustWork = FALSE)
}

run_all_write_parameters <- function(metadata, file) {
  if (is.null(metadata) || length(metadata) == 0) return("")
  values <- vapply(metadata, function(x) {
    if (is.null(x)) return("")
    x <- as.character(unlist(x, use.names = FALSE))
    x <- x[!is.na(x)]
    paste(x, collapse = ", ")
  }, character(1))
  out <- data.frame(
    Parameter = names(values),
    Value = unname(values),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  utils::write.csv(out, file, row.names = FALSE)
  normalizePath(file, winslash = "/", mustWork = FALSE)
}

run_all_plot_file <- function(output_dir, filename_base, format = "png") {
  format <- tolower(format %||% "png")
  if (!format %in% c("png", "svg", "pdf")) format <- "png"
  file.path(output_dir, paste0(run_all_safe_filename(filename_base), ".", format))
}

run_all_unique_file_path <- function(path) {
  if (!file.exists(path)) return(path)
  dir <- dirname(path)
  ext <- tools::file_ext(path)
  stem <- tools::file_path_sans_ext(basename(path))
  for (i in seq_len(999)) {
    candidate <- file.path(dir, paste0(stem, "_", i, if (nzchar(ext)) paste0(".", ext) else ""))
    if (!file.exists(candidate)) return(candidate)
  }
  path
}

run_all_prefix_output_files <- function(files, prefix = "") {
  prefix <- run_all_comparison_label(prefix)
  if (!nzchar(prefix)) return(files)
  vapply(files, function(file) {
    if (is.na(file) || !nzchar(file) || !file.exists(file) || dir.exists(file)) return(file)
    base <- basename(file)
    if (startsWith(base, paste0(prefix, "_"))) return(normalizePath(file, winslash = "/", mustWork = FALSE))
    target <- run_all_unique_file_path(file.path(dirname(file), paste0(prefix, "_", base)))
    ok <- tryCatch(file.rename(file, target), error = function(e) FALSE)
    normalizePath(if (isTRUE(ok)) target else file, winslash = "/", mustWork = FALSE)
  }, character(1), USE.NAMES = FALSE)
}

run_all_save_plot <- function(plot_obj, file, width_px = 700, height_px = 450, dpi = 96) {
  if (is.null(plot_obj)) stop("No plot is available.")
  width_in <- max(suppressWarnings(as.numeric(width_px)), 100, na.rm = TRUE) / dpi
  height_in <- max(suppressWarnings(as.numeric(height_px)), 100, na.rm = TRUE) / dpi
  ggplot2::ggsave(file, plot = plot_obj, width = width_in, height = height_in, units = "in", dpi = dpi, bg = "white")
  normalizePath(file, winslash = "/", mustWork = FALSE)
}

run_all_log_line <- function(lines, status, label, message = "") {
  st_sep <- if (nchar(status) == 4) {
    "   "
  } else if (nchar(status) == 5) {
    "  "
  } else if (nchar(status) == 6) {
    " "
  } else {
    NULL
  }

  return(c(lines, paste(format(Sys.time(), "%d-%m-%Y %H-%M-%S"), "", status, st_sep, label, message)))
}

run_all_task_folder <- function(id, task = NULL) {
  explicit <- task$folder %||% ""
  if (nzchar(explicit)) return(run_all_safe_filename(explicit))
  if (id %in% c("de_table", "de_summary", "normalized_counts", "volcano", "ma", "pca")) return("Core_outputs")
  if (grepl("^(go_|revigo_|go_offspring$|go_stress_)", id)) return("GO")
  if (id %in% c("kegg", "pathview")) return("KEGG")
  if (grepl("^msigdb_", id)) return("MSigDB_Hallmark")
  if (grepl("^(pmn_|pmn_pathway$)", id)) return("PMN")
  if (grepl("^gene_family_", id)) return("Gene_Families")
  if (grepl("^(te_|te_overlap_)", id)) return("TE_analysis")
  run_all_safe_filename(id)
}

run_all_report_template_path <- function() {
  candidates <- c(
    file.path(getwd(), "R", "run_all_report.Rmd"),
    file.path(getwd(), "app_140526", "R", "run_all_report.Rmd")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) stop("Could not find R Markdown report template: R/run_all_report.Rmd")
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

run_all_find_pandoc <- function() {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) return("")
  existing <- tryCatch(rmarkdown::find_pandoc(), error = function(e) list(dir = "", version = ""))
  if (!is.null(existing$dir) && nzchar(existing$dir) && file.exists(file.path(existing$dir, "pandoc.exe"))) {
    return(normalizePath(existing$dir, winslash = "/", mustWork = FALSE))
  }

  env_paths <- unique(c(
    Sys.getenv("RSTUDIO_PANDOC", unset = ""),
    Sys.getenv("QUARTO_PANDOC", unset = "")
  ))
  env_paths <- env_paths[nzchar(env_paths)]

  pf <- Sys.getenv("ProgramFiles", unset = "C:/Program Files")
  pfx86 <- Sys.getenv("ProgramFiles(x86)", unset = "C:/Program Files (x86)")
  local_app <- Sys.getenv("LOCALAPPDATA", unset = "")

  candidates <- unique(c(
    env_paths,
    file.path(pf, "RStudio", "resources", "app", "bin", "quarto", "bin", "tools"),
    file.path(pf, "RStudio", "resources", "app", "bin", "pandoc"),
    file.path(pf, "RStudio", "bin", "pandoc"),
    file.path(pfx86, "RStudio", "bin", "pandoc"),
    file.path(pf, "Quarto", "bin", "tools"),
    file.path(local_app, "Programs", "Quarto", "bin", "tools")
  ))
  candidates <- candidates[nzchar(candidates)]
  candidates <- candidates[file.exists(file.path(candidates, "pandoc.exe"))]
  if (length(candidates) == 0) return("")
  normalizePath(candidates[1], winslash = "/", mustWork = FALSE)
}

run_all_prepare_pandoc <- function() {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) return(FALSE)
  pandoc_dir <- run_all_find_pandoc()
  if (!nzchar(pandoc_dir)) return(FALSE)
  Sys.setenv(RSTUDIO_PANDOC = pandoc_dir)
  if (requireNamespace("rmarkdown", quietly = TRUE) && exists("find_pandoc", envir = asNamespace("rmarkdown"), mode = "function")) {
    tryCatch(rmarkdown::find_pandoc(cache = FALSE), error = function(e) NULL)
  }
  isTRUE(rmarkdown::pandoc_available())
}

run_all_report_available <- function() {
  requireNamespace("rmarkdown", quietly = TRUE) &&
    requireNamespace("knitr", quietly = TRUE) &&
    isTRUE(run_all_prepare_pandoc())
}

run_all_render_html_report <- function(run_dir, template = run_all_report_template_path(), comparison_name = NULL) {
  if (!run_all_report_available()) {
    stop("R Markdown report dependencies are not available. Install rmarkdown/knitr and make sure Pandoc is available.")
  }
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  comparison_label <- run_all_comparison_label(comparison_name)
  report_base <- if (nzchar(comparison_label)) paste0("RNAseq_Run_All_Report_", comparison_label) else "RNAseq_Run_All_Report"
  report_file <- file.path(run_dir, paste0(report_base, ".html"))
  rmarkdown::render(
    input = template,
    output_file = basename(report_file),
    output_dir = run_dir,
    params = list(run_dir = run_dir, comparison_name = comparison_name %||% ""),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )
  normalizePath(report_file, winslash = "/", mustWork = FALSE)
}

run_all_execute <- function(selected, output_dir, tasks, append_log = NULL, render_report = TRUE, metadata = NULL, comparison_name = NULL) {
  if (is.null(selected) || length(selected) == 0) stop("Select at least one analysis to run.")
  if (is.null(output_dir) || !nzchar(output_dir)) stop("Choose an output folder.")
  if (!dir.exists(output_dir)) stop("Output folder does not exist: ", output_dir)

  comparison_label <- run_all_comparison_label(comparison_name)
  run_prefix <- if (nzchar(comparison_label)) paste0("run_all_", comparison_label, "_") else "run_all_"
  run_dir <- file.path(output_dir, paste0(run_prefix, run_all_timestamp()))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(run_dir)) stop("Could not create output folder: ", run_dir)

  log_file <- file.path(run_dir, "run_all_log.txt")
  lines <- character()
  lines <- run_all_log_line(lines, "[START]", "Run All", paste("Output:", normalizePath(run_dir, winslash = "/", mustWork = FALSE)))
  writeLines(lines, log_file)

  if (!is.null(metadata) && length(metadata) > 0) {
    tryCatch(
      run_all_write_parameters(metadata, file.path(run_dir, "run_all_parameters.csv")),
      error = function(e) {
        lines <<- run_all_log_line(lines, "[WARN]", "Run parameters", e$message)
        writeLines(lines, log_file)
      }
    )
  }

  results <- list()
  for (id in selected) {
    task <- tasks[[id]]
    label <- if (!is.null(task$label)) task$label else id
    if (is.null(task) || !is.function(task$run)) {
      lines <- run_all_log_line(lines, "[SKIP]", id, "Task is not available in this app session.")
      writeLines(lines, log_file)
      next
    }

    if (is.function(append_log)) append_log("Run All:", label)
    lines <- run_all_log_line(lines, "[RUN]", label)
    writeLines(lines, log_file)

    result <- tryCatch({
      task_dir <- file.path(run_dir, run_all_task_folder(id, task))
      dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(task_dir)) stop("Could not create task output folder: ", task_dir)
      files <- task$run(task_dir)
      files <- files[!is.na(files) & nzchar(files)]
      files <- run_all_prefix_output_files(files, comparison_label)
      rel_files <- file.path(basename(dirname(files)), basename(files))
      msg <- if (length(files) > 0) paste("Files:", paste(rel_files, collapse = ", ")) else "Completed; no files produced."
      lines <- run_all_log_line(lines, "[OK]", label, msg)
      list(ok = TRUE, files = files, error = NULL)
    }, error = function(e) {
      lines <<- run_all_log_line(lines, "[ERROR]", label, e$message)
      list(ok = FALSE, files = character(), error = e$message)
    })

    results[[id]] <- result
    writeLines(lines, log_file)
  }

  report_file <- ""
  if (isTRUE(render_report)) {
    if (run_all_report_available()) {
      lines <- run_all_log_line(lines, "[RUN]", "HTML report")
      writeLines(lines, log_file)
      report_result <- tryCatch({
        report_file <- run_all_render_html_report(run_dir, comparison_name = comparison_name)
        lines <- run_all_log_line(lines, "[OK]", "HTML report", paste("File:", basename(report_file)))
        list(ok = TRUE, files = report_file, error = NULL)
      }, error = function(e) {
        lines <<- run_all_log_line(lines, "[ERROR]", "HTML report", e$message)
        list(ok = FALSE, files = character(), error = e$message)
      })
      results[["html_report"]] <- report_result
      writeLines(lines, log_file)
    } else {
      lines <- run_all_log_line(lines, "[SKIP]", "HTML report", "Install rmarkdown/knitr and make sure Pandoc is available to render the report.")
      writeLines(lines, log_file)
    }
  }

  lines <- run_all_log_line(lines, "[DONE]", "Run All", paste("Log:", normalizePath(log_file, winslash = "/", mustWork = FALSE)))
  writeLines(lines, log_file)

  list(
    output_dir = normalizePath(run_dir, winslash = "/", mustWork = FALSE),
    log_file = normalizePath(log_file, winslash = "/", mustWork = FALSE),
    report_file = if (nzchar(report_file)) normalizePath(report_file, winslash = "/", mustWork = FALSE) else "",
    log = lines,
    results = results
  )
}
