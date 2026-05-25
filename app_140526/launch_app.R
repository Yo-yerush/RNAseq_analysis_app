# Launch the local RNA-seq dashboard app
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
app_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
setwd(app_dir)
if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("Package 'shiny' is not installed. Run install_packages.R first.")
}
shiny::runApp(appDir = app_dir, launch.browser = TRUE)
