# Load required libraries
library(grid)
library(png)

# Function to create pathway comparison PDFs
create_pathway_comparison_pdf <- function(base_path, comparisons, output_path) {
  # change for methylome file names (will work only if needed)
  comparisons_2 <- gsub("_vs_EV|_vs_wt", "", comparisons)

  # Validate inputs
  if (!any(dir.exists(base_path))) {
    stop("Base path does not exist: ", base_path)
  }

  if (length(comparisons) < 2) {
    stop("Need at least 2 comparison folders to compare")
  }

  # Extract pathway IDs from filenames
  extract_pathway_id <- function(filename) {
    base_name <- basename(filename)
    matches <- regexpr("ath\\d{5}", base_name)
    if (matches != -1) {
      id <- regmatches(base_name, matches)[1]
      return(id)
    } else {
      return(NA)
    }
  }

  # Extract pathway name from filename
  extract_pathway_name <- function(filename) {
    base_name <- basename(filename)

    # Check if the filename matches our expected pattern
    if (grepl("_ath\\d{5}_", base_name)) {
      # Find position after "SSE_condition_"
      start_pos <- regexpr(paste0(paste(comparisons_2, collapse = "|"), "_"), base_name)[1]
      start_pos <- start_pos + attr(regexpr(paste0(paste(comparisons_2, collapse = "|"), "_"), base_name), "match.length")

      # Find position of "_ath" that follows
      end_pos <- regexpr("_ath\\d{5}_", base_name)[1]

      # Extract the pathway name
      if (start_pos > 0 && end_pos > start_pos) {
        return(substr(base_name, start_pos, end_pos - 1))
      }
    }

    # Fallback: just return the ID
    id_match <- regexpr("ath\\d{5}", base_name)
    if (id_match != -1) {
      return(regmatches(base_name, id_match)[1])
    }

    return("Unknown pathway")
  }

  # Prepare to collect pathway files and IDs from each comparison folder
  all_files <- list()
  all_ids <- list()

  # Collect files from each comparison folder
  for (comp in comparisons_2) {
    # Construct full path to the pathview folder
    full_path <- file.path(base_path, "pathview")

    if (!any(dir.exists(full_path))) {
      warning("Path does not exist, skipping: ", full_path)
      next
    }

    # Get PNG files
    png_files <- list.files(path = full_path, pattern = "\\.png$", full.names = TRUE)

    if (length(png_files) == 0) {
      warning("No PNG files found in: ", full_path)
      next
    }

    # Store files and extract IDs
    all_files[[comp]] <- png_files[grep(comp, png_files)]
    all_ids[[comp]] <- sapply(png_files[grep(comp, png_files)], extract_pathway_id)
  }

  # Find common pathway IDs across all comparison folders
  if (length(all_ids) < 2) {
    stop("Could not find at least 2 valid comparison folders")
  }

  common_ids <- Reduce(intersect, all_ids)

  if (length(common_ids) == 0) {
    stop("No common pathway IDs found across the comparison folders")
  }

  # Create mappings from pathway ID to file path for each comparison
  file_maps <- list()
  for (comp in names(all_files)) {
    # file_maps[[comp]] <- structure(all_files[[comp]], names = all_ids[[comp]]) #################
    file_maps[[comp]] <- structure(all_files[[comp]], names = all_ids[[comp]])
  }

  # Create a PDF file
  output_prefix <- paste0(gsub("_vs_.*", "", comparisons), collapse = "_vs_")
  pdf_filename <- paste0(output_prefix, "_pathways.pdf")
  pdf(paste0(output_path, pdf_filename), width = 16, height = 9, family = "serif")

  # Loop through each common pathway ID
  for (id in common_ids) {
    # Get the pathway name from the first available file
    first_comp <- names(file_maps)[1]
    pathway_name <- extract_pathway_name(file_maps[[first_comp]][id])

    # Create a new page
    grid.newpage()

    # Calculate layout dimensions based on number of comparisons
    n_comps <- length(file_maps)

    # Determine grid layout (1 row for title, 1 for subtitles, 1 for images)
    # For the columns, we use as many as we have comparisons
    layout_grid <- grid.layout(3, n_comps,
      heights = unit(c(0.1, 0.05, 0.85), "npc"),
      widths = unit(rep(1 / n_comps, n_comps), "npc")
    )

    pushViewport(viewport(layout = layout_grid))

    # Draw the title spanning all columns
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:n_comps))
    grid.text(paste0("Pathway:", gsub("_", " ", pathway_name), " (", id, ")"),
      gp = gpar(fontsize = 16, fontface = "bold")
    )
    popViewport()

    # Process each comparison
    col_index <- 1
    for (comp in names(file_maps)) {
      # Get the file for this comparison
      comp_file <- file_maps[[comp]][id]
      # comp_file <- grep(paste0(), file_maps) ######################################

      if (is.na(comp_file) || !file.exists(comp_file)) {
        warning(paste("Skipping", id, "for", comp, "- file not found"))
        col_index <- col_index + 1
        next
      }

      # change ouput title name (of the treatment)
      comp_name_plot <- ifelse(grepl("mto", comp), paste0(comp, "\t\tVS.\t\tWT"), paste0(comp, "\t\tVS.\t\tEV"))
      comp_name_plot <- gsub("_", "-", comp_name_plot)

      # Draw subtitle
      pushViewport(viewport(layout.pos.row = 2, layout.pos.col = col_index))
      grid.text(comp_name_plot, gp = gpar(fontsize = 15, fontface = "italic"))
      popViewport()

      # Read and draw image
      tryCatch(
        {
          img <- png::readPNG(comp_file)

          pushViewport(viewport(layout.pos.row = 3, layout.pos.col = col_index))
          grid.raster(img, width = unit(1, "npc"), height = unit(1, "npc"))
          popViewport()
        },
        error = function(e) {
          warning(paste("Error processing", id, "for", comp, ":", e$message))
        }
      )

      col_index <- col_index + 1
    }

    popViewport() # pop the layout viewport
  }

  # Close the PDF device
  dev.off()

  cat("PDF created:", pdf_filename, "\n")
}

# You can run multiple comparisons in one PDF:
# comparisons <- c("SSE_high_vs_EV", "SSE_low_vs_EV", "Another_comparison")
# create_pathway_comparison_pdf(base_path, comparisons, "multi_comparison_pathways")
