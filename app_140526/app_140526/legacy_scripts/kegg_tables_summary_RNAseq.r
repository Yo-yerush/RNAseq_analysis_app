# Script to analyze KEGG pathway gene expression data
library(tidyverse)
library(knitr)
library(stringr)
library(openxlsx)

# Function to process a single directory
process_directory <- function(dir_name) {
    # Get the base path
    base_path <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/EC_kegg_maps/by_DEseq2/with_EC"
    full_path <- file.path(base_path, dir_name)

    # Check if directory exists
    if (!dir.exists(full_path)) {
        warning(paste("Directory does not exist:", full_path))
        return(NULL)
    }

    # Get all CSV files in the directory
    csv_files <- list.files(full_path, pattern = "*.csv", full.names = TRUE)

    if (length(csv_files) == 0) {
        warning(paste("No CSV files found in:", full_path))
        return(NULL)
    }

    # Initialize results dataframe
    results <- data.frame(
        Condition = character(),
        Pathway_ID = character(),
        Pathway = character(),
        Total_Genes = integer(),
        Upregulated = integer(),
        Downregulated = integer(),
        Total_Significants = integer(),
        Percentage = numeric(),
        stringsAsFactors = FALSE
    )

    # Process each CSV file
    for (file in csv_files) {
        # Extract raw pathway name from filename
        raw_pathway_name <- basename(file) %>% tools::file_path_sans_ext()

        # Extract pathway ID using regex (format like 'ath00270')
        pathway_id <- str_extract(raw_pathway_name, "ath+\\d+")

        # Clean up the pathway name by removing prefixes and suffixes
        clean_pathway <- str_replace(raw_pathway_name, paste0(dir_name, "_values_"), "")
        clean_pathway <- str_replace(clean_pathway, paste0("_", pathway_id), "")
        clean_pathway <- str_replace(clean_pathway, "_with_EC_by_DEseq2", "")

        # Read the CSV file
        tryCatch(
            {
                data <- read.csv(file)

                # Count total number of genes
                total_genes <- nrow(data)

                # Count up/down-regulated and significant genes
                sig_genes <- sum(data$pValue < 0.05 & !is.na(data$pValue))
                up_genes <- sum(data$log2FC > 0 & data$pValue < 0.05 & !is.na(data$pValue))
                down_genes <- sum(data$log2FC < 0 & data$pValue < 0.05 & !is.na(data$pValue))

                # Calculate percentage
                percentage <- ifelse(total_genes > 0, round(sig_genes / total_genes * 100, 1), NA)

                # Add to results
                results <- rbind(results, data.frame(
                    Condition = dir_name,
                    Pathway_ID = pathway_id,
                    Pathway = clean_pathway,
                    Total_Annotated = total_genes,
                    Upregulated = up_genes,
                    Downregulated = down_genes,
                    Total_Significants = sig_genes,
                    Percentage = percentage
                ))
            },
            error = function(e) {
                warning(paste("Error processing file:", file, "-", e$message))
            }
        )
    }

    return(results)
}

# List of conditions to process
conditions <- c("mto1", "mto3", "dCGS", "SSE_high", "SSE_low", "high_vs_low")
conditions_2 <- c("mto1_vs_wt", "mto3_vs_wt", "dCGS_vs_EV", "SSE_high_vs_EV", "SSE_low_vs_EV", "SSE_high_vs_SSE_low")

# Process all conditions and combine results
all_results <- data.frame()
for (cond in conditions) {
    results <- process_directory(cond)
    if (!is.null(results)) {
        all_results <- rbind(all_results, results)
    }
}

# View the results
if (nrow(all_results) > 0) {
    # Sort by condition and pathway
    # all_results <- all_results %>% arrange(Condition, Pathway_ID)

    # Display the table
    kable(all_results, caption = "KEGG Pathway Gene Expression Analysis")
} else {
    print("No results were found. Please check the directory paths and file formats.")
}

# Write results to Excel with each condition on a separate sheet
uni_cond <- unique(all_results$Condition)
if (any(uni_cond == conditions)) {
    # Create a new workbook
    wb <- createWorkbook()

    # Process each condition into a separate sheet
    for (cond in 1:length(uni_cond)) {
        # Filter for the current condition
        condition_data <- all_results %>%
            filter(Condition == uni_cond[cond]) %>%
            mutate(Condition = conditions_2[cond]) %>%
        arrange(desc(Percentage)) # Order by Percentage (descending)

        # Add a worksheet for this condition
        addWorksheet(wb, uni_cond[cond])

        # Write the data to the worksheet
        writeData(wb, uni_cond[cond], condition_data)

        # Auto-size columns for better readability
        setColWidths(wb, uni_cond[cond], cols = 1:ncol(condition_data), widths = "auto")
    }

    # Save the workbook
    saveWorkbook(wb, "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/EC_kegg_maps/by_DEseq2/with_EC/kegg_pathway_summary.xlsx", overwrite = TRUE)
}
