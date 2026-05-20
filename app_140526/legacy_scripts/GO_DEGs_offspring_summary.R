offspring_fun <- function(go_id_i, xx = as.list(GOBPOFFSPRING)) { # 'GOBPCHILDREN' for child terms

    child_terms_0 <- as.character(xx[[go_id_i]])
    child_terms <- c(go_id_i, child_terms_0)

    # another round of children terms
    for (i in 1:length(child_terms_0)) {
        child_terms <- c(child_terms, as.character(xx[[child_terms[i + 1]]]))
    }

    return(child_terms[!is.na(child_terms)] %>% unique()) # %>% paste(collapse = "|"))
}

################

grep_position <- function(x, rna_df) {
    vec <- NULL
    for (terms_l in x) {
        vec <- c(vec, grep(terms_l, rna_df$Gene.Ontology..biological.process.))
    }
    return(unique(vec))
}

################

GO_offspring_summary_row <- function(RNAseq_df, go_id_i, go_title = NULL) {
    go_title <- ifelse(!is.null(go_title), go_title, GOTERM[[go_id_i]]@Term)
    RNAseq_df <- RNAseq_df %>% distinct(gene_id, , .keep_all = T)

    offspring_ids <- suppressMessages(offspring_fun(go_id_i))
    go_term <- RNAseq_df[grep_position(offspring_ids, RNAseq_df), ]

    # percentages
    go_term_total <- go_term
    go_term_sig <- filter(go_term, pValue < 0.05)
    go_term_up <- filter(go_term_sig, log2FoldChange > 0)
    go_term_down <- filter(go_term_sig, log2FoldChange < 0)
    go_term_persentage <- paste0(round((nrow(go_term_sig) / nrow(go_term_total)) * 100, 2), " %")

    summary_row <- data.frame(
        Parent_GO_ID = go_id_i,
        Category = go_title,
        Total = go_term_total %>% nrow(),
        Upregulated = go_term_up %>% nrow(),
        Downregulated = go_term_down %>% nrow(),
        Significant = go_term_sig %>% nrow(),
        Percentage = go_term_persentage,
        stringsAsFactors = FALSE
    )

    cat("GO ID:\t\t", go_id_i)
    cat("\nGO Term:\t", go_title)
    cat("\nOffspring:\t", length(offspring_ids), "\n")

    summary_row
}

rename_col <- function(x, old_name, new_name) {
    names(x)[names(x) == old_name] <- new_name
    return(x)
}

GO_offspring_summary <- function(RNAseq_df, go_id, go_title = NULL, gene_id_col=NULL, pvalue_col=NULL, log2fc_col=NULL, go_bp_col=NULL) {
    suppressMessages(library(dplyr))
    suppressMessages(library(GO.db))

    if (!is.null(gene_id_col)) {rename_col(RNAseq_df, "gene_id", gene_id_col)}
    if (!is.null(pvalue_col)) {rename_col(RNAseq_df, "pValue", pvalue_col)}
    if (!is.null(log2fc_col)) {rename_col(RNAseq_df, "log2FoldChange", log2fc_col)}
    if (!is.null(go_bp_col)) {rename_col(RNAseq_df, "Gene.Ontology..biological.process.", go_bp_col)}

    final_df <- data.frame()
    for (go_id_loop in go_id) {
        # if (go_id_loop != go_id[1]) {
            cat("----------------------------\n")
        # }
        row_result <- GO_offspring_summary_row(RNAseq_df, go_id_loop, go_title)
        final_df <- rbind(final_df, row_result)
    }
    return(final_df)
}
