grep_position <- function(x, df, go_bp_colname = "Gene.Ontology..biological.process.") {
    vec <- NULL
    for (terms_l in x) {
        vec <- c(vec, grep(terms_l, df[, go_bp_colname]))
    }
    return(unique(vec))
}

sig_tairs_from_terms <- function(x) {
    x0 <- x %>% distinct(gene_id)
    x1 <- x %>%
        filter(padj < 0.05) %>%
        distinct(gene_id)
    return(list(x0 = x0[, 1], x1 = x1[, 1]))
}

abiotic_stress_enrichment_test <- function(df_rna, dataset = NULL, sub_title_prefix = NULL, print_as_table = F) {
    # 'dataset' can be: 'DEGs', 'up' and 'down' (for up-/down- regulated)
    library(dplyr)
    library(ggplot2)

    source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/GO_from_tair_list.R")
    source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/GO_offspring_terms.R")

    df_rna <- df_rna %>% distinct(gene_id, .keep_all = T)
    df_sig <- df_rna %>% filter(padj < 0.05)

    if (!is.null(dataset)) {
        if (dataset == "up") {
            df_sig <- df_sig %>% filter(log2FoldChange > 0)
        } else if (dataset == "down") {
            df_sig <- df_sig %>% filter(log2FoldChange < 0)
        } else if (dataset != "up" | dataset != "down") {
            stop(error("<dataset> argument options: 'up', 'down' or keep it NULL for all significants (DEGs)"))
        }
    }

    datasetName <- ifelse(is.null(dataset), "DEGs", paste0(dataset, "DEGs"))
    sub_title <- paste0(ifelse(is.null(sub_title_prefix), "", paste0(sub_title_prefix, " - ")), datasetName)

    child_terms_cold <- offspring_fun("GO:0009409")
    child_terms_osmotic <- offspring_fun("GO:0006970")
    child_terms_salt <- offspring_fun("GO:1902074")
    child_terms_water_deprivation <- offspring_fun("GO:0009414")
    child_terms_DNA_damage <- offspring_fun("GO:0006974")
    child_terms_oxidative <- offspring_fun("GO:0006979")
    child_terms_uvb <- offspring_fun("GO:0010224")
    child_terms_wounding <- offspring_fun("GO:0009611")
    child_terms_heat <- offspring_fun("GO:0009408")

    # Test all stress terms
    stress_terms <- list(
        "Cold" = list(terms = child_terms_cold, name = "Cold stress"),
        "Osmotic" = list(terms = child_terms_osmotic, name = "Osmotic stress"),
        "Salt" = list(terms = child_terms_salt, name = "Salt stress"),
        "Water_deprivation" = list(terms = child_terms_water_deprivation, name = "Water deprivation"),
        "DNA_damage" = list(terms = `child_terms_DNA_damage`, name = "DNA damage"),
        "Oxidative" = list(terms = child_terms_oxidative, name = "Oxidative stress"),
        "UVB" = list(terms = child_terms_uvb, name = "UV-B stress"),
        "Wounding" = list(terms = child_terms_wounding, name = "Wounding"),
        "Heat" = list(terms = child_terms_heat, name = "Heat stress")
    )

    # Enrichment test function
    enrichment_test <- function(term_sig, term_total, bg_sig, bg_total, test_name = "GO term") {
        # Create contingency table
        # Rows: in term / not in term
        # Cols: significant / not significant

        a <- length(intersect(term_sig, bg_sig)) # sig genes in term
        b <- length(term_total) - a # non-sig genes in term
        c <- length(bg_sig) - a # sig genes not in term
        d <- length(bg_total) - length(term_total) - c # non-sig genes not in term

        # Contingency table
        contingency_table <- matrix(c(a, b, c, d),
            nrow = 2,
            dimnames = list(
                c("In_term", "Not_in_term"),
                c("Significant", "Not_significant")
            )
        )

        # Fisher's exact test
        fisher_result <- fisher.test(contingency_table, alternative = "greater")

        # Calculate enrichment metrics
        fold_enrichment <- (a / length(term_total)) / (length(bg_sig) / length(bg_total))

        # Results summary
        results <- list(
            test_name = test_name,
            contingency_table = contingency_table,
            sig_in_term = a,
            total_in_term = length(term_total),
            sig_in_background = length(bg_sig),
            total_in_background = length(bg_total),
            fold_enrichment = fold_enrichment,
            p_value = fisher_result$p.value,
            odds_ratio = fisher_result$estimate,
            confidence_interval = fisher_result$conf.int,
            significant = fisher_result$p.value < 0.05
        )

        return(results)
    }

    # Run enrichment for all stress terms
    enrichment_results <- list()

    for (stress_type in names(stress_terms)) {
        # Get genes for this stress type
        stress_tair <- sig_tairs_from_terms(df_rna[grep_position(stress_terms[[stress_type]]$terms, df_rna), ])

        # Run enrichment test
        enrichment_results[[stress_type]] <- enrichment_test(
            term_sig = stress_tair$x1,
            term_total = stress_tair$x0,
            bg_sig = df_sig$gene_id,
            bg_total = df_rna$gene_id,
            test_name = stress_terms[[stress_type]]$name
        )
    }

    # Summary table of all enrichments
    enrichment_summary <- data.frame(
        # Stress_type = names(enrichment_results),
        Test_name = sapply(enrichment_results, function(x) x$test_name),
        Sig_in_term = sapply(enrichment_results, function(x) x$sig_in_term),
        Total_in_term = sapply(enrichment_results, function(x) x$total_in_term),
        Fold_enrichment = round(sapply(enrichment_results, function(x) x$fold_enrichment), 2),
        P_value = sapply(enrichment_results, function(x) {
            p_val <- x$p_value
            if (p_val < 0.01) {
                format(p_val, scientific = TRUE, digits = 3)
            } else {
                sprintf("%.3f", p_val)
            }
        }),
        Odds_ratio = round(sapply(enrichment_results, function(x) x$odds_ratio), 2),
        Significant = sapply(enrichment_results, function(x) x$significant),
        stringsAsFactors = FALSE
    )

    if (print_as_table) {
        cat("=== Summary of All Stress Enrichments ===\n")
        rownames(enrichment_summary) <- seq(nrow(enrichment_summary))
        return(enrichment_summary)
    } else {
        enrichment_summary$neg_log_p <- -log10(as.numeric(gsub("e.*", "", enrichment_summary$P_value)) *
            10^as.numeric(gsub(".*e", "", enrichment_summary$P_value)))

        ggplot(enrichment_summary, aes(x = reorder(Test_name, Fold_enrichment), y = Fold_enrichment)) +
            geom_col(aes(fill = Significant), alpha = 0.85) +
            geom_text(aes(label = paste0("p=", P_value)), hjust = 1.15, size = 3) +
            geom_text(aes(label = paste(Sig_in_term, "/", Total_in_term)), hjust = -0.15, size = 3) +
            coord_flip() +
            # theme_minimal() +
            # theme(panel.grid = element_blank()) +
            ggthemes::theme_base() +
            theme(
                plot.title = element_text(size = 14, face = "bold"),
                plot.subtitle = element_text(size = 12, face = "bold"),
            ) +
            labs(
                title = "Stress Response Enrichment",
                # subtitle = bquote(italic("mto1") ~ "vs WT -" ~ .(set_title)),
                subtitle = sub_title,
                x = "Stress Response Type",
                y = "Fold Enrichment",
                fill = "Significant\n(Fisher, p<0.05)"
            ) +
            scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
            scale_y_continuous(limits = c(0, max(enrichment_summary$Fold_enrichment) * 1.15))
    }
}
