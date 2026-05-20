output_path_0 <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/NGS_merged_results/KEGG_combined_pathView_plots/"
output_path_rna <- paste0(output_path_0, "RNAseq/")
output_path_meth <- paste0(output_path_0, "methylome/")
dir.create(output_path_0, showWarnings = FALSE)
dir.create(output_path_rna, showWarnings = FALSE)
dir.create(output_path_meth, showWarnings = FALSE)

source("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/scripts/methylome/kegg_maps_pathview_combined_function.r")

################################################
### RNAseq
comparisons_list <- list(
    c("SSE_high", "SSE_low"),
    c("mto1", "mto3"),
    c("mto1", "dCGS")
)

base_path_0 <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/rnaseq_23/EC_kegg_maps/by_DEseq2/with_EC"

for (i.comparisons in 1:length(comparisons_list)) {
    base_path <- paste0(base_path_0, "/", comparisons_list[[i.comparisons]])
    create_pathway_comparison_pdf(
        base_path = base_path,
        comparisons = comparisons_list[[i.comparisons]],
        output_path = output_path_rna
    )
}


################################################
### methylome
comparisons_list <- list(
    c("SSE_high_vs_EV", "SSE_low_vs_EV"),
    c("mto1_vs_wt", "mto3_vs_wt"),
    c("mto1_vs_wt", "dCGS_vs_EV")
)

base_path_0 <- "C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/methionine/methylome_23/BSseq_results/KEGG_pathway/with_EC"

for (i.comparisons in 1:length(comparisons_list)) {
    for (ann in c("Genes", "Promoters")) {
        base_path <- paste0(base_path_0, "/", comparisons_list[[i.comparisons]], "/", ann, "/all")

        output_path_ann <- paste0(output_path_meth, ann, "/")
        dir.create(output_path_ann, showWarnings = FALSE)

        create_pathway_comparison_pdf(
            base_path = base_path,
            comparisons = comparisons_list[[i.comparisons]],
            output_path = output_path_ann
        )
    }
}
