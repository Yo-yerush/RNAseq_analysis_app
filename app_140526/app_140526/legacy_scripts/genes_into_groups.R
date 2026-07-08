################################
### Helper utilities
################################

if (!exists("plot_theme_choice", mode = "function")) {
  plot_theme_choice <- function(plot_theme = "classic", base_size = 12, font_family = "serif") {
    if (is.null(plot_theme) || !nzchar(plot_theme)) plot_theme <- "classic"
    if (is.null(font_family) || !nzchar(font_family)) font_family <- "serif"
    switch(tolower(plot_theme),
      bw = ggplot2::theme_bw(base_size = base_size, base_family = font_family),
      minimal = ggplot2::theme_minimal(base_size = base_size, base_family = font_family),
      linedraw = ggplot2::theme_linedraw(base_size = base_size, base_family = font_family),
      light = ggplot2::theme_light(base_size = base_size, base_family = font_family),
      gray = ggplot2::theme_gray(base_size = base_size, base_family = font_family),
      dark = ggplot2::theme_dark(base_size = base_size, base_family = font_family),
      void = ggplot2::theme_void(base_size = base_size, base_family = font_family),
      classic = ggplot2::theme_classic(base_size = base_size, base_family = font_family),
      ggplot2::theme_classic(base_size = base_size, base_family = font_family)
    )
  }
}

offspring_fun <- function(go_id, xx = as.list(GO.db::GOBPOFFSPRING)) {
  child_terms_0 <- as.character(xx[[go_id]])
  child_terms   <- child_terms_0
  for (i in seq_along(child_terms_0)) {
    child_terms <- c(child_terms, as.character(xx[[child_terms[i]]]))
  }
  unique(child_terms[!is.na(child_terms)])
}

grep_position <- function(go_ids_vec, df) {
  if (!"GO.biological.process" %in% names(df) || length(go_ids_vec) == 0) return(integer(0))
  vec <- NULL
  for (term in go_ids_vec) {
    vec <- c(vec, grep(term, df$GO.biological.process, fixed = TRUE))
  }
  unique(vec)
}

safe_merge <- function(keys, data) {
  tryCatch(
    merge.data.frame(data.frame(gene_id = as.character(keys$gene_id)),
                     data, by = "gene_id"),
    error = function(e) data[0, , drop = FALSE]
  )
}

################################
### Static group names (used to populate dropdown before computation)
################################

GENE_GROUP_NAMES <- c(
  "RdDM_pathway",
  "Histone_Lysine_MTs",
  "Royal_Family_Proteins",
  "DNA_deMTs",
  "histone_deMTs",
  "REM_TFs",
  "Cohen_SSE",
  "Ash_SSE",
  "seed_specific_genes",
  "primary_metabolism",
  "secondary_metabolism",
  "methionine_biosynthesis",
  "sulfur_biosynthesis",
  "sulfur_responsive",
  "sulfur_pathway_related",
  "glucosinolate_biosynthesis",
  "hormone_signal",
  "response_to_stress",
  "response_to_biotic",
  "response_to_abiotic",
  "defense_response"
)

################################
### STEP 1 — load all reference data (run once on button click)
### Returns a "context" list stored in rv$gg_context
################################

setup_gene_groups <- function(
    de_df,
    alpha        = 0.05,
    lfc_cutoff   = 1,
    datasets_dir = "https://raw.githubusercontent.com/Yo-yerush/RA_lab_db/refs/heads/main/Arabidopsis/Groups_genes_list"
) {
  read_tab <- function(fname) {
    tryCatch(
      read.table(paste0(datasets_dir, "/", fname), sep = "\t", header = TRUE),
      error = function(e) data.frame(gene_id = character())
    )
  }
  read_csv2 <- function(fname) {
    tryCatch(
      read.csv(paste0(datasets_dir, "/", fname)),
      error = function(e) data.frame(gene_id = character())
    )
  }

  rddm <- tryCatch({
    rbind(read_tab("rddm_Matzke_et_al.txt"),
          read_tab("rddm_Cuerda_n_Slotkin.txt")) |> dplyr::distinct()
  }, error = function(e) data.frame(gene_id = character()))
  names(rddm)[1] <- "gene_id"

  HLM <- tryCatch({
    rbind(read_tab("histone_lysin_MTs_PONTVIANNE_et_al.txt"),
          read_tab("set_domain_danny_wd_et_al.txt")) |> dplyr::distinct()
  }, error = function(e) data.frame(gene_id = character()))
  names(HLM)[1] <- "gene_id"

  RF <- read_tab("At_Agenet_Tudor_family_brazil_et_al.txt")
  names(RF)[1] <- "gene_id"
  RF$gene_id <- gsub("DUF\\d+", "", RF$gene_id)

  Ash_leaf <- read_tab("cgs_ash_et_al.txt")
  if (ncol(Ash_leaf) >= 1) names(Ash_leaf)[1] <- "gene_id"

  sse <- read_tab("SSE_cohen_et_al.txt")
  if (ncol(sse) >= 1) { names(sse)[1] <- "gene_id"; sse$gene_id <- toupper(sse$gene_id) }

  primary_metabolism_0   <- read_csv2("primary_metabolism_mukherjee_et_al.csv")
  secondary_metabolism_0 <- read_csv2("secondary_metabolism_mukherjee_et_al.csv")

  gene_list_seed_specific <- tryCatch(
    read_csv2("seed_specific_TAIR_Dry_seed_n_all_Stages.csv") |>
      dplyr::rename(gene_id = tair_id),
    error = function(e) data.frame(gene_id = character())
  )

  # GO offspring terms (the slow part)
  go_offspring <- tryCatch({
    list(
      stress   = offspring_fun("GO:0006950"),
      biotic   = offspring_fun("GO:0009607"),
      abiotic  = offspring_fun("GO:0009628"),
      defence  = offspring_fun("GO:0006952")
    )
  }, error = function(e) list(stress=character(), biotic=character(), abiotic=character(), defence=character()))

  list(
    de_df                  = de_df,
    alpha                  = alpha,
    lfc_cutoff             = lfc_cutoff,
    rddm                   = rddm,
    HLM                    = HLM,
    RF                     = RF,
    Ash_leaf               = Ash_leaf,
    sse                    = sse,
    primary_metabolism_0   = primary_metabolism_0,
    secondary_metabolism_0 = secondary_metabolism_0,
    gene_list_seed_specific = gene_list_seed_specific,
    go_offspring           = go_offspring,
    sulfur_responsive_ids  = c(
      "AT3G28740","AT4G01870","AT1G34670","AT4G21990","AT4G08620","AT1G66760",
      "AT5G13580","AT1G78000","AT5G16770","AT3G27150","AT3G44320","AT4G39950",
      "AT1G47400","AT5G43780","AT3G19710","AT4G09820","AT5G24520","AT1G24100","AT3G49680"
    ),
    sulfur_pathway_ids = toupper(c(
      "At4g08620","At1g78000","At1g22150","At5g10180","At1g77990","At3g51895",
      "At4g02700","At1g23090","At3g15990","At5g19600","At5g13550","At3g12520",
      "At3g22890","At1g19920","At4g14680","At5g43780","At4g04610","At1g62180",
      "At4g21990","At2g14750","At4g39940","At3g03900","At5g67520","At5g04590",
      "At5g56760","At1g55920","At3g13110","At2g17640","At4g35640","At4g14880",
      "At3g59760","At2g43750","At3g22460","At3g03630","At3g04940","At3g61440","At5g28030"
    ))
  )
}

################################
### STEP 2 — compute a single group on demand
### Returns list(data = df, plot = ggplot2 object)
################################

make_gene_group_volcano_plot <- function(group_df, group_name, alpha = 0.05, lfc_cutoff = 1,
                                         color_up = "#B2182B", color_down = "#2166AC", color_ns = "grey70",
                                         plot_theme = "classic", font_family = "serif") {
  plot_df <- group_df[!is.na(group_df$log2FoldChange) & !is.na(group_df$padj), , drop = FALSE]

  sig_up   <- !is.na(plot_df$padj) & plot_df$padj < alpha & plot_df$log2FoldChange >= lfc_cutoff
  sig_down <- !is.na(plot_df$padj) & plot_df$padj < alpha & plot_df$log2FoldChange <= -lfc_cutoff
  plot_df$DE_class <- "not_significant"
  plot_df$DE_class[sig_up]   <- "up"
  plot_df$DE_class[sig_down] <- "down"

  plot_df$neg_log10_padj <- -log10(pmax(plot_df$padj, .Machine$double.xmin))

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = log2FoldChange, y = neg_log10_padj, color = DE_class)
  ) +
    ggplot2::geom_point(alpha = 0.8, size = 1.4) +
    ggplot2::scale_color_manual(values = c(
      up              = color_up,
      down            = color_down,
      not_significant = color_ns
    )) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
                        linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(alpha),
                        linetype = "dashed", linewidth = 0.3) +
    plot_theme_choice(plot_theme, base_size = 12, font_family = font_family) +
    ggplot2::labs(
      x     = "log2 fold change",
      y     = "-log10 adjusted p-value",
      title = gsub("_", " ", group_name),
      color = "Class"
    )
}

compute_single_group <- function(group_name, ctx, color_up = "#B2182B", color_down = "#2166AC", color_ns = "grey70",
                                 plot_theme = "classic", font_family = "serif") {
  df         <- ctx$de_df
  alpha      <- ctx$alpha
  lfc_cutoff <- ctx$lfc_cutoff

  # Ensure RNA_pvalue alias
  if (!"RNA_pvalue" %in% names(df)) {
    df$RNA_pvalue <- if ("pValue" %in% names(df)) df$pValue else df$padj
  }

  kegg_grep <- function(code) {
    if (!"KEGG_pathway" %in% names(df)) return(df[0, ])
    df[grep(code, df$KEGG_pathway), ]
  }
  col_grep <- function(pattern, col) {
    if (!col %in% names(df)) return(df[0, ])
    df[grep(pattern, df[[col]], ignore.case = TRUE), ]
  }

  group_df <- switch(group_name,

    RdDM_pathway = safe_merge(ctx$rddm, df) |> dplyr::arrange(padj),

    Histone_Lysine_MTs = {
      rbind(
        safe_merge(ctx$HLM, df),
        col_grep("SET domain|SET-domain", "short_description"),
        col_grep("atx[1-9]|atxr[1-9]|SDG[1-9]", "Symbol"),
        col_grep("atx[1-9]|atxr[1-9]|SDG[1-9]", "Protein.names"),
        col_grep("class v-like sam-binding methyltransferase", "Protein.families")
      ) |> dplyr::distinct(gene_id, .keep_all = TRUE) |> dplyr::arrange(padj)
    },

    Royal_Family_Proteins = safe_merge(ctx$RF, df) |> dplyr::arrange(padj),

    DNA_deMTs = df[grep("AT4G34060|AT5G04560|AT2G36490|AT3G10010", df$gene_id), ] |>
      dplyr::distinct(gene_id, .keep_all = TRUE) |> dplyr::arrange(padj),

    histone_deMTs = col_grep("LDL|FLD|ELF|IBM|JMJ|REF", "Symbol") |>
      dplyr::distinct(gene_id, .keep_all = TRUE) |> dplyr::arrange(padj),

    REM_TFs = df[grep(paste(c(
      "AT4G31610","AT2G24700","AT2G24690","AT2G24680","AT2G24650","AT2G24630",
      "AT4G33280","AT1G26680","AT2G46730","AT1G49480","AT4G31620","AT3G53310",
      "AT3G06220","AT3G46770","AT5G09780","AT4G31630","AT4G31640","AT4G31650",
      "AT4G31660","AT4G31680","AT4G31690","AT5G18000","AT5G60140"
    ), collapse = "|"), df$gene_id), ] |>
      dplyr::distinct(gene_id, .keep_all = TRUE) |> dplyr::arrange(padj),

    Cohen_SSE = safe_merge(ctx$sse, df) |> dplyr::arrange(padj),

    Ash_SSE = safe_merge(ctx$Ash_leaf, df) |> dplyr::arrange(padj),

    seed_specific_genes = safe_merge(ctx$gene_list_seed_specific, df) |>
      dplyr::filter(RNA_pvalue < 0.05) |> dplyr::arrange(padj),

    primary_metabolism = {
      pm <- ctx$primary_metabolism_0
      if (nrow(pm) > 0 && "gene_id" %in% names(pm))
        safe_merge(dplyr::distinct(pm, gene_id), df) |> dplyr::arrange(padj)
      else df[0, ]
    },

    secondary_metabolism = {
      sm <- ctx$secondary_metabolism_0
      if (nrow(sm) > 0 && "gene_id" %in% names(sm))
        safe_merge(dplyr::distinct(sm, gene_id), df) |> dplyr::arrange(padj)
      else df[0, ]
    },

    methionine_biosynthesis  = kegg_grep("ath00270"),
    sulfur_biosynthesis      = kegg_grep("ath00920"),
    glucosinolate_biosynthesis = kegg_grep("ath00966"),
    hormone_signal           = kegg_grep("ath04075"),

    sulfur_responsive = safe_merge(
      data.frame(gene_id = ctx$sulfur_responsive_ids), df) |> dplyr::arrange(padj),

    sulfur_pathway_related = safe_merge(
      data.frame(gene_id = ctx$sulfur_pathway_ids), df) |> dplyr::arrange(padj),

    response_to_stress = {
      idx <- grep_position(ctx$go_offspring$stress, df)
      if (length(idx) == 0) df[0, ] else df[idx, ] |> dplyr::filter(RNA_pvalue < 0.05)
    },

    response_to_biotic = {
      idx <- grep_position(ctx$go_offspring$biotic, df)
      if (length(idx) == 0) df[0, ] else df[idx, ] |> dplyr::filter(RNA_pvalue < 0.05)
    },

    response_to_abiotic = {
      idx <- grep_position(ctx$go_offspring$abiotic, df)
      if (length(idx) == 0) df[0, ] else df[idx, ] |> dplyr::filter(RNA_pvalue < 0.05)
    },

    defense_response = {
      idx <- grep_position(ctx$go_offspring$defence, df)
      if (length(idx) == 0) df[0, ] else df[idx, ] |> dplyr::filter(RNA_pvalue < 0.05)
    },

    df[0, ]  # fallback
  )

  if (is.null(group_df) || nrow(group_df) == 0) return(list(data = NULL, plot = NULL))

  # Build volcano – plot only the genes in the group, colored by DE status (same as main volcano)
  plot_df <- group_df[!is.na(group_df$log2FoldChange) & !is.na(group_df$padj), , drop = FALSE]

  sig_up   <- !is.na(plot_df$padj) & plot_df$padj < alpha & plot_df$log2FoldChange >= lfc_cutoff
  sig_down <- !is.na(plot_df$padj) & plot_df$padj < alpha & plot_df$log2FoldChange <= -lfc_cutoff
  plot_df$DE_class <- "not_significant"
  plot_df$DE_class[sig_up]   <- "up"
  plot_df$DE_class[sig_down] <- "down"

  plot_df$neg_log10_padj <- -log10(pmax(plot_df$padj, .Machine$double.xmin))

  volcano_plot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = log2FoldChange, y = neg_log10_padj, color = DE_class)
  ) +
    ggplot2::geom_point(alpha = 0.8, size = 1.4) +
    ggplot2::scale_color_manual(values = c(
      up              = color_up,
      down            = color_down,
      not_significant = color_ns
    )) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
                        linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(alpha),
                        linetype = "dashed", linewidth = 0.3) +
    plot_theme_choice(plot_theme, base_size = 12, font_family = font_family) +
    ggplot2::labs(
      x     = "log2 fold change",
      y     = "-log10 adjusted p-value",
      title = gsub("_", " ", group_name),
      color = "Class"
    )

  list(data = group_df, plot = volcano_plot)
}
