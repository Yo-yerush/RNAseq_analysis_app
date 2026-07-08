deseq_fc <- function(A.B_VS_c, # DE design. <"." for "&" (and)> <"_" for " " (space-bar)>
                     experiment, # experiment folder name. example: "peel"
                     samples.number, # samples row numbers from 'col.data' file
                     exp.treatment = NA, # treated character value from 'exp' column in 'col.data' file
                     control = NA, # control character value from 'exp' column in 'col.data' file
                     group.A = NA, # samples row number for treated samples from 'col.data' file
                     path, # path for DEseq folder
                     genes_results_path = path,
                     colData_file,
                     description_file = NA,
                     lfc_shrink = FALSE) {
  if (is.na(group.A)[1] == T) {
    if (is.na(exp.treatment) == T) {
      stop("experiment treatment name from coldata file: <exp.treatment> VS <control> (mut1 VS WT)")
    }
    if (is.na(control) == T) {
      stop("control treatment name from coldata file: <exp.treatment> VS <control> (mut1 VS WT)")
    }
    contrast <- c("exp", exp.treatment, control)
  } else {
    contrast <- c("group.for.stat", "A", "B")
  }

  #################################

  # upload libraries
  lib_packages <- c(
    "dplyr", "DESeq2", "apeglm", "ggplot2", "tximport", "RColorBrewer", "pheatmap", "ashr"
  )
  for (n.pkg in seq(lib_packages)) {
    suppressMessages(library(lib_packages[n.pkg], character.only = TRUE))
    perc_val <- (n.pkg / length(lib_packages)) * 100
    cat(paste0("\rloading libraries [", round(perc_val, 1), "%] "))
  }
  cat("\n")


  # coldata file
  # info <- read.table(paste0(path, "/coldata.", experiment, ".txt"), header = T, sep = "\t")
  info <- colData_file
  names(info) <- c("x", "sample", "exp")
  
  info$sample <- as.character(info$sample)

  if (is.na(group.A)[1] == T) {
    info <- info[samples.number, ]
  } else {
    info$group.for.stat <- "B"
    info$group.for.stat[group.A] <- "A"
    info <- info[samples.number, ]
  }

  files <- c(paste0(genes_results_path, "/", info$x, ".genes.results"))
  names(files) <- c(info$x)

  ########################## DESeq2 ######################

  txData <- tximport(files, "rsem")
  txData$length[txData$length == 0] <- 1
  if (contrast[1] == "exp") {
    dds <- DESeqDataSetFromTximport(txData, info, ~exp)
  } else {
    dds <- DESeqDataSetFromTximport(txData, info, ~group.for.stat)
  }

  keep <- rowSums(counts(dds)) >= 10 # remove low expressed genes
  dds <- dds[keep, ]

  ddsDE <- DESeq(dds)
  res <- results(ddsDE, contrast = contrast, alpha = 0.05)
  if (lfc_shrink) {
    res <- lfcShrink(ddsDE, contrast = contrast, res = res, type = "ashr")
  }
  
  samples.deseq <- data.frame(rownames(res), res$log2FoldChange, res$padj, res$pvalue)
  names(samples.deseq) <- c("gene_id", "log2FoldChange", "padj", "pValue")

  samples.deseq$gene_id <- sub("\\..*", "", samples.deseq$gene_id) 

  ################ merge names and samples by ACCESSION column, save csv file #############

  # gene_id and description
  if (!is.na(description_file)) { 
    gene_id.and.description <- read.csv(description_file)
    df_merge <- merge(samples.deseq, gene_id.and.description, by = "gene_id", all.x = T) %>%
      arrange(padj)
  }
  
  # norm counts file
  normCounts <- as.data.frame(counts(ddsDE, normalized = T)) %>%
    mutate(gene_id = row.names(.)) %>%
    dplyr::select(gene_id, everything())

  # save .csv files
  new_path_1 <- paste(path, experiment, sep = "/")
  new_path_2 <- paste(path, experiment, A.B_VS_c, sep = "/")
  dir.create(new_path_1, showWarnings = F)
  dir.create(new_path_2, showWarnings = F)

  write.csv(df_merge, paste0(new_path_2, "/all_genes_results_", A.B_VS_c, ".csv"), row.names = FALSE)
  write.csv(normCounts, paste0(new_path_2, "/norm_counts_", A.B_VS_c, ".csv"), row.names = FALSE)

  # save results summary
  sink(paste0(new_path_2, "/results_summary_", A.B_VS_c, ".txt"))
  print(summary(res))
  sink()

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  ###################### plots ###############################

  ############################################################
  ######### volcano plot
  x_max <- max(abs(df_merge$log2FoldChange))
  x_max_pos <- ifelse(x_max > 15, 15, x_max)
  y_max_pos <- max(-log10(df_merge$padj))

  source("https://raw.githubusercontent.com/Yo-yerush/general_scripts/main/scripts/volcano.R")
  volcano_plot <- plot_volcano(res, xlim_both_sides = x_max_pos, ylim_up = y_max_pos, alpha.yo = 0.4)

  png(paste0(new_path_2, "/volcano_plot_", A.B_VS_c, ".png"), width = 3.5, height = 2, units = "in", res = 300, family = "serif")
  print(volcano_plot)
  dev.off()

  ############################################################
  ######### PCA
  # Label <- info$sample
  # pca_col <- info$exp
  vsd <- varianceStabilizingTransformation(dds)
  dds.est <- estimateSizeFactors(dds) # also possible to perform custom transformation
  se <- SummarizedExperiment(log2(counts(dds.est, normalized = TRUE) + 1), colData = colData(dds)) # shifted log of normalized counts
  pca_data <- plotPCA(DESeqTransform(se), intgroup = "sample", returnData = T)
  percentVar <- round(100 * attr(pca_data, "percentVar"))

  pca_data$Genotype <- info$exp[match(pca_data$name, info$x)]

  if (!is.na(control) & !is.na(exp.treatment)) {
    pca_data <- rbind(
      pca_data[pca_data$Genotype == control, ],
      pca_data[pca_data$Genotype == exp.treatment, ]
    )
    pca_data$col[pca_data$Genotype == control] <- "gray40"
    pca_data$col[pca_data$Genotype == exp.treatment] <- "#4a833a"

  } else if (!any(is.na(group.A))) {
    pca_data <- rbind(
      pca_data[pca_data$group.for.stat == "B", ],
      pca_data[pca_data$group.for.stat == "A", ]
    )
    pca_data$col[pca_data$group.for.stat == "B"] <- "gray40"
    uni_A <- pca_data$Genotype[pca_data$group.for.stat == "A"]
    pca_data$col[pca_data$group.for.stat == "A"] <- brewer.pal(length(unique(uni_A)), "Dark2")[as.numeric(as.factor(pca_data$Genotype[pca_data$group.for.stat == "A"]))]
  }

  write.csv(pca_data, paste0(new_path_2, "/PCA_table_", A.B_VS_c, ".csv"), row.names = F)

  pca_p <- ggplot(data = pca_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Genotype), size = 3.5, alpha = 0.85) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 9),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      legend.title = element_blank()
    ) +
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
    geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") +
    scale_color_manual(values = pca_data$col, breaks = pca_data$Genotype)

  svg(paste0(new_path_2, "/PCA_plot_", A.B_VS_c, ".svg"), width = 3.25, height = 2, family = "serif")
  print(pca_p)
  dev.off()

  ############################################################
  ######### MA plot
  if (!lfc_shrink) {
    resApeT <- lfcShrink(ddsDE, contrast = contrast, res = res, type = "ashr", lfcThreshold = 1)
  } else {
    resApeT <- res
  }
  
  svg(paste0(new_path_2, "/MA_plot_", A.B_VS_c, ".svg"), width = 4.65, height = 2.5, family = "serif")

  # layout(matrix(c(1, 2), nrow = 1))#, widths = c(2,2))
  par(mar = c(3, 3, 1, 8))
  par(cex.axis = 0.75) # change axis text size
  par(mgp = c(1.25, 0.1, 0)) # Move axis text closer to axis lines
  par(tck = -0.02) # Shorter tick marks

  plotMA(resApeT,
    ylim = c(-6, 6), alpha = 0.05, main = "",
    xlab = "Mean Of Normalized Counts",
    ylab = "log2 (Fold Change)",
    colSig = "white", cex = 0.25, colLine = "#ffffff00"
  )

  # add horizontal line at y=0 with specified line width
  abline(h = 0, lwd = 1, col = "gray40")
  abline(h = 1, lwd = 1, col = "gray30", lty = 2)
  abline(h = -1, lwd = 1, col = "gray30", lty = 2)

  # add points for upregulated or downregulated genes
  with(
    subset(as.data.frame(resApeT), svalue < 0.05 & log2FoldChange > 0 & log2FoldChange < 6),
    points(baseMean, log2FoldChange, col = "#a84848", cex = 0.35, pch = 20)
  )
  with(
    subset(as.data.frame(resApeT), svalue < 0.05 & log2FoldChange < 0 & log2FoldChange > -6),
    points(baseMean, log2FoldChange, col = "#5d60ba", cex = 0.35, pch = 20)
  )

  # triangle for above/bellow 6 or -6 log2FC
  with(
    subset(as.data.frame(resApeT), svalue < 0.05 & log2FoldChange >= 6),
    points(baseMean, rep(6, length(baseMean)), col = "#a84848", cex = 0.35, pch = 24)
  )
  with(
    subset(as.data.frame(resApeT), svalue < 0.05 & log2FoldChange <= -6),
    points(baseMean, rep(-6, length(baseMean)), col = "#5d60ba", cex = 0.35, pch = 25)
  )

  par(cex.axis = 1)

  # plot.new()
  # par(mar=c(0, 0, 1, 0))
  par(xpd = TRUE) # Allow drawing outside of plot region
  legend("topright",
    inset = c(-0.45, 0),
    legend = c("Upregulated", "Downregulated", "Non-DE"),
    col = c("#a84848", "#5d60ba", "gray60"),
    pch = 20,
    cex = 0.75,
    pt.cex = 1.5,
    title = "Gene Expression",
    bty = "n",
    title.font = 2
  )
  par(xpd = FALSE)
  dev.off()

  ############################################################

  message("\n\t*                           *\n\t*                           *\n\t* files saved! check errors *\n\t*                           *\n\t*                           *\n")
}
