GO_from_tairList <- function(interesting.TAIRs, background.TAIRs = "all", treatment = "TAIR list", ontology_type = "BP", pValue.threshold = 0.01, topGO.algorithm = "weight01", topGO.statistic = "fisher", netPlot = FALSE) {
    # this function run GO analysis for arabidopsis TAIR IDs
    # just insert the gene list to test (interesting.TAIRs)
    # and the beckground TAIR list (background.TAIRs)
    #
    # if using 'netPlot=TRUE', 'background.TAIRs' argument is not needed
    #
    # GO_from_tairList(interesting.TAIRs, background.TAIRs, pValue.threshold = 0.001)
    # GO_from_tairList(interesting.TAIRs, netPlot = TRUE)

    suppressMessages({
        library(topGO)
        library(org.At.tair.db)
        library(ggplot2)
    })

    msg <- paste(
        "\n***************************\nontology type:", ontology_type, "\npadj threshold:", pValue.threshold, "\ntopGO algorithm:", topGO.algorithm, "\ntopGO statistic:", topGO.statistic,
        ifelse(
            any(background.TAIRs == "all") & !netPlot,
            "\nbackground genes: all TAIR IDs\n\n",
            ifelse(netPlot, "\nNetPlot TAIRs count:", "\n\n")
        )
    )

    # all TAIR IDs (if not inserted beckground TAIR list)
    if (any(background.TAIRs == "all")) {
        background.TAIRs <- keys(org.At.tair.db, keytype = "TAIR")
    }


    if (!netPlot) {
        ###### GO analysis table
        # indicating if a gene is interesting (1) or not (0)
        geneList <- factor(as.integer(background.TAIRs %in% interesting.TAIRs))
        names(geneList) <- background.TAIRs

        # topGOdata object
        GOdata <- new("topGOdata",
            ontology = ontology_type,
            allGenes = geneList,
            geneSelectionFun = function(x) (x == 1),
            # nodeSize = 10,
            annot = annFUN.org,
            mapping = "org.At.tair.db"
        )

        # run and output significants terms
        resultClassic <- runTest(GOdata, algorithm = topGO.algorithm, statistic = topGO.statistic)
        allResClassic <- GenTable(GOdata, pValue = resultClassic, topNodes = length(resultClassic@score))
        allResClassic_sig <- allResClassic[allResClassic$pValue <= pValue.threshold, ]

        cat(msg)
        return(allResClassic_sig)
    } else {
        ###### Net-Plot
        suppressMessages(library(clusterProfiler))

        cat("netPlot...\n")
        ego <- enrichGO(
            gene = interesting.TAIRs,
            OrgDb = org.At.tair.db,
            keyType = "TAIR",
            ont = ontology_type,
            pAdjustMethod = "BH",
            pvalueCutoff = pValue.threshold,
            qvalueCutoff = 0.2
        )
        ego_sig <- ego@result[ego@result$p.adjust < pValue.threshold, ]
        cat(msg, nrow(ego_sig), "\n\n")

        title_size <- if (nrow(ego_sig) <= 10) {
            20
        } else if (nrow(ego_sig) > 100) {
            nrow(ego_sig) / 10
        } else {
            nrow(ego_sig) * 2
        }

        if (nrow(ego_sig) != 0) {
            cnetplot(ego, showCategory = nrow(ego_sig)) +
                ggtitle(paste0(gsub("_", " ", treatment), " (enrichment p.adj<", pValue.threshold, ")")) +
                theme(plot.title = element_text(hjust = 0.5, size = title_size))
        } else {
            stop("there is no significant enrichments")
        }
    }
}
