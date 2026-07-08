############################################################
# Demo preview tables for RNA-seq input file types
############################################################

count_data_preview_table <- function(input_type) {
    switch(
        input_type,
        rsem = data.frame(
            gene_id = c("AT1G01010", "AT1G01020", "AT1G01030"),
            `transcript_id(s)` = c(
                "AT1G01010.1,AT1G01010.2",
                "AT1G01020.1",
                "AT1G01030.1"
            ),
            length = c(1688, 1320, 2201),
            effective_length = c(1510.2, 1142.7, 2020.5),
            expected_count = c(350.21, 0.00, 142.88),
            TPM = c(12.45, 0.00, 5.61),
            FPKM = c(8.32, 0.00, 3.75),
            check.names = FALSE
        ),
        salmon = data.frame(
            Name = c("AT1G01010.1", "AT1G01010.2", "AT1G01020.1", "AT1G01030.1"),
            Length = c(1688, 1450, 1320, 2201),
            EffectiveLength = c(1510.2, 1272.5, 1142.7, 2020.5),
            TPM = c(8.90, 3.55, 0.00, 5.61),
            NumReads = c(250.12, 100.09, 0.00, 142.88),
            check.names = FALSE
        ),
        kallisto = data.frame(
            target_id = c("AT1G01010.1", "AT1G01010.2", "AT1G01020.1", "AT1G01030.1"),
            length = c(1688, 1450, 1320, 2201),
            eff_length = c(1510.2, 1272.5, 1142.7, 2020.5),
            est_counts = c(250.12, 100.09, 0.00, 142.88),
            tpm = c(8.90, 3.55, 0.00, 5.61),
            check.names = FALSE
        ),
        featurecounts = data.frame(
            Geneid = c("AT1G01010", "AT1G01020", "AT1G01030"),
            Chr = c("Chr1", "Chr1", "Chr1"),
            Start = c(3631, 5928, 11649),
            End = c(5899, 8737, 13714),
            Strand = c("+", "-", "+"),
            Length = c(2269, 2810, 2066),
            `sample1.bam` = c(350, 0, 143),
            `sample2.bam` = c(420, 3, 160),
            `sample3.bam` = c(900, 10, 80),
            check.names = FALSE
        ),
        countmatrix = data.frame(
            gene_id = c("AT1G01010", "AT1G01020", "AT1G01030"),
            sample1 = c(350, 0, 143),
            sample2 = c(420, 3, 160),
            sample3 = c(900, 10, 80),
            sample4 = c(870, 12, 75),
            check.names = FALSE
        ),
        count_data_preview_table("rsem")
    )
}

count_data_preview_notes <- function(input_type) {
    switch(
        input_type,
        rsem = paste(
            "RSEM - *.genes.results",
            "----------------------",
            "Load a folder with one `*.genes.results` file per sample.",
            "For RSEM transcript-level input, check the transcript-ID option and use `*.transcripts.results` with tx2gene.",
            "",
            "Example folder:",
            "RSEM_results/",
            "  sample1.genes.results",
            "  sample2.genes.results",
            "  sample3.genes.results",
            "  sample4.genes.results",
            "",
            "Main DESeq2 column: expected_count",
            sep = "\n"
        ),
        salmon = paste(
            "Salmon - */quant.sf",
            "-------------------",
            "Load the main Salmon results folder.",
            "Each sample should have its own folder containing `quant.sf`.",
            "",
            "Example folder:",
            "salmon_results/",
            "  sample1/quant.sf",
            "  sample2/quant.sf",
            "  sample3/quant.sf",
            "  sample4/quant.sf",
            "",
            "Requires tx2gene. Main DESeq2 column through tximport: NumReads",
            sep = "\n"
        ),
        kallisto = paste(
            "Kallisto - */abundance.tsv",
            "--------------------------",
            "Load the main Kallisto results folder.",
            "Each sample should have its own folder containing `abundance.tsv`.",
            "",
            "Example folder:",
            "kallisto_results/",
            "  sample1/abundance.tsv",
            "  sample2/abundance.tsv",
            "  sample3/abundance.tsv",
            "  sample4/abundance.tsv",
            "",
            "Requires tx2gene. Main DESeq2 column through tximport: est_counts",
            sep = "\n"
        ),
        featurecounts = paste(
            "featureCounts output",
            "--------------------",
            "Load one featureCounts output table.",
            "The file usually contains annotation columns followed by one count column per sample.",
            "The app uses columns after gene_biotype when present, otherwise columns after Length.",
            sep = "\n"
        ),
        countmatrix = paste(
            "Count matrix - gene x sample",
            "----------------------------",
            "Load one ready count matrix.",
            "The first column must contain gene IDs.",
            "All other columns must contain raw sample counts, with sample names as column headers.",
            sep = "\n"
        ),
        count_data_preview_notes("rsem")
    )
}
