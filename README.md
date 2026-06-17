# RNA-seq Analysis Dashboard

A local R Shiny dashboard for RNA-seq differential expression analysis, count/quantification import, annotation building, dependency checks, and downstream visualization.

**Author:** Yonatan Yerushalmy  
Plant's metabolism and molecular genetic lab, Prof. Rachel Amir group

**Repository:**  
[https://github.com/Yo-yerush/RNAseq_analysis_app.git](https://github.com/Yo-yerush/RNAseq_analysis_app.git)

## How To Run On Windows

1. Install R for Windows, preferably R 4.2 or newer.
2. Extract the project folder.
3. First time: double-click `install.bat`.
4. Later runs: double-click `RNAseq_analysis_app.bat`.

If `Rscript.exe was not found in PATH`, add your R `bin` folder to PATH, usually:

```text
C:\Program Files\R\R-4.x.x\bin
```

## Data Input

### Upload DE Results

Upload CSV, TSV, TXT, XLS, or XLSX differential expression results. The app auto-detects comma and tab delimiters for text files.

Must include columns:

| Column | Description |
|--------|-------------|
| `gene_id` | Gene identifier, for example TAIR, Ensembl, Entrez, or symbol |
| `log2FoldChange` | log2 fold change |
| `padj` | adjusted p-value |

If standard column names are not found, the app treats the first three columns as `gene_id`, `log2FoldChange`, and `padj`.

After loading the table, the app auto-detects common `gene_id` formats such as TAIR, Ensembl, RefSeq, Entrez, UniProt, and gene symbols, then updates the Gene ID type in the Organism annotations tab. You can still change it manually.

### Run DESeq2 From Quantification Or Count Files

The app can run DESeq2 directly from several count/quantification input types:

| Input type | Expected files | Notes |
|------------|----------------|-------|
| RSEM gene mode | `*.genes.results` | Default RSEM mode. Input rows are treated as gene IDs. |
| RSEM transcript mode | `*.transcripts.results` | Check **RSEM files contain transcript IDs** and provide a two-column tx2gene table. |
| Salmon | `*/quant.sf` | One sample folder per `quant.sf`; requires tx2gene. |
| Kallisto | `*/abundance.tsv` | One sample folder per `abundance.tsv`; requires tx2gene. |
| featureCounts | one uploaded TXT/TSV/CSV table | Uses columns after `gene_biotype` when present, otherwise after `Length`. |
| Count matrix | one uploaded CSV/TSV/TXT/XLS/XLSX table | First column is gene IDs; remaining columns are raw sample counts. |

The app scans sample IDs and creates editable colData.

colData supports CSV, TSV, TXT, XLS, and XLSX files. Text files can use comma or tab delimiters.

Supported column patterns:

| Column | Column positions options | Is unique? |
|----------|--------------|--------------|
| `sample_id` | `1st` / `2nd` | Yes |
| `condition` | `1st` / `2nd` / `3rd` | No |
| `sample_label` optional | `2nd` / `3rd` / `4th` | Yes |
| `effect` optional | `≥ 3rd` | NO |

After colData is loaded, use the **Condition / group column** selector above the colData preview to choose which column defines treatment/control groups.

tx2gene options support:

- Upload a two-column tx2gene table.
- Build tx2gene from GTF/GFF by choosing transcript and gene ID attributes.
- Build Arabidopsis TAIR-style tx2gene by stripping the final `.number` suffix from transcript IDs.
- Preview and download the generated tx2gene table before using it.

Extra colData columns are preserved. You can choose one extra effect column for DESeq2:

- No extra effect: `~ condition`
- Adjusted model: `~ condition + effect`
- Interaction model: `~ condition + effect + condition:effect`

The Data tab prints the exact design formula and contrast used.

When running DESeq2, the selected treatment/control comparison still drives the main DE table and all downstream analyses. The app also extracts every condition level versus the selected control from the same DESeq2 model. The Data tab shows a Venn diagram of shared significant gene IDs across comparisons, with display controls for the Venn colors, and provides a combined all-comparisons CSV download.

In DESeq2 mode, tables with gene rows include a small plot button column. Clicking it opens a normalized-count boxplot for that gene when normalized counts are available.

## Organism And Gene ID Settings

The **Organism annotations** tab controls organism-wide settings:

- Search/select organism name or NCBI taxonomy ID.
- Choose **Gene ID type**: common OrgDb key types such as `TAIR`, `ENTREZID`, `SYMBOL`, `ALIAS`, `ENSEMBL`, `REFSEQ`, `ACCNUM`, and `UNIPROT` where supported by the selected OrgDb.
- The Gene ID type is used by annotation builders, GO, KEGG, PMN, and MSigDB/Hallmark where relevant.

Default annotation buttons are available for Arabidopsis, human, and E. coli K-12 MG1655 when those organisms are selected.

The UniProt builder scans available UniProt ID columns for the selected organism, lets you choose an ID source, and builds a description table using that source as the new `gene_id`. For human Ensembl IDs, for example `ENSG00000141510`, the builder can use the selected OrgDb, such as `org.Hs.eg.db`, to bridge Ensembl IDs to UniProt annotations.

The RefSeq GTF builder accepts NCBI-downloaded GTF files, scans GTF attributes and `db_xref` values, then lets you choose which ID source becomes the new `gene_id`. This is useful for switching from locus tags such as `b0001` to RefSeq protein IDs such as `NP_414542`.

When a builder replaces `gene_id`, the original loaded ID is kept in `original_gene_id` for traceability.

Annotation table downloads include the organism name and tax ID in the filename.

## Tabs And Features

### Data Input

- Summary of loaded genes.
- Editable colData in DESeq2 count/quantification mode.
- All condition-vs-control DESeq2 results with a shared significant-gene Venn diagram.
- PCA plot after running DESeq2 from quantification/count files.
- PCA condition selector and sample-label toggle.
- DE table preview.
- Download DE table and normalized counts.

### Organism Annotations

- Load a manual annotation CSV/TSV/TXT.
- Build annotation table from UniProt.
- Build annotation table from RefSeq/NCBI GTF.
- Select organism and Gene ID type.
- Preview and download current annotation source.

### DE Results

- Volcano plot.
- MA plot, requires `baseMean`.
- Full-text search across annotation columns.
- Supports AND-style space-separated terms and OR with `|`.

### GO Analysis

- topGO enrichment with `weight01`, `classic`, or `elim`.
- Fisher or KS statistic.
- **GO genes** sub-tab: choose one or more GO IDs and view matched loaded genes as a volcano plot and table.
- REVIGO-like local semantic reduction using `rrvgo`.
- GO offspring summaries for custom parent GO terms.
- Abiotic stress (plants) GO enrichment.

GO display cutoff, ontology, and top-N controls appear in the sidebar only while on the GO tab.

### KEGG Analysis

- KEGG enrichment using `KEGGREST`.
- The app downloads KEGG pathway gene sets by KEGG organism code, for example `ath` or `hsa`, and caches them locally.
- For human Ensembl data, KEGG uses Entrez/numeric IDs, so the app maps selected Gene ID type to Entrez IDs through the selected OrgDb before matching pathways.
- For E. coli `eco`, KEGG pathway genes are usually `b####` locus tags. RefSeq protein IDs such as `NP_414542` can work for GO with `REFSEQ`, but they do not directly match the current KEGG `eco` pathway gene IDs unless mapped back to KEGG-compatible locus IDs.
- Bubble plot and enrichment table.
- Pathview pathway maps colored by log2FoldChange. Pathview also uses the same ID mapping, so Ensembl human genes can be colored correctly after mapping.

### MSigDB/Hallmark

- Hallmark over-representation analysis using `msigdbr`.
- Supports up, down, or all significant genes.
- Uses the selected organism/species and Gene ID type where available.
- **Hallmark genes** sub-tab: choose one or more Hallmark set codes and view matched loaded genes as a volcano plot and table.
- Bubble plot and downloadable results table.

`msigdbr` may download/cache MSigDB data on first use.

### PMN (plants)

- Plant Metabolic Network pathway enrichment for plant Cyc databases, for example `AraCyc`, `OryzaCyc`, `CornCyc`, and `TomatoCyc`.
- The PMN Cyc DB is auto-selected from the organism in the Organism annotations tab when a known plant mapping exists.
- Downloads PMN tab-delimited pathway tables from the public PMN file bucket when PMN analysis or pathway lookup is run.
- PMN matching usually expects the organism locus IDs used by that Cyc database, such as TAIR locus IDs for AraCyc.
- Bubble plot, downloadable enrichment table, pathway-gene lookup table, and selected-pathway volcano plot.

### Genes Groups

- **Gene families - enrichment** sub-tab: tests significant genes against families and shows a top-10 dotplot plus a downloadable table.
- **Gene families** sub-tab: selected-family volcano plots and tables. Arabidopsis uses the RA lab database file:
  - `https://raw.githubusercontent.com/Yo-yerush/RA_lab_db/refs/heads/main/description_files/gene_families_sep_29_09_update.txt`
- Human uses HGNC `gene_has_family.csv` plus `family.csv`; `family.csv$name` is used as the family name. Input gene IDs are mapped to `hgnc_id`, while the original `gene_id` is kept in the preview table.
- Arabidopsis matching uses `Genomic_Locus_Tag` as uppercase `gene_id`, `Gene_Family` for selectable families, and shows `Sub_Family` in the preview table.
- **Custom groups - RA lab (At)** sub-tab: curated built-in gene sets with group-specific volcano plots and tables.
- Reference data is downloaded from GitHub on first use and cached in-session.

### TE Analysis

- Arabidopsis-specific TE workflows for transposable-element genes and DE genes overlapping nearby/gene-body TEs.
- **TEGs Enrichment Analysis**: TE superfamily enrichment for transposable-element genes.
- **TEGs Volcano Plot**: volcano plot for selected TE superfamilies.
- **Overlapped TEs**: finds DE genes whose upstream, downstream, both-side, or gene-body ranges overlap TAIR10 TEs.
- The Overlapped TEs tab shows:
  - overlapped-gene volcano, filterable by TE family
  - TE family count plot and downloadable family-count table
  - TE family enrichment plot and downloadable enrichment table
  - gene table with overlapped TE IDs, families, superfamilies, DE columns, and available annotation columns
- TE family enrichment uses Fisher exact tests and can compare against either all TAIR10 TEs or region-aware TEs that overlap the same selected gene-region windows across all TAIR genes.
- Uses default Arabidopsis annotation, TAIR10 TE metadata, and TAIR gene ranges from GitHub:
  - `https://github.com/Yo-yerush/RA_lab_db/raw/refs/heads/main/description_files/At_custom_description_file.csv.gz`
  - `https://raw.githubusercontent.com/Yo-yerush/RA_lab_db/refs/heads/main/description_files/TAIR10_Transposable_Elements.txt`
  - `https://raw.githubusercontent.com/Yo-yerush/RA_lab_db/refs/heads/main/description_files/TAIR_genes_short.csv`

Human or other organism TE analysis requires a compatible TE-level annotation table or TE quantification output, such as RepeatMasker/TEtranscripts-derived families. Standard human Ensembl gene-level DE tables do not directly contain TE family assignments.

### Run All

- Batch runner for exporting selected analyses into a timestamped `run_all_YYYYMMDD_HHMMSS` folder.
- Select all available outputs, clear selections, or choose specific outputs by section:
  - Core outputs: DE table, DE summary, normalized counts, Volcano, MA, and PCA.
  - Enrichment: GO enrichment, REVIGO-like GO reduction, GO offspring summary, abiotic-stress GO enrichment, KEGG enrichment, Pathview, MSigDB/Hallmark, PMN enrichment, and PMN pathway-gene lookup.
  - TE analysis: TEG enrichment, TEG volcano, and Overlapped TEs when Arabidopsis is selected.
- Direction selector for enrichment-style analyses: upregulated only, downregulated only, up and down, all DE genes, or up/down/all.
- Unavailable analyses are hidden based on the selected organism and current app state.
- Outputs are organized into subfolders such as `Core_outputs`, `GO`, `KEGG`, `MSigDB_Hallmark`, `PMN`, and `TE_analysis`.
- Each run writes `run_all_log.txt`. If one selected analysis fails, the batch continues and the error is recorded in the log.
- Optional HTML report: when `rmarkdown`, `knitr`, and Pandoc are available, Run All creates `RNAseq_Run_All_Report.html` with plot images, CSV previews, file links, and short explanations for each analysis section.
- Gene family enrichment is included for Arabidopsis and human when available, with outputs saved under `Gene_Families`.
- The HTML report focuses on plots, analysis explanations, interpretation notes, and a parameter summary. Result tables are saved as CSV files in the output folders but are not embedded in the report.

### Log / Help

- Session log.
- Notes and app usage help.
- Parameters sub-tab showing the current analysis settings. It starts from defaults and updates when controls are changed.
- Dependencies sub-tab showing required R packages, installed status, versions, and which analyses are affected if a package is missing.

## Sidebar Controls

- Data input mode and DESeq2 count/quantification controls.
- Significance thresholds: `padj` and `|log2FC|`.
- Tab-specific plot size controls.
- GO and Hallmark filters only appear on their relevant tabs.
- Global point size, alpha, and colors.

The upregulated/downregulated/not-significant colors are used across DE result plots and enrichment visualizations where relevant, including GO, KEGG, MSigDB/Hallmark, Pathview, and gene group plots.

## Notes

- The bundled TE workflow is Arabidopsis-specific.
- GO requires the selected OrgDb package to be installed, for example `org.At.tair.db` or `org.Hs.eg.db`.
  - Some analysis tabs require extra Bioconductor packages when used, especially `GO.db`, `KEGGREST`, `pathview`, and the selected `org.*.db` organism package. The app can install supported OrgDb packages when GO is run.
- E. coli K-12 GO analysis supports `b####` locus tags through the app's alias normalization and supports RefSeq accessions through OrgDb `REFSEQ`/`ACCNUM` where available.
- KEGG and Pathview require internet access when KEGG data or pathway images are not cached.
- UniProt annotation building requires internet access.
- RefSeq GTF annotation building works from a local GTF file and does not require internet access after the GTF is downloaded.
- PCA requires count data and is available only after running DESeq2 from quantification/count files.
- Run All HTML reports require `rmarkdown`, `knitr`, and Pandoc. The installer installs the R packages, but Pandoc usually comes from RStudio, Quarto, or a standalone Pandoc installation. If Pandoc is missing, Run All still completes and skips the HTML report.
- No results are saved automatically. Use the download buttons.

## Troubleshooting R Not Found

The launcher searches for `Rscript.exe` in PATH and common Windows install locations:

- `C:\Program Files\R\...`
- `C:\Program Files (x86)\R\...`
- `%LOCALAPPDATA%\Programs\R\...`

If it still cannot find R, run `app_140526/diagnose_R_installation.bat`, then edit the launcher path manually if needed:

```bat
set "RSCRIPT=C:\Program Files\R\R-4.5.0\bin\Rscript.exe"
```

## File Structure

- `app_140526/` - main app folder:
  - `app.R` - Shiny dashboard UI/server.
  - `R/helpers.R` - DESeq2 pipeline, DE result plots, GO/MSigDB helpers, annotation utilities.
  - `R/run_all.R` - Batch Run output helpers, task execution, organized output folders, logging, and optional HTML report rendering.
  - `R/run_all_report.Rmd` - R Markdown template used to render the optional Run All HTML report.
  - `R/build_uniprot_description_file.R` - UniProt annotation builder, including Ensembl/SYMBOL/ENTREZID/TAIR ID support through OrgDb where available.
  - `R/build_refseq_gftf_description_file.R` - RefSeq GTF annotation builder for NCBI-downloaded GTF files, including selectable `db_xref`/attribute ID sources.
  - `legacy_scripts/kegg_analysis.R` - KEGG enrichment and KEGG ID mapping helpers.
  - `legacy_scripts/volcano_TEG_overlap_with_TE_families_RNAseq.R` - Arabidopsis TE superfamily enrichment and volcano helpers.
  - `install_packages.R` - installs the core CRAN and Bioconductor packages used by the app.
  - `launch_app.R` - launches the Shiny app from R.
  - `diagnose_R_installation.bat` - checks common Windows R installation paths if the launchers cannot find `Rscript.exe`.
  - `description_files/` - optional local data files. Default Arabidopsis, human, E. coli K-12 MG1655 annotation tables and TAIR10 TE metadata are loaded from GitHub when internet is available.
- `install.bat` - Windows first-run launcher that installs missing packages.
- `RNAseq_analysis_app.bat` - faster launcher for later runs.
