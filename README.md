# RNA-seq Analysis Dashboard

A local R Shiny dashboard for Arabidopsis RNA-seq differential expression analysis and downstream visualization.

**Author:** Yonatan Yerushalmy  
Plant's metabolism and molecular genetic lab, Prof. Rachel Amir group

**Source Code & Repository:**  
You can download, edit, and access the raw files and source code for this project from its GitHub repository:  
[https://github.com/Yo-yerush/RNAseq_analysis_app.git](https://github.com/Yo-yerush/RNAseq_analysis_app.git)

## File structure

- `app_140526/` – contains all the R scripts and app data:
  - `app.R` – main Shiny dashboard.
  - `R/helpers.R` – core functions: DESeq2 pipeline, volcano/MA/PCA plots, GO enrichment, GO offspring summaries, abiotic stress enrichment, REVIGO-like semantic reduction.
  - `install_packages.R` – installs all required CRAN and Bioconductor packages.
  - `launch_app.R` – launches the Shiny app from R.
  - `example_data/all_genes_results_mto1_vs_wt.csv` – example DE results table (mto1 vs WT), loaded by default.
  - `legacy_scripts/` – original analysis scripts kept as reference and loaded by the app.
  - `description_files/` – annotation files (Methylome.At, TAIR10 transposable elements).
- `install.bat` – Windows launcher that installs missing packages on first run and then opens the app.
- `RA_RNAseq_analysis_app.bat` – fast launcher for subsequent runs (no install step).

## How to run on Windows

1. Install R for Windows (≥ 4.2 recommended).
2. Extract the folder, ensuring the `.bat` files are in the main directory and the app files are inside the `app_140526` folder.
3. **First time:** double-click `install.bat`. This installs all required packages and launches the app. May take several minutes the first time.
4. **Subsequent runs:** double-click `RA_RNAseq_analysis_app.bat` (skips the install step, starts faster).

If `Rscript.exe was not found in PATH`, add your R `bin` folder to the Windows PATH. It is usually:

```text
C:\Program Files\R\R-4.x.x\bin
```

## Data input

### Option 1: Load a DE results CSV / TSV

Upload any CSV or TSV file. The app will auto-detect columns by name. If the standard names (`gene_id`, `log2FoldChange`, `padj`) are not found, the app automatically treats the **first three columns** as gene ID, log2 fold change, and adjusted p-value respectively.

**Recommended columns:**

| Column | Description |
|--------|-------------|
| `gene_id` | TAIR gene identifier |
| `log2FoldChange` | log2 fold change |
| `padj` | adjusted p-value (BH/FDR) |
| `pValue` | raw p-value (optional) |
| `baseMean` | mean normalized counts – required for MA plot |
| `GO.biological.process` | GO BP terms – used by GO offspring/stress tabs |

### Option 2: Run DESeq2 from RSEM `.genes.results` files

Point the app at a folder containing files ending in `*.genes.results`. The app:

1. Scans the folder and creates an editable `colData` table.
2. You can edit the `condition` column directly in the app, or upload a separate colData CSV.
3. Select treatment vs. control, then click **Run DESeq2**.

colData file formats supported:
- New: `sample_id`, `condition`, optional `sample_label`
- Legacy: `x`, `sample`, `exp`

`sample_label` values are used as labels in the PCA plot.

## Tabs and features

### Data
- Summary of loaded genes (up/down/NS counts, presence of baseMean/PCA).
- Editable colData table (RSEM mode).
- Preview of the DE table (up to 5000 rows).
- Download DE table and normalized counts.

### DE plots
- **Volcano plot** – all genes colored by DE class (up/down/NS).
- **MA plot** – requires `baseMean` column.
- **PCA plot** – available after running DESeq2 from RSEM files.
  - Default shows only the treatment and control conditions used for the comparison.
  - Add/remove conditions using the "PCA conditions" selector in the sidebar.
  - Toggle sample name labels on/off.
  - Labels use `sample_label` column if available; otherwise falls back to sample ID.

### GO analysis
- **Enrichment** – topGO enrichment (weight01, classic, or elim algorithms; Fisher or KS statistic).
- **REVIGO-like reduction** – local semantic reduction using `rrvgo` (no internet upload needed).
- **GO offspring** – summarizes custom parent GO term sets with offspring genes.
- **Abiotic stress** – Fisher's exact test against curated stress GO parent terms.

### KEGG analysis
- Wilcoxon rank-sum enrichment against Arabidopsis KEGG pathways.
- Bubble plot and enrichment table.
- **Pathview** – generates pathway maps colored by log2FC (requires internet access).

### TE analysis
- **Enrichment** – Wilcoxon rank-sum test to identify overrepresented TE superfamilies in DE genes.
- **Volcano plot** – volcano showing genes overlapping selected TE superfamilies.

### Gene groups
- Built-in curated gene sets (RdDM pathway, histone methyltransferases, seed-specific genes, metabolic pathways, and more).
- Volcano plot shows only the genes in the selected group, colored by DE status (same color scheme as main volcano).
- Reference lists are downloaded from GitHub on first click; subsequent groups load from cache.

### Search Annotations
- Full-text search across all annotation columns in the loaded DE table.
- Supports AND (space) and OR (|) logic. Column-level filters available in the table.

### Log / Help
- **Log** tab: timestamped record of all actions performed in the session (data loaded, analyses run, errors).
- **Help** tab: usage notes and AI developer instructions for future feature additions and debugging.

## Sidebar controls

### 1. Data input
Choose between the example CSV, upload a custom CSV/TSV, or run DESeq2 from RSEM files.

### 2. Thresholds
- **padj cutoff** – significance threshold for DE classification.
- **|log2FC| cutoff** – fold-change threshold.
- **GO p-value cutoff** – filters the GO bubble plot and table (without re-running).
- **GO ontology** – BP, MF, or CC.
- **Top GO terms** – number of terms shown in the bubble plot.

### 3. Plot Settings
- Per-tab width/height sliders for each plot.
- **PCA options** (shown when on the DE plots tab): label toggle, condition multi-select.
- **Point size** and **Point alpha** (global).
- **Point colors** – click the color swatches to open a visual color picker for Upregulated, Downregulated, and Not significant points. Applied to volcano, MA, and gene group plots.

## Notes

- Designed primarily for Arabidopsis TAIR IDs.
- GO analysis uses `topGO`, `org.At.tair.db`, and `GO.db`.
- REVIGO-like analysis uses the `rrvgo` package locally (no data is uploaded externally).
- PCA requires the count matrix and is only available after running DESeq2 from RSEM files.
- MA plot requires `baseMean`; the example CSV does not include it.
- No results are saved automatically — use the download buttons on each tab.

## Troubleshooting: R not found by the launcher

The launcher searches for `Rscript.exe` in PATH and in common Windows install locations:

- `C:\Program Files\R\...`
- `C:\Program Files (x86)\R\...`
- `%LOCALAPPDATA%\Programs\R\...`

If it still cannot find R, run `diagnose_R_installation.bat`, then open `install_RA_RNAseq_analysis_app.bat` in Notepad and set the path manually near the top:

```bat
set "RSCRIPT=C:\Program Files\R\R-4.5.0\bin\Rscript.exe"
```

Save and run again.
