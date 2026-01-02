# GEO Data Ingestion Pipeline

This R script automates the process of searching, downloading, and organizing transcriptomic datasets from NCBI GEO (Gene Expression Omnibus). It's designed to support early-stage in silico analyses by providing a clean, standardized input layer.

We use it internally at **PharmaToxAI** to streamline access to public gene expression data relevant to toxicology and drug development, but it's general enough to adapt to other domains.

## What it does

* Searches GEO Series (GSE) based on a custom keyword (e.g., *"Semaglutide"*, *"cosmetics ingredients"*) and organism (e.g., *Homo sapiens*, *Mus musculus*).
* Downloads SeriesMatrix files and supplementary raw data (when available).
* Organizes everything into a structured folder tree for downstream processing.
* Generates:

  * A metadata summary (CSV) of all hits.
  * A download log.
  * A quick scan of expression stats (min/median/max, guess on scale, candidate grouping columns, etc.).

* **NEW: Advanced Quality Analysis & Visualizations:**
  * Comprehensive quality control reports for each dataset
  * Multiple visualization types (distribution plots, PCA, correlation heatmaps, boxplots)
  * Statistical calculations (CV, correlations, outlier detection, variance analysis)
  * Automated detection of problematic samples and low-quality data

No differential expression analysis is done at this stage. This is focused on ingestion, organization, and quality assessment.

## Features

### Data Ingestion
- Automated GEO database search
- Batch download of SeriesMatrix files
- Download of supplementary raw data files
- Structured organization of downloaded datasets

### Quality Control & Analysis
- **Expression Distribution Analysis**: Histograms and density plots of expression values
- **Sample Correlation Analysis**: Heatmaps showing inter-sample correlations
- **Principal Component Analysis (PCA)**: Dimensionality reduction and visualization
- **Sample Boxplots**: Distribution of expression values across samples
- **Coefficient of Variation (CV)**: Analysis of gene-level variability
- **Outlier Detection**: Identification of problematic samples based on correlation
- **Missing Data Assessment**: Statistics on data completeness
- **Low-Variance Gene Detection**: Identification of uninformative features

### Output Files
- `geo_search_results.csv`: Metadata for all found datasets
- `geo_download_log.csv`: Download status for each dataset
- `geo_eset_scan_summary.csv`: Quick expression statistics
- `quality_analysis_summary.csv`: Consolidated quality metrics for all datasets
- `quality_reports/[GSE_ID]/`: Individual quality reports and visualizations
  - `*_quality_report.txt`: Text-based quality summary
  - `*_expression_distribution.pdf`: Expression value distributions
  - `*_sample_correlation.pdf`: Sample correlation heatmap
  - `*_pca.pdf`: PCA plots with variance explained
  - `*_sample_boxplots.pdf`: Per-sample expression boxplots
  - `*_cv_distribution.pdf`: Coefficient of variation distribution

## Dependencies

Tested on R ≥ 4.2 with the following packages:

* `rentrez`
* `GEOquery`
* `dplyr`
* `purrr`
* `readr`
* `stringr`
* `Biobase` (optional, used in the summary step)

## Example use case

```r
SEARCH_TERM <- "cosmetics ingredients"
ORGANISM    <- "Homo sapiens"
```

This query fetched several GSEs related to human response to common cosmetic agents. Matrix and raw files were retrieved, structured, and summarized for quick inspection before differential expression analysis.

## Why it matters

Getting reliable, properly structured public data is a recurring bottleneck in many bioinformatics workflows. This tool helps remove friction and enforces a basic level of reproducibility from the start.

The new quality control features help identify potential issues early:
- **Outlier samples** that might skew downstream analyses
- **Low-quality datasets** with excessive missing data
- **Batch effects** visible in PCA plots
- **Scale and normalization issues** detected through distribution analysis
- **Low-variance genes** that provide little information

This early QC step can save significant time by highlighting datasets that require additional preprocessing or should be excluded from further analysis.
