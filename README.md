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

No analysis is done at this stage. This is just ingestion and organization.

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
