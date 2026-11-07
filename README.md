# Bulk RNA-seq differential expression analysis in R

An R analysis script that imports Kallisto transcript abundance estimates, summarizes expression using a human transcript annotation, and compares **disease versus healthy** samples with edgeR and limma-voom.

This repository documents a transcriptomics learning project. It contains analysis code, not a validated biological finding or a complete raw-read processing workflow. The study metadata and Kallisto files are not included, so the original analysis cannot currently be reproduced from this repository alone.

## Workflow

1. Read sample metadata and check the expected input files.
2. Import Kallisto `abundance.tsv` files with `tximport` and `lengthScaledTPM` counts.
3. Map transcripts to gene symbols using `EnsDb.Hsapiens.v86`.
4. Retain genes with CPM > 1 in at least five samples and apply TMM normalization.
5. Explore expression distributions, hierarchical clustering, and PCA.
6. Fit a limma-voom model with the contrast `disease - healthy`.
7. Export differential expression tables and explore GO Biological Process enrichment.

## Repository contents

| File | Purpose |
| --- | --- |
| `my first rna  _seq analysis.R` | Analysis script |
| `studydesign.example.tsv` | Illustrative metadata format; not real study data |
| `README.md` | Setup, methods, outputs, and limitations |

## Requirements

Install R, then install the packages below. Package versions are not yet locked; the script exports `sessionInfo.txt` to record the environment of a completed run.

```r
install.packages(c("BiocManager", "tidyverse", "matrixStats", "cowplot",
                   "plotly", "gt", "pheatmap", "RColorBrewer"))
BiocManager::install(c("tximport", "EnsDb.Hsapiens.v86", "ensembldb",
                      "edgeR", "limma", "statmod", "clusterProfiler", "org.Hs.eg.db"))
```

## Input layout

Run from the repository directory. Supply a tab-separated `studydesign.txt` with columns `sample` and `group`. Each sample identifier must be unique and match its folder name; group labels must be exactly `healthy` and `disease`, with biological replication in each group.

```text
studydesign.txt
mappedReads/
  sample_01/abundance.tsv
  sample_02/abundance.tsv
  ...
```

Use `studydesign.example.tsv` only as a format reference. It does not supply expression data. The current filtering rule requires expression in at least five samples; review that rule before applying this script to another study.

Kallisto quantification is an upstream prerequisite. Confirm that the transcriptome used for quantification is compatible with the Ensembl release 86 human annotation. The script ignores transcript version suffixes during matching.

## Run

Open R or RStudio with the repository as the working directory:

```r
source("my first rna  _seq analysis.R")
```

The original script was written for interactive exploration. In a sourced or noninteractive session, explicitly print ggplot objects to display them. Plots are not automatically exported as image files.

## Outputs and interpretation

| Output | Description |
| --- | --- |
| `DGE_all_genes_results.tsv` | Tested genes with limma statistics and BH-adjusted p-values |
| `up_regulated_genes.csv` | Adjusted p-value < 0.01 and log2 fold change >= 2 |
| `down_regulated_genes.csv` | Adjusted p-value < 0.01 and log2 fold change <= -2 |
| `sessionInfo.txt` | R and package versions for the completed run |

Positive model log2 fold changes indicate higher expression in disease. The exploratory mean log-CPM difference is separate from the model estimate and must not replace it. Heatmap and enrichment selection use `decideTests` with adjusted p-value 0.01 and an absolute log2 fold-change threshold of 2.

GO enrichment uses gene symbols and the filtered, tested gene set as its background. Enrichment is exploratory; pathway significance does not establish a causal mechanism.

## Current limitations

- Dataset accession, sample provenance, original sample metadata, and quantification inputs still need to be supplied before a reproducible biological report can be produced.
- The model includes group only; pairing, batches, and other covariates need study-specific assessment.
- Gene-symbol mapping can be ambiguous. Stable gene identifiers and annotation coverage should be reviewed before research use.
- A fixed expression filter may not suit a different sample size or design.
- Very small significant gene sets, no significant genes, or empty enrichment results can require additional plotting guards.
- No executed analysis report, benchmark, or independently verified result is provided here.

## Next research milestones

Add the verified dataset accession and experimental design, document data retrieval, lock the software environment, run the full analysis, and publish a concise report with QC figures, effect sizes, uncertainty, and biological limitations. Any study or tutorial used as a source should be acknowledged when its provenance is confirmed.
