# RNA-Seq Differential Expression Pipeline

An end-to-end bioinformatics pipeline written in R for processing transcript-level expression data and identifying key biological pathways.

##  Key Features
* **Data Import:** Streamlined parsing of Kallisto abundance outputs using `tximport` and `EnsDb.Hsapiens.v86`.
* **Quality Control & Normalization:** Implements TMM normalization and filtering thresholds via `edgeR`.
* **Exploratory Data Analysis:** Generates interactive interactive PCA plots and hierarchical clustering dendrograms.
* **Differential Expression:** Utilizes `limma-voom` linear modeling for robust statistical testing.
* **Functional Annotation:** Performs Gene Ontology (GO) enrichment analysis to uncover biological insights.
