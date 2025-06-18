# Introduction:
This document describes the RNA-seq data analysis performed for the STARProtocol (2025) manuscript titled **Protocol to Investigate Metabolically Altered Pathways  in Small Cell Lung Cancer: From RNA-Seq Analysis to Seahorse-Based Functional Validation** published in journal **"STAR Protocols"** in 2025.

# Operating System:
R version 4.3.1 (2023-06-16 ucrt) on Microsoft Windows 11 Education Version 10.0.26100 Build 26100

# Installation of R-packages
Install R through CRAN and then install packages listed above through CRAN or Bioconductor. Typical install time is ~1 hour.

To install Bioconductor packages, do:
```
install.packages("BiocManager")
BiocManager::install("packageName")
```
To install R packages through CRAN:
```
install.packages("packageName")
```

# Installation of bioinformatics tools
Install bioinformatics tools on Linux based operating system.  

Providing guidance to install bioinformatics tools is beyond the scope of this manuscript.


# Code availability:
Please follow [this link](https://sagarutturkar.github.io/STARProtocol_TranLab_2025/) to view the document describing the Bioinformatics commands for RNAseq data analysis (raw FASTQ to counts matrix) and R-code (from counts matrix to Differentially Expressed Genes [DEGs]) and generate the figures described in the manuscript.

Complete code as R Markdown file (index.Rmd) is available via [this GitHub Repository](https://github.com/sagarutturkar/STARProtocol_TranLab_2025/blob/main/docs/index.Rmd).
