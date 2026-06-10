## Overview

This repository contains an R-based analytical pipeline developed for the analysis of mouse brain cortex synaptosome datasets.

The source code is publicly available to ensure the reproducibility of our analyses and to facilitate independent validation of the findings reported in the manuscript. The repository includes scripts required to reproduce most of the analyses and figures described in the manuscript.

Please note that certain auxiliary scripts, such as those used solely for generating expression pattern plots of selected individual molecules, are not included in this repository.

To execute the pipeline successfully, users must install several required R packages, which are listed in the **Requirements** section below.

## Feedback and Support

We welcome feedback from the research community. If you have any questions, encounter issues, or would like to provide comments or suggestions, please feel free to contact us.\

## Workflow

### 1. Generate identification and quantification lists

#### The identification and quantification lists used as input files for the downstream analyses are available from Zenodo.

**Glycopeptide data**

Input files:

-   `whole_gpep.xlsx`: whole-brain homogenate glycopeptide data

-   `sps_gpep.xlsx`: synaptosome fraction glycopeptide data

These files were generated from glycopeptide identification results obtained by analyzing the raw LC-MS data using Glyco-Decipher.

**IGOT data**

Input file:

-   `igot_v3_restrict.xlsx`: whole-brain homogenate and synaptosome fraction IGOT data

This file was generated from isotope-coded glycosylation-site-specific tagging (IGOT) analysis of glycopeptide samples, followed by analysis of the resulting raw data using MASCOT.

### 2. Data preprocessing and normalization

#### Purpose

-   Data selection

-   Normalization

-   Missing value imputation

-   Glycopeptide/IGOT calculation (IGOT-normalized glycopeptides)

#### Scripts

-   pretreat_v2.Rmd

-   utils_pretreat_v2.R

#### Output

-   sps_gpep_igot_ratio_4.1.xlsx

-   whole_gpep_igot_ratio_4.1.xlsx

-   sps_all_igot_imp.xlsx

-   whole_all_igot_imp.xlsx

### 3. Temporal pattern analysis of glycopeptides and IGOT-normalized glycopeptides

#### Purpose

-   Identification of glycopeptides showing significant temporal changes during brain development

-   Classification of temporal expression patterns

-   Correlation analysis between glycopeptide abundance and IGOT abundance

#### Input

-   sps_gpep_igot_ratio_4.1.xlsx

-   whole_gpep_igot_ratio_4.1.xlsx

#### Scripts

-   sps_pattern_analysis.4.0.Rmd

-   utils_pattern_analysis.4.0.R

-   sps_correlation_analysis.4.0.Rmd

-   utils_correlation_analysis.R

#### output

-   sps_pattern_result.4.1.xlsx

-   whole_pattern_result.4.1.xlsx

### 4. Temporal pattern analysis of IGOT

#### Purpose

-   IGOT data pattern analysis

-   IGOT data GO enrichment analysis

-   Comparison of glycosylation sites identified by IGOT and glycopeptide analyses

#### Input

-   sps_all_igot_imp.xlsx

-   whole_all_igot_imp.xlsx

-   sps_pattern_result.4.1.xlsx

-   whole_pattern_result.4.1.xlsx

#### Scripts

-   all_igot_analysis.Rmd

-   utils_all_igot_analysis.R

### 5. Glycan remodeling analysis

#### Input

-   sps_pattern_result.4.1.xlsx

-   whole_pattern_result.4.1.xlsx

### a)Abundance pattern Heatmap analysis

#### purpose

-   Categorize Glycan-driven protein and Protein-driven protein

-   Visualize abundance pattern of each Glycopeptide by Heatmap

#### scripts

-   sps_heatmaply_modified_final.Rmd

-   utils_heatmap_modified.R

### b)Glycan remodeling protein analysis

#### purpose

-   Analysis of four representative glycan-core peptides

-   Vizualization of abundance patterns of target Glycopeptide

#### scripts

-   CFB_line_plot.Rmd

-   utils_CFB_line_plot.R

### c)Differential branch HexNAc analysis

#### purpose

-   Identification of proteins exhibiting differential HexNAc modification during brain development

-   Functional characterization of proteins exhibiting decreased Branch HexNAc glycosylation

#### scripts

-   HexNAc_freq_analysis.Rmd

-   utils_HeNAc_freq_analysis.R

-   Domain_analysis.Rmd

-   utils_Domain_analysis.R
