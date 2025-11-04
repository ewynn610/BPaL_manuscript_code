# Code for *Deconvoluting drug interactions based on M. tuberculosis physiologic processes: Transcriptional disaggregation of the BPaL regimen in vivo*

**Manuscript Status**: [Now published in Antimicrobial Agents and Chemotherapy](https://journals.asm.org/doi/full/10.1128/aac.00492-25)

**Title**: *Deconvoluting drug interactions based on M. tuberculosis physiologic processes: Transcriptional disaggregation of the BPaL regimen in vivo* 

## Repository Contents

This repository includes code used to generate the analyses in the manuscript. The main scripts are described below:

- **`vst_calculation`**  
  Calculates VST-normalized expression data with batch correction using `limma::removeBatchEffect`.

- **`cfu_analysis`**  
  Scripts for visualization and statistical analysis of CFU over time across different drug treatment conditions.

- **`pca_analysis`**  
  Performs principal component analysis (PCA) and generates corresponding plots.

- **`de_analysis`**  
  Runs differential expression analyses and statistical tests comparing drug treatments and timepoints.

- **`clustering/heatmap_generation`**  
  Example code for generating heatmaps of average gene expression over time. Provided for bedaquiline (B); the same process was used for pretomanid (P), linezolid (L), and B vs. BPaL comparisons.

- **`gene_set_average_expression_analysis`**  
  Visualizes and tests for differences in average expression of predefined gene sets across treatment conditions.

Data files are not included in this repository due to size and sharing limitations. Please contact the corresponding author for access or details.

## Contact
For questions or issues, please contact Elizabeth Wynn at wynne@njhealth.org.

