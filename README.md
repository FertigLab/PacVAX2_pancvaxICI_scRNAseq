# MHC class II-specific neoantigens enhance anti-tumor immunity of a personalized vaccine in murine pancreatic cancer

## Abstract

Personalized cancer vaccines aim to activate and expand cytotoxic anti-tumor CD8$^{+}$ T cells to recognize and kill tumor cells. However, the role of CD4$^{+}$ T cell activation in the clinical benefit of these vaccines is not well defined. We previously established a personalized neoantigen vaccine (PancVAX) for the pancreatic cancer cell line Panc02, which activates tumor-specific CD8$^{+}$ T cells but required combinatorial checkpoint modulators to achieve therapeutic efficacy. To determine the effects of neoantigen-specific CD4+ T cell activation, we generated a new vaccine (PancVAX2) targeting both MHCI- and MHCII-specific neoantigens. Tumor-bearing mice vaccinated with PancVAX2 had significantly improved control of tumor growth and long-term survival benefit without concurrent administration of checkpoint inhibitors. PancVAX2 significantly enhanced priming and recruitment of neoantigen-specific CD8+ T into the tumor with lower PD1 expression after reactivation compared to the CD8$^{+}$ vaccine alone. Vaccine-induced neoantigen-specific Th1 CD4$^{+}$ T cells in the tumor were associated with decreased T regulatory cells (Tregs). Consistent with this, PancVAX2 was associated with more pro-immune myeloid-derived suppressor cells and M1-like macrophages in the tumor demonstrating a less immunosuppressive tumor microenvironment. This study demonstrates the biological importance of prioritizing and including CD4 T cell-specific neoantigens for personalized cancer vaccine modalities.

## PacVAX2_pancvaxICI_scRNAseq
scRNAseq and survival analysis in Panc02-bearing mice treated with PancVAX and immune checkpoint inhibition.
 
This repository contains the scripts necessary to replicate the results in Figures 1 and S1 of Huff et al. 

## Script Order

### Starting Raw Data

Raw data for single-cell RNA-seq and TCR-seq can be downloaded from -----. Sequencing alignment was conducted with 10x Genomics CellRanger v4.0.0 using the mm10 reference genome. The resulting counts matrix should be stored in 

```
data/filtered_bc_matrix/barcodes.tsv.gz
data/filtered_bc_matrix/features.tsv.gz
data/filtered_bc_matrix/matrix.mtx.gz
```

Aligned TCR sequences should be stored in

```
data/filtered_contigs/1T1_filtered_contig_annotations.csv
data/filtered_contigs/1T2_filtered_contig_annotations.csv
data/filtered_contigs/1T3_filtered_contig_annotations.csv
data/filtered_contigs/1T4_filtered_contig_annotations.csv
data/filtered_contigs/1T5_filtered_contig_annotations.csv
data/filtered_contigs/1T6_filtered_contig_annotations.csv
data/filtered_contigs/1T7_filtered_contig_annotations.csv
data/filtered_contigs/1T8_filtered_contig_annotations.csv
```

Survival data should be stored in

```
data/20230828_complete/PancVAX + ICI survival tidy.csv
```

### Directory Creation 

Before beginning analysis, directories should be created to store processed data results, figures, and intermediate data structures. Sub-directories are named after the scripts that generate the files stored within.

#### data

Unprocessed data used in the initialization of the analysis pipeline.

#### processed_data

Intermediate data structures and resulting data tables.

#### figures

Graphical figures generated during analysis.

#### reports

Logs and reports from batch processing of scripts

### Scripts

#### 00_runCellranger_4.0.0_aggr.sh \& 00_runCellranger_4.0.0_mm10_vdj.sh

Parameters for running CellRanger v4.0.0 to assemble sequencing results to obtain the aggregated expression matrix and TCR contigs.

#### 01_read_10x_merge_with_TCR.Rmd

Load the RNA counts matrix into R as a Seurat object. Integrate TCR sequences as cell meta data using Scirpy in Python.

#### 02_seurat_preprocessing_qc.Rmd

Quality control filtering of cells and genes in the RNA counts matrix. Genes were removed if they did not have non-zero expression in at least 3 cells. Cells were maintained in the counts matrix if they expressed between 200 and 3000 unique genes, had between 1000 and 10,000 UMIs, had counts from mitochondrial genes comprising less than 10% of total RNA counts, and had no expression of hemoglobin-encoding genes. 

#### 03_cell_cycle_scoring.Rmd

Score cell cycle in cells using the module score function in Seurat and mouse orthologs of gene lists associated with S phase and G2/M phase from Kowalczyk et al., Genome Res., 2015 (10.1101/gr.192237.115).

#### 04_Normalization_Scaling_PCA_and_Harmony_Batch_Correction.Rmd

Normalize expression expression within cells, scale expression of the top 2000 variable genes across cells with regression of cell cycle S phase and G2/M phase scores, calculation of principle components on the scaled expression of the top 2000 variable genes, and batch correction of the low-dimensional PCA embedding of cells using Harmony where samples derived from whole tumors and percoll-enriched immune cells were considered separate batches.

#### 05_UMAP_Embedding_and_Clustering.Rmd

Calculate UMAP embedding and cell clusters based on Harmony dimensionality reduction. Cell neighborhoods and the UMAP embedding were calculated based on the first 25 dimensions. Cell clusters were identified using Louvain clustering at resolutions of 0.6, 1.0, 1.4, and 1.8 where higher resolutions result in separating the cells into a higher number of clusters. 

#### 06_Marker_Genes.Rmd

Generation of UMAP plots and violin plots by cluster at 0.6 resolution to assess expression of marker genes for cell types expected to be present in the Panc02 tumor microenvironment.

#### 07_##-##_Cluster_Differential_Expression.R

Identifies differentially expressed genes in each cell cluster relative to other clusters considering all gene features using a MAST test. The scripts were split for analysis of three clusters at a time to facilitate parallelization of this computationally slow step. The scripts were run on the Rockfish computing cluster operated by Advanced Research Computing at Hopkins.

- ##-##: clusters considered each script.
- 07_##-##_slurm_job.sh: Bash script for submitting each R script for batch computing using the SLURM job management system.

#### 10_Cell_Type_Labeling.Rmd

Final annotation of cell types based on Louvain clustering at resolution 0.6 and 1.0. Clusters at 0.6 resolution containing fewer than 200 cells (clusters 25, 26, 27, 28, and 29) were removed. Cells in cluster 24 were removed due to the cluster consisting of cells entirely from percoll separated cells from PancVAX + PD-1 treated mice that did not express marker genes corresponding to expected cell types, suggesting that they represented a sample processing artifact. Cells were annotated as types:

- cytotoxic_CD8_T
- exhausted_CD8_T
- cycling_CD8_T
- effector_CD4_T
- naive_CD8_T
- effector_CD4_T
- exhausted_CD8_T
- natural_killer
- regulatory_CD4_T
- B (cell)
- myeloid_derived_supressor
- monocyte
- tumor_associated_macrophage
- dendritic_cell
- endothelial
- cancer

#### 15_A_PancVAX2_Prelim_Analysis.Rmd

Comparison of functional T cell type compositions in tumors from mice treated with Isotype, Anti-PD-1 + Anti-CTLA-4, PancVAX, PancVAX + Anti-PD-1 + Anti-CTLA-4. Subsets the data set to cells from these four treatment groups. Scaled expression data is rescaled to consider the decreased number of cells. PCA is recomputed and Harmony batch-corrected embeddings are calculated. UMAP embedding is recomputed maintaining the same cell type annotations. Generates the UMAP figure of cell types colored at hatched for all cells (figure S1A). Generates the dotplot of cell type marker expression (figure S1B). The full data set is subset to T cell clusters, and the UMAP plot coloring T cells by type and treatment (figure 1C) is generated. T cell type composition as proportion of T cells belonging to each T cell functional type is compared between treatments (figure 1D). Generates UMAPs of T cells colored by expression of T cell functional type markers (figure S1C).

#### 15_E_Differential_Expression_Batch_Control.R

Differential expression analysis across treatments for T cell and myeloid cell types. Differential gene expression by MAST test is assessed within cell types comparing

- Untreated vs PancVAX
- PancVAX vs PancVAX + Anti-PD-1 + Anti-CTLA-4

Batch is included as a control variable in the test model. 

#### 15_F_Differential_Expression_Batch-Corrected_Plotting.Rmd

Generates volcano plots of differential gene expression results from the 15_E script.

#### 15_H_Figure_PDF_Copies.Rmd

Generates the figures from 15_A as pdf files.

#### 16_Complete_Survival_Analysis.Rmd

Statistical comparisons of survival probabilities between mice treated with PancVAX, Anti-PD-1 + Anti-CTLA-4, PancVAX + Anti-PD-1 + Anti-CTLA-4, and isotope treated. Survival curves were estimated by the Kaplan-Meier methods. Survival functions amongst all four groups was compared by Mantel-Haenszel test and pairwise comparisons between pairs of groups were compared by log-rank test. All tests used their implementation in R package survival v3.5-5.

#### 17_Generate_per-well_Processed_Data.Rmd

Generates rds files containing Seurat objects for cells from each sample well that underwent sequencing and analysis.

