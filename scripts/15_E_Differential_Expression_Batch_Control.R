# 15_E_Differential_Expression_Batch_Control
# Jacob Mitchell
# 05/28/23

# differential expression is run for the full feature set in order to use the results for
# gene set enrichment analysis

library(Seurat)
library(dplyr)

set.seed(123)
sessionInfo()

# create directories
result_dir <- "processed_data/15_E_Differential_Expression_Batch_control"
if(!dir.exists(result_dir)){
  dir.create(result_dir)
}
figure_dir <- "figures/15_E_Differential_Expression_Batch_control"
if(!dir.exists(figure_dir)){
  dir.create(figure_dir)
}

# read seurat object
pancvax <- readRDS(file = "processed_data/15_A_PancVAX2_Prelim_Analysis/seurat_untreated_PancVAX_PancVAXICI.rds")

mast.test <- function(seurat, group.by, ident.1, ident.2, filename){
  de <- FindMarkers(object = seurat,
                    test.use = "MAST",
                    group.by = group.by,
                    ident.1 = ident.1, ident.2 = ident.2,
                    min.pct = -Inf,
                    logfc.threshold = -Inf,
                    min.cells.feature = 1,
                    min.cells.group = 1,
                    latent.vars = c("batch"))
  de$Gene.name.unique <- rownames(de)
  write.csv(x = de, file = filename)
  return(de)
}

cell_groupings <- list(
  "CD8_cells" = c("cycling_CD8_T", "cytotoxic_CD8_T", "exhausted_CD8_T", "naive_CD8_T"),
  "CD8_nonExhausted" = c("cycling_CD8_T", "cytotoxic_CD8_T", "naive_CD8_T"),
  "CD8_effector" = c("cycling_CD8_T", "cytotoxic_CD8_T"),
  "CD8_cytotoxic" = c("cytotoxic_CD8_T"),
  "CD8_cycling" = c("cycling_CD8_T"),
  "CD8_exhausted" = c("exhausted_CD8_T"),
  "CD8_naive" = c("naive_CD8_T"),
  "CD4_cells" = c("effector_CD4_T", "regulatory_CD4_T"),
  "CD4_effector" = c("effector_CD4_T"),
  "CD4_Treg" = c("regulatory_CD4_T"),
  "Myeloid" = c("dendritic_cell", "monocyte", "myeloid_derived_supressor", "tumor_associated_macrophage"),
  "Monocyte+TAM" = c("monocyte", "tumor_associated_macrophage"),
  "TAM" = c("tumor_associated_macrophage"),
  "Monocyte" = c("monocyte"),
  "Dendritic_cell" = c("dendritic_cell"),
  "MDSC" = c("myeloid_derived_supressor"),
  "B_cell" = c("B"),
  "NK_cell" = c("natural_killer")
)

cell_de_list <- list()

for(cell in names(cell_groupings)){
  print(cell)
  ser <- pancvax[, pancvax$cell_type %in% cell_groupings[[cell]] ]
  
  # test untreated vs PancVAX
  testname_PvU <- paste0(cell, "_PvU")
  de_PvU <- mast.test(seurat = ser, group.by = "treatment",
                      ident.1 = "PancVAX", ident.2 = "Untreated", 
                      filename = paste0(result_dir, "/MAST_test_BatchCorrection_", testname_PvU, ".csv"))
  cell_de_list[[testname_PvU]] <- de_PvU
  
  # test PancVAX vs PancVAX + ICI
  testname_IvP <- paste0(cell, "_IvP")
  de_IvP <- mast.test(seurat = ser, group.by = "treatment",
                      ident.1 = "PancVAX + Anti-PD-1 + Anti-CTLA-4", ident.2 = "PancVAX", 
                      filename = paste0(result_dir, "/MAST_test_BatchCorrection_", testname_IvP, ".csv"))
  cell_de_list[[testname_IvP]] <- de_IvP
}

# clear memory
rm(list = ls())
