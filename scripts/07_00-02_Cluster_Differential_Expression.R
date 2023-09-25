# 07_Cluster_Differential_Expression
# Jacob Mitchel
# 08/10/2022

# Differential expression of each cluster vs the other cells cannot run locally
# due to memory constraints and will be conducted on Rockfish

library(Seurat)
library(parallel)
library(dplyr)

sessionInfo()

results_dir <- "processed_data/07_Cluster_Differential_Expression"
if(!dir.exists(results_dir)){
  dir.create(results_dir)
}
figure_dir <- "figures/07_Cluster_Differential_Expression"
if(!dir.exists(figure_dir)){
  dir.create(figure_dir)
}

writeLines(capture.output(sessionInfo()), paste0(results_dir, "/sessionInfo.txt"))

seurat <- readRDS("processed_data/05_UMAP_Embedding_and_Clustering/05_UMAP_Cluster_mouse_ici_scRNAseq.rds")

target_clusters <- levels(seurat$Louvain_cluster_0.6)[1:3]
mcFindMarkers <- function(i){
  ident1 <- i
  df <- FindMarkers(seurat, group.by = "Louvain_cluster_0.6",
                    ident.1 = ident1,
                    test.use = "MAST",
                    min.pct = -Inf,
                    logfc.threshold = -Inf,
                    min.cells.feature = 1,
                    min.cells.group = 1)
  df$Gene.name.unique <- rownames(df)
  df$cluster <- rep(i, nrow(df))
  write.csv(df, file = paste0(results_dir, "/Louvain_0.6_cluster_",
                              i, "_MAST_de_results.csv"))
  return(df)
}

marker_results <- list()[target_clusters]
ptm <- proc.time()
marker_results <- parallel::mclapply(target_clusters, mcFindMarkers, mc.cores = 16)
time_diff <- proc.time() - ptm
print(time_diff)
