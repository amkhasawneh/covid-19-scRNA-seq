library(monocle3)
library(SeuratWrappers)

ex <- covid@assays$RNA@counts
meta <- covid@meta.data
anot <- rownames(covid)

gc()
cds <- as.cell_data_set(readRDS("04-covid-clustered.rds"))


cds <- cluster_cells(cds = cds, reduction_method = "UMAP")

cds <- learn_graph(cds, use_partition = TRUE)

cds <- order_cells(cds,reduction_method = "UMAP")

cds <- reduce_dimension(cds, reduction_method="tSNE")

plot_cells(cds,color_cells_by = "azimuthNames", show_trajectory_graph = F)


