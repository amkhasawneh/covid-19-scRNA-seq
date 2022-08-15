################################Loading#########################################
library(Seurat)
library(presto)
library(msigdbr)
library(fgsea)
library(tidyverse)

library(escape) #devtools::install_github("ncborcherding/escape@dev")
library(dittoSeq) #BiocManager::install("dittoSeq")
library(SingleCellExperiment)

library(SeuratObject)

library(enrichR)



library(DOSE)
library(clusterProfiler)

BCR <- readRDS("05-BCR-combined.rds")

################################Gene Info - Wilcoxon############################

BCR.genes <- wilcoxauc(BCR, 'azimuthNames')

BCR.genes
dplyr::count(BCR.genes, group)

msigdbr_species()
m_df<- msigdbr(species = "Homo sapiens", category = "H")
m_df

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

fgsea_sets




















