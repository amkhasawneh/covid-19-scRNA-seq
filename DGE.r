################################Loading#########################################
library(escape) #devtools::install_github("ncborcherding/escape@dev")
library(dittoSeq) #BiocManager::install("dittoSeq")
library(SingleCellExperiment)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(enrichR)
library(presto)
library(msigdbr)
library(fgsea)
library(tidyverse)
library(DOSE)
library(clusterProfiler)

BCR <- readRDS("05-BCR-combined.rds")
BCR <- BCR[,BCR$patient != "Control 4"]
BCR@meta.data$sample <- droplevels(BCR@meta.data$sample)

colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
              "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
              "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))
################################Gene Set Enrichment Analysis####################
#Importing gene sets:
gene.sets <- getGeneSets(library = "H")

#Performing the enrichment on the RNA count data:
ES <- enrichIt(obj = BCR, 
               gene.sets = gene.sets, 
               groups = 1000, cores = 8)

BCR <- AddMetaData(BCR, ES)

dittoHeatmap(BCR, genes = NULL, metas = names(ES), 
             annot.by = "sample", 
             fontsize = 7, 
             cluster_cols = F,
             heatmap.colors = rev(colorblind_vector(50)))

dittoHeatmap(BCR, genes = NULL, 
             metas = c("HALLMARK_APOPTOSIS", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_P53_PATHWAY", 
                       "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
                       "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING", 
                       "HALLMARK_MYC_TARGETS_V1", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"), #, "HALLMARK_TNFA_SIGNALING_VIA_NFKB" 
             annot.by = "sample", 
             fontsize = 7,
             heatmap.colors = rev(colorblind_vector(50)))

dittoPlot(BCR, "HALLMARK_DNA_REPAIR", group.by = "sample") + 
    scale_fill_manual(values = colorblind_vector(17))

################################Marker genes####################################

#Printing out the most abundant V-J combinations for each sample:
hiclo <- NULL
for (i in levels(factor(BCR@meta.data$sample))) {
  BCR@meta.data[BCR$sample == i,] %>%
    group_by(v_gene, j_gene, sample)  %>% dplyr::count() %>% na.omit() %>%
    arrange(desc(n)) %>% as.data.frame() -> matrix
    hiclo <- rbind(hiclo, matrix[1,])
  
}



BCR@meta.data %>% group_by(v_gene, j_gene) %>% na.omit() %>% count() %>% 
  arrange(desc(n)) %>% head(40) 


BCR <- readRDS("05-BCR-combined.rds")

#For Patient 1:
pt1 <- BCR[,BCR$patient == "Patient 1" & BCR$v_gene == "IGHV1-18"]
pt1$sample <- droplevels(pt1$sample)
DefaultAssay(pt1) <- "RNA"
Idents(pt1)
# perform a fast Wilcoxon rank sum test with Seurat:

mrk.pt1 <- FindAllMarkers(pt1, logfc.threshold = 0.58, min.pct = 0.25, only.pos = F) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)

head(mrk.pt1)

# we have all the genes for each cluster
dplyr::count(mrk.pt1, cluster)


write.table(mrk.pt1, file = "top-markers-pt1.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- mrk.pt1 %>%
  mutate(rank = -log10(p_val_adj) * sign(avg_log2FC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank
#head(n = 10)

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-seurat-pt1.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUpmod <- ranked.genes %>% 
  dplyr::filter(cluster == "moderate272_Patient1") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDownmod <- ranked.genes %>%
  dplyr::filter(cluster == "moderate272_Patient1") %>%
  filter(rank < 0) %>%
  top_n(50, wt=-p_val_adj)

topUpcrit <- ranked.genes %>% 
  dplyr::filter(cluster == "critical293_Patient1") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDowncrit <- ranked.genes %>% 
  dplyr::filter(cluster == "critical293_Patient1") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-p_val_adj)


topPathways <- bind_rows(topUpmod, topDownmod, topUpcrit, topDowncrit) #%>% arrange(-rank)

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-seurat-pt1.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(cluster) %>% subset(p_val_adj < 0.05)

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = pt1, 
                     features = as.character(unique(top$gene), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = "graphs/seurat-rrho-pt1.jpeg", dpi = "print",
       width = 15, height = 10,
       plot = DoHeatmap(alldata, 
                        features = unique(top$gene), group.by = "sample",
                        size = 2, raster=TRUE, label=FALSE) + 
         theme(axis.text.y = element_text(size = 7)))

# Perform enrichment on top 10 genes enriched in Patient 1:
enrich.topUpcrit1 <- enrichr(genes = topUpcrit$gene, 
                             databases = "GO_Biological_Process_2021")
enrich.topUpmod1 <- enrichr(genes = topUpmod$gene, 
                            databases = "GO_Biological_Process_2021")


ggsave(filename = "graphs/enrichment-pt1-seurat-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod1[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 1 moderate (IGHV1-18)"))

ggsave(filename = "graphs/enrichment-pt1-seurat-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit1[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 1 critical (IGHV1-18)"))


rm(list=setdiff(ls(), c("BCR", "hiclo")))
#For Patient 2:
pt2 <- BCR[,BCR$patient == "Patient 2" & BCR$v_gene == "IGHV4-34"]
pt2$sample <- droplevels(pt2$sample)
DefaultAssay(pt2) <- "RNA"

# perform a fast Wilcoxon rank sum test with Seurat

mrk.pt2 <- FindAllMarkers(pt2, logfc.threshold = 0.58, min.pct = 0.25, only.pos = F) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj) 

head(mrk.pt2)

# we have all the genes for each cluster
dplyr::count(mrk.pt2, cluster)


write.table(mrk.pt2, file = "top-seurat-markers-pt2.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- mrk.pt2 %>%
  mutate(rank = -log10(p_val_adj) * sign(avg_log2FC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-seurat-pt2.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUpmod <- ranked.genes %>% 
  dplyr::filter(cluster == "moderate303_Patient2") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDownmod <- ranked.genes %>%
  dplyr::filter(cluster == "moderate303_Patient2") %>%
  filter(rank < 0) %>%
  top_n(50, wt=-p_val_adj)

topUpcrit <- ranked.genes %>% 
  dplyr::filter(cluster == "critical308_Patient2") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDowncrit <- ranked.genes %>% 
  dplyr::filter(cluster == "critical308_Patient2") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-p_val_adj)


topPathways <- bind_rows(topUpmod, topDownmod, topUpcrit, topDowncrit) #%>% arrange(-rank)

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-seurat-pt2.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(cluster)

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = pt2, 
                     features = as.character(unique(top$gene), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = "graphs/seurat-rrho-pt2.jpeg", dpi = "print",
       width = 15, height = 10,
       plot = DoHeatmap(alldata, 
                        features = unique(top$gene), group.by = "sample",
                        size = 2, raster=TRUE, label=FALSE) + 
         theme(axis.text.y = element_text(size = 7)))

# Perform enrichment on top genes enriched in Patient 2
enrich.topUpcrit2 <- enrichr(genes = topUpcrit$gene, 
                             databases = "GO_Biological_Process_2021")#none significant
enrich.topUpmod2 <- enrichr(genes = topUpmod$gene, 
                            databases = "GO_Biological_Process_2021")


ggsave(filename = "graphs/enrichment-pt2-seurat-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit2[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 2 critical (IGHV4-34)"))

ggsave(filename = "graphs/enrichment-pt2-seurat-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod2[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 2 moderate (IGHV4-34)"))


rm(list=setdiff(ls(), "BCR"))

#For Patient 3:
pt3 <- BCR[,BCR$patient == "Patient 3" & BCR$v_gene == "IGHV3-48"]
pt3$sample <- droplevels(pt3$sample)
DefaultAssay(pt3) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto

mrk.pt3 <- FindAllMarkers(pt3, logfc.threshold = 0.58, min.pct = 0.25, only.pos = F) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)

head(mrk.pt3)

# we have all the genes for each cluster
dplyr::count(mrk.pt3, cluster)


write.table(mrk.pt3, file = "top-seurat-markers-v-3-48-pt3.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- mrk.pt3 %>%
  mutate(rank = -log10(p_val_adj) * sign(avg_log2FC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-seurat-v-3-48-pt3.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUpmod <- ranked.genes %>% 
  dplyr::filter(cluster == "mild186_Patient3") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDownmod <- ranked.genes %>%
  dplyr::filter(cluster == "mild186_Patient3") %>%
  filter(rank < 0) %>%
  top_n(50, wt=-p_val_adj)

topUpcrit <- ranked.genes %>% 
  dplyr::filter(cluster == "critical213_Patient3") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDowncrit <- ranked.genes %>% 
  dplyr::filter(cluster == "critical213_Patient3") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-p_val_adj)


topPathways <- bind_rows(topUpmod, topDownmod, topUpcrit, topDowncrit) #%>% arrange(-rank)

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-seurat-v-3-48-pt3.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(cluster) 

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = pt3, 
                     features = as.character(unique(top$gene), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = "graphs/seurat-rrho-v-3-48-pt3.jpeg", dpi = "print",
       width = 15, height = 10,
       plot = DoHeatmap(alldata, 
                        features = unique(top$gene), group.by = "sample",
                        size = 2, raster=TRUE, label=FALSE) + 
         theme(axis.text.y = element_text(size = 7)))

# Perform enrichment on top 10 genes enriched in Patient 3:
enrich.topUpcrit3 <- enrichr(genes = topUpcrit$gene, 
                             databases = "GO_Biological_Process_2021")
enrich.topUpmod3 <- enrichr(genes = topUpmod$gene, 
                            databases = "GO_Biological_Process_2021")


ggsave(filename = "graphs/enrichment-pt3-seurat-crit-v-3-48.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit3[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 3 critical (IGHV3-48)"))

ggsave(filename = "graphs/enrichment-pt3-seurat-mod-v-3-48.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod3[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 3 moderate (IGHV3-48)"))


rm(list=setdiff(ls(), c("BCR", "hiclo")))

#For Patient 4:
pt4 <- BCR[,BCR$patient == "Patient 4" & BCR$v_gene == "IGHV3-23"]
pt4$sample <- droplevels(pt4$sample)
DefaultAssay(pt4) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto

mrk.pt4 <- FindAllMarkers(pt4, logfc.threshold = 0.58, min.pct = 0.25, only.pos = F) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)

head(mrk.pt4)

# we have all the genes for each cluster
dplyr::count(mrk.pt4, cluster)


write.table(mrk.pt4, file = "top-wilcox-markers-v-3-23-pt4.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- mrk.pt4 %>%
  mutate(rank = -log10(p_val_adj) * sign(avg_log2FC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank


# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-v-3-23-pt4.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUpmod <- ranked.genes %>% 
  dplyr::filter(cluster == "mild227_Patient4") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDownmod <- ranked.genes %>%
  dplyr::filter(cluster == "mild227_Patient4") %>%
  filter(rank < 0) %>%
  top_n(50, wt=-p_val_adj)

topUpcrit <- ranked.genes %>% 
  dplyr::filter(cluster == "critical238_Patient4") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDowncrit <- ranked.genes %>% 
  dplyr::filter(cluster == "critical238_Patient4") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-p_val_adj)


topPathways <- bind_rows(topUpmod, topDownmod, topUpcrit, topDowncrit) #%>% arrange(-rank)

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-v-3-23-pt4.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(cluster) 

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = pt4, 
                     features = as.character(unique(top$gene), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = "graphs/seurat-rrho-v-3-23-pt4.jpeg", dpi = "print",
       width = 15, height = 10,
       plot = DoHeatmap(alldata, 
                        features = unique(top$gene), group.by = "sample",
                        size = 2, raster=TRUE, label=FALSE) + 
         theme(axis.text.y = element_text(size = 7)))

# Perform enrichment on top 10 genes enriched in Patient 4:
enrich.topUpcrit4 <- enrichr(genes = topUpcrit$gene, 
                             databases = "GO_Biological_Process_2021")
enrich.topUpmod4 <- enrichr(genes = topUpmod$gene, 
                            databases = "GO_Biological_Process_2021")


ggsave(filename = "graphs/enrichment-pt4-seurat-v-3-23-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit4[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 4 critical (IGHV3-23)"))

ggsave(filename = "graphs/enrichment-pt4-seurat-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod4[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 4 moderate (IGHV3-23)"))

rm(list=setdiff(ls(), c("BCR", "hiclo")))
#For Patient 5:
pt5 <- BCR[,BCR$patient == "Patient 5" & BCR$v_gene == "IGHV3-23"]
pt5$sample <- droplevels(pt5$sample)
DefaultAssay(pt5) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto

mrk.pt5 <- FindAllMarkers(pt5, logfc.threshold = 0.58, min.pct = 0.25, only.pos = F) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)

head(mrk.pt5)

# we have all the genes for each cluster
dplyr::count(mrk.pt5, cluster)


write.table(mrk.pt5, file = "top-seurat-markers-pt5.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- mrk.pt5 %>%
  mutate(rank = -log10(p_val_adj) * sign(avg_log2FC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-seurat-pt5.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUpcrit <- ranked.genes %>% 
  dplyr::filter(cluster == "critical119_Patient5") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDowncrit <- ranked.genes %>% 
  dplyr::filter(cluster == "critical119_Patient5") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-p_val_adj)

topUpsev <- ranked.genes %>% 
  dplyr::filter(cluster == "severe123_Patient5") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDownsev <- ranked.genes %>% 
  dplyr::filter(cluster == "severe123_Patient5") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-p_val_adj)

topUpmod <- ranked.genes %>% 
  dplyr::filter(cluster == "moderate138_Patient5") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDownmod <- ranked.genes %>%
  dplyr::filter(cluster == "moderate138_Patient5") %>%
  filter(rank < 0) %>%
  top_n(50, wt=-p_val_adj)

topPathways <- bind_rows(topUpcrit, topDowncrit, topUpsev, topDownsev, topUpmod, topDownmod)

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-seurat-pt5.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(cluster) 

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = pt5, 
                     features = as.character(unique(top$gene), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = "graphs/seurat-rrho-pt5.jpeg", dpi = "print",
       width = 15, height = 10,
       plot = DoHeatmap(alldata, 
                        features = unique(top$gene), group.by = "sample",
                        size = 4, raster=TRUE, label=FALSE) + 
         theme(axis.text.y = element_text(size = 4)))

# Perform enrichment on top 10 genes enriched in Patient 5:
enrich.topUpcrit5 <- enrichr(genes = topUpcrit$gene, 
                             databases = "GO_Biological_Process_2021")
enrich.topUpsev5 <- enrichr(genes = topUpsev$gene, 
                            databases = "GO_Biological_Process_2021")
enrich.topUpmod5 <- enrichr(genes = topUpmod$gene, 
                            databases = "GO_Biological_Process_2021")


ggsave(filename = "graphs/enrichment-pt5-seurat-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit5[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 5 critical (IGHV3-23)"))

ggsave(filename = "graphs/enrichment-pt5-seurat-sev.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpsev5[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 5 severe (IGHV3-23)"))

ggsave(filename = "graphs/enrichment-pt5-seurat-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod5[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 5 moderate (IGHV3-23)"))


rm(list=setdiff(ls(), c("BCR", "hiclo")))
#For Patient 6:
pt6 <- BCR[,BCR$patient == "Patient 6" & BCR$v_gene == "IGHV3-23"]
pt6$sample <- droplevels(pt6$sample)
DefaultAssay(pt6) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto

mrk.pt6 <- FindAllMarkers(pt6, logfc.threshold = 0.58, min.pct = 0.25, only.pos = F) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)

head(mrk.pt6)

# we have all the genes for each cluster
dplyr::count(mrk.pt6, cluster)


write.table(top_markers(wlx.mrk.pt6), file = "top-seurat-markers-pt6.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- mrk.pt6 %>%
  mutate(rank = -log10(p_val_adj) * sign(avg_log2FC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-v-3-23-pt6.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUpcrit <- ranked.genes %>% 
  dplyr::filter(cluster == "critical120_Patient6") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDowncrit <- ranked.genes %>% 
  dplyr::filter(cluster == "critical120_Patient6") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-p_val_adj)

topUpsev <- ranked.genes %>% 
  dplyr::filter(cluster == "severe122_Patient6") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDownsev <- ranked.genes %>% 
  dplyr::filter(cluster == "severe122_Patient6") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-p_val_adj)

topUpmod <- ranked.genes %>% 
  dplyr::filter(cluster == "moderate124_Patient6") %>%
  filter(rank > 0) %>% 
  top_n(60, wt=-p_val_adj)

topDownmod <- ranked.genes %>%
  dplyr::filter(cluster == "moderate124_Patient6") %>%
  filter(rank < 0) %>%
  top_n(50, wt=-p_val_adj)

topPathways <- bind_rows(topUpcrit, topDowncrit, topUpsev, topDownsev, topUpmod, topDownmod)

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-v-3-23-pt6.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(cluster) 

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = pt6, 
                     features = as.character(unique(top$gene), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = "graphs/seurat-rrho-pt6.jpeg", dpi = "print",
       width = 15, height = 10,
       plot = DoHeatmap(alldata, 
                        features = unique(top$gene), group.by = "sample",
                        size = 4, raster=TRUE, label=FALSE) + 
         theme(axis.text.y = element_text(size = 4)))

# Perform enrichment on top 10 genes enriched in Patient 6:
enrich.topUpcrit6 <- enrichr(genes = topUpcrit$gene, 
                             databases = "GO_Biological_Process_2021")
enrich.topUpsev6 <- enrichr(genes = topUpsev$gene, 
                            databases = "GO_Biological_Process_2021")
enrich.topUpmod6 <- enrichr(genes = topUpmod$gene, 
                            databases = "GO_Biological_Process_2021")


ggsave(filename = "graphs/enrichment-pt6-seurat-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit6[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 6 critical (IGHV3-23)"))

ggsave(filename = "graphs/enrichment-pt6-seurat-sev.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpsev6[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 6 severe (IGHV3-23)"))

ggsave(filename = "graphs/enrichment-pt6-seurat-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod6[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 6 moderate (IGHV3-23)"))


rm(list=setdiff(ls(),c("BCR","hiclo")))
gc()

#Patient 1 vs. all:
DefaultAssay(BCR) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto
BCR$comp <- 0
BCR$comp[BCR$patient == "Patient 1"] <- 1

wlx.mrk <- wilcoxauc(BCR, "comp", c(1,0),
                     seurat_assay='RNA', assay = "data") %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)

head(wlx.mrk)

# we have all the genes for each cluster
dplyr::count(wlx.mrk, group)

# summarize the top abundant marker genes for each group
top_markers(wlx.mrk)

write.table(top_markers(wlx.mrk), file = "top-wilcox-markers-pt1-vs-all.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk %>%
  mutate(rank = -log10(p_val_adj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank
#head(n = 10)

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-pt1-vs-all.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUppt1 <- ranked.genes %>% 
  dplyr::filter(group == 1) %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDownpt1 <- ranked.genes %>%
  dplyr::filter(group == 1) %>%
  filter(rank < 0) %>%
  top_n(50, wt=-p_val_adj)

topUpother <- ranked.genes %>% 
  dplyr::filter(group == 0) %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDownother <- ranked.genes %>% 
  dplyr::filter(group == 0) %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-p_val_adj)


topPathways <- bind_rows(topUppt1, topDownpt1, topUpother, topDownother) #%>% arrange(-rank)

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-pt1-vs-all.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(group) 

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = BCR, 
                     features = as.character(unique(top$gene), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = "graphs/wilcox-rrho-pt1-vs-all.jpeg", dpi = "print",
       width = 15, height = 10,
       plot = DoHeatmap(alldata, 
                        features = unique(top$gene), group.by = "patient",
                        size = 2, raster=TRUE, label=FALSE) + 
         theme(axis.text.y = element_text(size = 7)))

# Perform enrichment on top genes enriched in Patient 1:
enrich.topUppt1 <- enrichr(genes = topUppt1$gene, 
                           databases = "GO_Biological_Process_2021")
enrich.topUpother <- enrichr(genes = topUpother$gene, 
                             databases = "GO_Biological_Process_2021")


ggsave(filename = "graphs/enrichment-pt1-vs-all-wlx.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUppt1[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 1 vs Others"))

ggsave(filename = "graphs/enrichment-all-vs-pt1-wlx.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpother[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Others vs Patient 1"))


#IGHV1-18 vs. all:
DefaultAssay(BCR) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto
BCR$comp <- 0
BCR$comp[BCR$patient == "Patient 1"] <- 1

wlx.mrk <- wilcoxauc(BCR, "comp", c(1,0),
                     seurat_assay='RNA', assay = "data") %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)

head(wlx.mrk)

# we have all the genes for each cluster
dplyr::count(wlx.mrk, group)

# summarize the top abundant marker genes for each group
top_markers(wlx.mrk)

write.table(top_markers(wlx.mrk), file = "top-wilcox-markers-pt1-vs-all.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk %>%
  mutate(rank = -log10(p_val_adj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank
#head(n = 10)

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-pt1-vs-all.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUppt1 <- ranked.genes %>% 
  dplyr::filter(group == 1) %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDownpt1 <- ranked.genes %>%
  dplyr::filter(group == 1) %>%
  filter(rank < 0) %>%
  top_n(50, wt=-p_val_adj)

topUpother <- ranked.genes %>% 
  dplyr::filter(group == 0) %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDownother <- ranked.genes %>% 
  dplyr::filter(group == 0) %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-p_val_adj)


topPathways <- bind_rows(topUppt1, topDownpt1, topUpother, topDownother) #%>% arrange(-rank)

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-pt1-vs-all.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(group) 

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = BCR, 
                     features = as.character(unique(top$gene), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = "graphs/wilcox-rrho-pt1-vs-all.jpeg", dpi = "print",
       width = 15, height = 10,
       plot = DoHeatmap(alldata, 
                        features = unique(top$gene), group.by = "patient",
                        size = 2, raster=TRUE, label=FALSE) + 
         theme(axis.text.y = element_text(size = 7)))

# Perform enrichment on top genes enriched in Patient 1:
enrich.topUppt1 <- enrichr(genes = topUppt1$gene, 
                           databases = "GO_Biological_Process_2021")
enrich.topUpother <- enrichr(genes = topUpother$gene, 
                             databases = "GO_Biological_Process_2021")


ggsave(filename = "graphs/enrichment-pt1-vs-all-wlx.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUppt1[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 1 vs Others"))

ggsave(filename = "graphs/enrichment-all-vs-pt1-wlx.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpother[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Others vs Patient 1"))




#For IGHV1-18 in Patient 1:
pt1 <- BCR[,BCR$patient == "Patient 1" & BCR$v_gene == "IGHV1-18"]
pt1$sample <- droplevels(pt1$sample)
DefaultAssay(pt1) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto

wlx.mrk.pt1 <- wilcoxauc(pt1, 'sample', c("moderate272_Patient1", "critical293_Patient1"),
                         seurat_assay='RNA', assay = "data") %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)

head(wlx.mrk.pt1)

# we have all the genes for each cluster
dplyr::count(wlx.mrk.pt1, group)

# summarize the top abundant marker genes for each group
top_markers(wlx.mrk.pt1)

write.table(top_markers(wlx.mrk.pt1), file = "top-wilcox-markers-v-1-18-pt1.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk.pt1 %>%
  mutate(rank = -log10(p_val_adj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank
#head(n = 10)

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-v-1-18-pt1.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUpmod <- ranked.genes %>% 
  dplyr::filter(group == "moderate272_Patient1") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDownmod <- ranked.genes %>%
  dplyr::filter(group == "moderate272_Patient1") %>%
  filter(rank < 0) %>%
  top_n(50, wt=-p_val_adj)

topUpcrit <- ranked.genes %>% 
  dplyr::filter(group == "critical293_Patient1") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-p_val_adj)

topDowncrit <- ranked.genes %>% 
  dplyr::filter(group == "critical293_Patient1") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-p_val_adj)


topPathways <- bind_rows(topUpmod, topDownmod, topUpcrit, topDowncrit) #%>% arrange(-rank)

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-v-1-18-pt1.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(group) %>% subset(p_val_adj < 0.05)

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = pt1, 
                     features = as.character(unique(top$gene), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = "graphs/wilcox-rrho-v-1-18-pt1.jpeg", dpi = "print",
       width = 15, height = 10,
       plot = DoHeatmap(alldata, 
                        features = unique(top$gene), group.by = "sample",
                        size = 2, raster=TRUE, label=FALSE) + 
         theme(axis.text.y = element_text(size = 7)))

# Perform enrichment on top 10 genes enriched in Patient 1:
enrich.topUpcrit1 <- enrichr(genes = topUpcrit$gene, 
                             databases = "GO_Biological_Process_2021")
enrich.topUpmod1 <- enrichr(genes = topUpmod$gene, 
                            databases = "GO_Biological_Process_2021")


ggsave(filename = "graphs/enrichment-v-1-18-pt1-wlx-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod1[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 1 moderate (IGHV1-18)"))

ggsave(filename = "graphs/enrichment-v-1-18-pt1-wlx-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit1[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 1 critical (IGHV1-18)"))



################################Wilcoxon Test Method############################

BCR <- readRDS("05-BCR-combined.rds")

#For Patient 1:
pt1 <- BCR[,BCR$patient == "Patient 1" & BCR$v_gene == "IGHV1-18"]
pt1$sample <- droplevels(pt1$sample)
DefaultAssay(pt1) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto

wlx.mrk.pt1 <- wilcoxauc(pt1, 'sample', c("moderate272_Patient1", "critical293_Patient1"),
                               seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(padj)

head(wlx.mrk.pt1)

# we have all the genes for each cluster
dplyr::count(wlx.mrk.pt1, group)

# summarize the top abundant marker genes for each group
top_markers(wlx.mrk.pt1)

write.table(top_markers(wlx.mrk.pt1), file = "top-wilcox-markers-pt1.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk.pt1 %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank
#head(n = 10)

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-pt1.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUpmod <- ranked.genes %>% 
  dplyr::filter(group == "moderate272_Patient1") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDownmod <- ranked.genes %>%
  dplyr::filter(group == "moderate272_Patient1") %>%
  filter(rank < 0) %>%
  top_n(50, wt=-padj)

topUpcrit <- ranked.genes %>% 
  dplyr::filter(group == "critical293_Patient1") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDowncrit <- ranked.genes %>% 
  dplyr::filter(group == "critical293_Patient1") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-padj)


topPathways <- bind_rows(topUpmod, topDownmod, topUpcrit, topDowncrit) #%>% arrange(-rank)

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-pt1.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(group) %>% subset(padj < 0.05)

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = pt1, 
                     features = as.character(unique(top$feature), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = "graphs/wilcox-rrho-pt1.jpeg", dpi = "print",
       width = 15, height = 10,
       plot = DoHeatmap(alldata, 
          features = unique(top$feature), group.by = "sample",
          size = 2, raster=TRUE, label=FALSE) + 
  theme(axis.text.y = element_text(size = 7)))

# Perform enrichment on top 10 genes enriched in Patient 1:
enrich.topUpcrit1 <- enrichr(genes = topUpcrit$feature, 
                               databases = "GO_Biological_Process_2021")
enrich.topUpmod1 <- enrichr(genes = topUpmod$feature, 
                               databases = "GO_Biological_Process_2021")


ggsave(filename = "graphs/enrichment-pt1-wlx-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod1[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 1 moderate"))

ggsave(filename = "graphs/enrichment-pt1-wlx-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit1[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 1 critical"))

#For Patient 2:
pt2 <- BCR[,BCR$patient == "Patient 2"]
pt2$sample <- droplevels(pt2$sample)
DefaultAssay(pt2) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto

wlx.mrk.pt2 <- wilcoxauc(pt2, 'sample', c("moderate303_Patient2", "critical308_Patient2"),
                               seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(padj)

head(wlx.mrk.pt2)

# we have all the genes for each cluster
dplyr::count(wlx.mrk.pt2, group)

# summarize the top abundant marker features for each group
top_markers(wlx.mrk.pt2)

write.table(top_markers(wlx.mrk.pt2), file = "top-wilcox-markers-pt2.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk.pt2 %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank
#head(n = 10)

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-pt2.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUpmod <- ranked.genes %>% 
  dplyr::filter(group == "moderate303_Patient2") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDownmod <- ranked.genes %>%
  dplyr::filter(group == "moderate303_Patient2") %>%
  filter(rank < 0) %>%
  top_n(50, wt=-padj)

topUpcrit <- ranked.genes %>% 
  dplyr::filter(group == "critical308_Patient2") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDowncrit <- ranked.genes %>% 
  dplyr::filter(group == "critical308_Patient2") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-padj)


topPathways <- bind_rows(topUpmod, topDownmod, topUpcrit, topDowncrit) #%>% arrange(-rank)

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-pt2.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(group)

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = pt2, 
                     features = as.character(unique(top$feature), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = "graphs/wilcox-rrho-pt2.jpeg", dpi = "print",
       width = 15, height = 10,
       plot = DoHeatmap(alldata, 
          features = unique(top$feature), group.by = "sample",
          size = 2, raster=TRUE, label=FALSE) + 
  theme(axis.text.y = element_text(size = 7)))

# Perform enrichment on top genes enriched in Patient 2
enrich.topUpcrit2 <- enrichr(genes = topUpcrit$feature, 
                               databases = "GO_Biological_Process_2021")
enrich.topUpmod2 <- enrichr(genes = topUpmod$feature, 
                               databases = "GO_Biological_Process_2021")


ggsave(filename = "graphs/enrichment-pt2-wlx-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit2[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 2 critical"))

ggsave(filename = "graphs/enrichment-pt2-wlx-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod2[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 2 moderate"))


rm(list=setdiff(ls(), "BCR"))
#For Patient 3:
pt3 <- BCR[,BCR$patient == "Patient 3"]
pt3$sample <- droplevels(pt3$sample)
DefaultAssay(pt3) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto

wlx.mrk.pt3 <- wilcoxauc(pt3, 'sample', c("mild186_Patient3", "critical213_Patient3"),
                               seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(padj)

head(wlx.mrk.pt3)

# we have all the genes for each cluster
dplyr::count(wlx.mrk.pt3, group)

# summarize the top abundant marker features for each group
top_markers(wlx.mrk.pt3)

write.table(top_markers(wlx.mrk.pt3), file = "top-wilcox-markers-pt3.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk.pt3 %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank
#head(n = 10)

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-pt3.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUpmod <- ranked.genes %>% 
  dplyr::filter(group == "mild186_Patient3") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDownmod <- ranked.genes %>%
  dplyr::filter(group == "mild186_Patient3") %>%
  filter(rank < 0) %>%
  top_n(50, wt=-padj)

topUpcrit <- ranked.genes %>% 
  dplyr::filter(group == "critical213_Patient3") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDowncrit <- ranked.genes %>% 
  dplyr::filter(group == "critical213_Patient3") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-padj)


topPathways <- bind_rows(topUpmod, topDownmod, topUpcrit, topDowncrit) #%>% arrange(-rank)

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-pt3.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(group) 

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = pt3, 
                     features = as.character(unique(top$feature), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = "graphs/wilcox-rrho-pt3.jpeg", dpi = "print",
       width = 15, height = 10,
       plot = DoHeatmap(alldata, 
          features = unique(top$feature), group.by = "sample",
          size = 2, raster=TRUE, label=FALSE) + 
  theme(axis.text.y = element_text(size = 7)))

# Perform enrichment on top 10 genes enriched in Patient 3:
enrich.topUpcrit3 <- enrichr(genes = topUpcrit$feature, 
                               databases = "GO_Biological_Process_2021")
enrich.topUpmod3 <- enrichr(genes = topUpmod$feature, 
                               databases = "GO_Biological_Process_2021")


ggsave(filename = "graphs/enrichment-pt3-wlx-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit3[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 3 critical"))

ggsave(filename = "graphs/enrichment-pt3-wlx-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod3[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 3 moderate"))


rm(list=setdiff(ls(), "BCR"))
#For Patient 4:
pt4 <- BCR[,BCR$patient == "Patient 4"]
pt4$sample <- droplevels(pt4$sample)
DefaultAssay(pt4) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto

wlx.mrk.pt4 <- wilcoxauc(pt4, 'sample', c("mild227_Patient4", "critical238_Patient4"),
                               seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(padj)

head(wlx.mrk.pt4)

# we have all the genes for each cluster
dplyr::count(wlx.mrk.pt4, group)

# summarize the top abundant marker features for each group
top_markers(wlx.mrk.pt4)

write.table(top_markers(wlx.mrk.pt4), file = "top-wilcox-markers-pt4.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk.pt4 %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank
#head(n = 10)

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-pt4.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUpmod <- ranked.genes %>% 
  dplyr::filter(group == "mild227_Patient4") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDownmod <- ranked.genes %>%
  dplyr::filter(group == "mild227_Patient4") %>%
  filter(rank < 0) %>%
  top_n(50, wt=-padj)

topUpcrit <- ranked.genes %>% 
  dplyr::filter(group == "critical238_Patient4") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDowncrit <- ranked.genes %>% 
  dplyr::filter(group == "critical238_Patient4") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-padj)


topPathways <- bind_rows(topUpmod, topDownmod, topUpcrit, topDowncrit) #%>% arrange(-rank)

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-pt4.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(group) 

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = pt4, 
                     features = as.character(unique(top$feature), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = "graphs/wilcox-rrho-pt4.jpeg", dpi = "print",
       width = 15, height = 10,
       plot = DoHeatmap(alldata, 
          features = unique(top$feature), group.by = "sample",
          size = 2, raster=TRUE, label=FALSE) + 
  theme(axis.text.y = element_text(size = 7)))

# Perform enrichment on top 10 genes enriched in Patient 4:
enrich.topUpcrit4 <- enrichr(genes = topUpcrit$feature, 
                               databases = "GO_Biological_Process_2021")
enrich.topUpmod4 <- enrichr(genes = topUpmod$feature, 
                               databases = "GO_Biological_Process_2021")


ggsave(filename = "graphs/enrichment-pt4-wlx-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit4[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 4 critical"))

ggsave(filename = "graphs/enrichment-pt4-wlx-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod4[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 4 moderate"))

#For Patient 5:
pt5 <- BCR[,BCR$patient == "Patient 5"]
pt5$sample <- droplevels(pt5$sample)
DefaultAssay(pt5) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto

wlx.mrk.pt5 <- wilcoxauc(pt5, 'sample',
                         c("critical119_Patient5", "severe123_Patient5", "moderate138_Patient5"),
                               seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(padj)

head(wlx.mrk.pt5)

# we have all the genes for each cluster
dplyr::count(wlx.mrk.pt5, group)

# summarize the top abundant marker features for each group
top_markers(wlx.mrk.pt5)

write.table(top_markers(wlx.mrk.pt5), file = "top-wilcox-markers-pt5.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk.pt5 %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank
#head(n = 10)

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-pt5.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUpcrit <- ranked.genes %>% 
  dplyr::filter(group == "critical119_Patient5") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDowncrit <- ranked.genes %>% 
  dplyr::filter(group == "critical119_Patient5") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-padj)

topUpsev <- ranked.genes %>% 
  dplyr::filter(group == "severe123_Patient5") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDownsev <- ranked.genes %>% 
  dplyr::filter(group == "severe123_Patient5") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-padj)

topUpmod <- ranked.genes %>% 
  dplyr::filter(group == "moderate138_Patient5") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDownmod <- ranked.genes %>%
  dplyr::filter(group == "moderate138_Patient5") %>%
  filter(rank < 0) %>%
  top_n(50, wt=-padj)

topPathways <- bind_rows(topUpcrit, topDowncrit, topUpsev, topDownsev, topUpmod, topDownmod)

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-pt5.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(group) 

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = pt5, 
                     features = as.character(unique(top$feature), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = "graphs/wilcox-rrho-pt5.jpeg", dpi = "print",
       width = 15, height = 10,
       plot = DoHeatmap(alldata, 
          features = unique(top$feature), group.by = "sample",
          size = 4, raster=TRUE, label=FALSE) + 
  theme(axis.text.y = element_text(size = 4)))

# Perform enrichment on top 10 genes enriched in Patient 5:
enrich.topUpcrit5 <- enrichr(genes = topUpcrit$feature, 
                               databases = "GO_Biological_Process_2021")
enrich.topUpsev5 <- enrichr(genes = topUpsev$feature, 
                               databases = "GO_Biological_Process_2021")
enrich.topUpmod5 <- enrichr(genes = topUpmod$feature, 
                               databases = "GO_Biological_Process_2021")


ggsave(filename = "graphs/enrichment-pt5-wlx-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit5[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 5 critical"))

ggsave(filename = "graphs/enrichment-pt5-wlx-sev.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpsev5[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 5 severe"))

ggsave(filename = "graphs/enrichment-pt5-wlx-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod5[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 5 moderate"))

#For Patient 6:
pt6 <- BCR[,BCR$patient == "Patient 6"]
pt6$sample <- droplevels(pt6$sample)
DefaultAssay(pt6) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto

wlx.mrk.pt6 <- wilcoxauc(pt6, 'sample',
                         c("critical120_Patient6", "severe122_Patient6", "moderate124_Patient6"),
                         seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(padj)

head(wlx.mrk.pt6)

# we have all the genes for each cluster
dplyr::count(wlx.mrk.pt6, group)

# summarize the top abundant marker features for each group
top_markers(wlx.mrk.pt6)

write.table(top_markers(wlx.mrk.pt6), file = "top-wilcox-markers-pt6.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk.pt6 %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank
#head(n = 10)

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-pt6.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUpcrit <- ranked.genes %>% 
  dplyr::filter(group == "critical120_Patient6") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDowncrit <- ranked.genes %>% 
  dplyr::filter(group == "critical120_Patient6") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-padj)

topUpsev <- ranked.genes %>% 
  dplyr::filter(group == "severe122_Patient6") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDownsev <- ranked.genes %>% 
  dplyr::filter(group == "severe122_Patient6") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-padj)

topUpmod <- ranked.genes %>% 
  dplyr::filter(group == "moderate124_Patient6") %>%
  filter(rank > 0) %>% 
  top_n(60, wt=-padj)

topDownmod <- ranked.genes %>%
  dplyr::filter(group == "moderate124_Patient6") %>%
  filter(rank < 0) %>%
  top_n(50, wt=-padj)

topPathways <- bind_rows(topUpcrit, topDowncrit, topUpsev, topDownsev, topUpmod, topDownmod)

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-pt6.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(group) 

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = pt6, 
                     features = as.character(unique(top$feature), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = "graphs/wilcox-rrho-pt6.jpeg", dpi = "print",
       width = 15, height = 10,
       plot = DoHeatmap(alldata, 
                        features = unique(top$feature), group.by = "sample",
                        size = 4, raster=TRUE, label=FALSE) + 
         theme(axis.text.y = element_text(size = 4)))

# Perform enrichment on top 10 genes enriched in Patient 6:
enrich.topUpcrit6 <- enrichr(genes = topUpcrit$feature, 
                             databases = "GO_Biological_Process_2021")
enrich.topUpsev6 <- enrichr(genes = topUpsev$feature, 
                            databases = "GO_Biological_Process_2021")
enrich.topUpmod6 <- enrichr(genes = topUpmod$feature, 
                            databases = "GO_Biological_Process_2021")


ggsave(filename = "graphs/enrichment-pt6-wlx-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit6[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 6 critical"))

ggsave(filename = "graphs/enrichment-pt6-wlx-sev.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpsev6[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 6 severe"))

ggsave(filename = "graphs/enrichment-pt6-wlx-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod6[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 6 moderate"))


rm(list=setdiff(ls(), "BCR"))
gc()

#Patient 1 vs. all:
DefaultAssay(BCR) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto
BCR$comp <- 0
BCR$comp[BCR$patient == "Patient 1"] <- 1

wlx.mrk <- wilcoxauc(BCR, "comp", c(1,0),
                               seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(padj)

head(wlx.mrk)

# we have all the genes for each cluster
dplyr::count(wlx.mrk, group)

# summarize the top abundant marker features for each group
top_markers(wlx.mrk)

write.table(top_markers(wlx.mrk), file = "top-wilcox-markers-pt1-vs-all.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank
#head(n = 10)

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-pt1-vs-all.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUppt1 <- ranked.genes %>% 
  dplyr::filter(group == 1) %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDownpt1 <- ranked.genes %>%
  dplyr::filter(group == 1) %>%
  filter(rank < 0) %>%
  top_n(50, wt=-padj)

topUpother <- ranked.genes %>% 
  dplyr::filter(group == 0) %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDownother <- ranked.genes %>% 
  dplyr::filter(group == 0) %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-padj)


topPathways <- bind_rows(topUppt1, topDownpt1, topUpother, topDownother) #%>% arrange(-rank)

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-pt1-vs-all.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(group) 

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = BCR, 
                     features = as.character(unique(top$feature), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = "graphs/wilcox-rrho-pt1-vs-all.jpeg", dpi = "print",
       width = 15, height = 10,
       plot = DoHeatmap(alldata, 
          features = unique(top$feature), group.by = "patient",
          size = 2, raster=TRUE, label=FALSE) + 
  theme(axis.text.y = element_text(size = 7)))

# Perform enrichment on top genes enriched in Patient 1:
enrich.topUppt1 <- enrichr(genes = topUppt1$feature, 
                               databases = "GO_Biological_Process_2021")
enrich.topUpother <- enrichr(genes = topUpother$feature, 
                               databases = "GO_Biological_Process_2021")


ggsave(filename = "graphs/enrichment-pt1-vs-all-wlx.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUppt1[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 1 vs Others"))

ggsave(filename = "graphs/enrichment-all-vs-pt1-wlx.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpother[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Others vs Patient 1"))


#IGHV1-18 vs. all:
DefaultAssay(BCR) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto
BCR$comp <- 0
BCR$comp[BCR$patient == "Patient 1"] <- 1

wlx.mrk <- wilcoxauc(BCR, "comp", c(1,0),
                               seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(padj)

head(wlx.mrk)

# we have all the genes for each cluster
dplyr::count(wlx.mrk, group)

# summarize the top abundant marker features for each group
top_markers(wlx.mrk)

write.table(top_markers(wlx.mrk), file = "top-wilcox-markers-pt1-vs-all.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank
#head(n = 10)

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-pt1-vs-all.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUppt1 <- ranked.genes %>% 
  dplyr::filter(group == 1) %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDownpt1 <- ranked.genes %>%
  dplyr::filter(group == 1) %>%
  filter(rank < 0) %>%
  top_n(50, wt=-padj)

topUpother <- ranked.genes %>% 
  dplyr::filter(group == 0) %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDownother <- ranked.genes %>% 
  dplyr::filter(group == 0) %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-padj)


topPathways <- bind_rows(topUppt1, topDownpt1, topUpother, topDownother) #%>% arrange(-rank)

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-pt1-vs-all.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(group) 

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = BCR, 
                     features = as.character(unique(top$feature), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = "graphs/wilcox-rrho-pt1-vs-all.jpeg", dpi = "print",
       width = 15, height = 10,
       plot = DoHeatmap(alldata, 
          features = unique(top$feature), group.by = "patient",
          size = 2, raster=TRUE, label=FALSE) + 
  theme(axis.text.y = element_text(size = 7)))

# Perform enrichment on top genes enriched in Patient 1:
enrich.topUppt1 <- enrichr(genes = topUppt1$feature, 
                               databases = "GO_Biological_Process_2021")
enrich.topUpother <- enrichr(genes = topUpother$feature, 
                               databases = "GO_Biological_Process_2021")


ggsave(filename = "graphs/enrichment-pt1-vs-all-wlx.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUppt1[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 1 vs Others"))

ggsave(filename = "graphs/enrichment-all-vs-pt1-wlx.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpother[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Others vs Patient 1"))






################################Gene Ontology Analysis##########################

#Checking available databases to perform enrichment (then choosing one (or more)):
dbs <- listEnrichrDbs()


#Enrichment on marker genes:
enriched <- enrichr(databases = c("GO_Molecular_Function_2021",
                                  "GO_Cellular_Component_2021",
                                  "GO_Biological_Process_2021"),
                    genes = mrk$gene)


# barplot

par(mfrow = c(1, 1), mar = c(3, 25, 2, 1))
barplot(height = -log10(enriched[[1]]$Adjusted.P.value[10:1]), names.arg = enriched[[1]]$Term[10:1],
        horiz = TRUE, las = 1, border = FALSE, cex.names = 0.55)
abline(v = c(-log10(0.05)), lty = 2)
abline(v = 0, lty = 1)


