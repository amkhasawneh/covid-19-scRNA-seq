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
library(DOSE)
library(clusterProfiler)
library(circlize)
library(ComplexHeatmap)

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
  BCR@meta.data[BCR$sample == i & !is.na(BCR$v_gene),] %>%
    group_by(v_gene, j_gene, sample)  %>% dplyr::count() %>% na.omit() %>%
    arrange(desc(n)) %>% as.data.frame() -> matrix
    hiclo <- rbind(hiclo, matrix[1,])
  
}



BCR@meta.data %>% group_by(v_gene, j_gene) %>% na.omit() %>% count() %>% 
  arrange(desc(n)) %>% as.data.frame() %>% head(40) 


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

rm(list=setdiff(ls(), "BCR"))

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
  top_n(50, wt=-padj) %>%
  head(10)

topDownmod <- ranked.genes %>%
  dplyr::filter(group == "moderate272_Patient1") %>%
  filter(rank < 0) %>%
  top_n(50, wt=-padj)

topUpcrit <- ranked.genes %>% 
  dplyr::filter(group == "critical293_Patient1") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj) %>%
  head(10)

topDowncrit <- ranked.genes %>% 
  dplyr::filter(group == "critical293_Patient1") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-padj)


topPathways <- bind_rows(topUpmod, topDownmod, topUpcrit, topDowncrit) #%>% arrange(-rank)

cat(topUpmod$feature, sep = ", ")
cat(topUpcrit$feature, sep = ", ")

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

cat(vapply(str_split(enrich.topUpmod1[["GO_Biological_Process_2021"]]$Term[1:20] ,"[(GO:*)]"), "[", "", 1), sep = ", ")
cat(vapply(str_split(enrich.topUpcrit1[["GO_Biological_Process_2021"]]$Term[1:20] ,"[(GO:*)]"), "[", "", 1), sep = ", ")



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


rm(list=setdiff(ls(), "BCR"))

BCR$comp <- 0
BCR$comp[BCR$v_gene == "IGHV4-34" & BCR$patient == "Patient 2"] <- 1
BCR$comp <- as.factor(BCR$comp)

#For Patient 2:
pt2 <- BCR[,BCR$patient == "Patient 2" & BCR$v_gene == "IGHV4-59"]
pt2 <- BCR[,BCR$patient == "Patient 2" & BCR$v_gene == "IGHV4-34"]
pt2$sample <- droplevels(pt2$sample)
DefaultAssay(pt2) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto

wlx.mrk.pt2 <- wilcoxauc(BCR, 'comp', c(1,0),
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
  top_n(50, wt=-padj) %>%
  head(10)

topDownmod <- ranked.genes %>%
  dplyr::filter(group == "moderate303_Patient2") %>%
  filter(rank < 0) %>%
  top_n(50, wt=-padj)

topUpcrit <- ranked.genes %>% 
  dplyr::filter(group == "critical308_Patient2") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj) %>%
  head(10)

topDowncrit <- ranked.genes %>% 
  dplyr::filter(group == "critical308_Patient2") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-padj)

cat(topUpmod$feature, sep = ", ")
cat(topUpcrit$feature, sep = ", ")

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
alldata <- ScaleData(object = BCR, 
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
enrich.topDowncrit2 <- enrichr(genes = topDowncrit$feature, 
                               databases = "GO_Biological_Process_2021")

cat(vapply(str_split(enrich.topUpmod2[["GO_Biological_Process_2021"]]$Term[1:20] ,"[(GO:*)]"), "[", "", 1), sep = ", ")
cat(vapply(str_split(enrich.topUpcrit2[["GO_Biological_Process_2021"]]$Term[1:20] ,"[(GO:*)]"), "[", "", 1), sep = ", ")
cat(vapply(str_split(enrich.topDowncrit2[["GO_Biological_Process_2021"]]$Term[1:6] ,"[(GO:*)]"), "[", "", 1), sep = ", ")
enrich.topDowncrit2[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topDowncrit2[["GO_Biological_Process_2021"]]$Term ,"[(GO:*)]"), "[", "", 1)

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

ggsave(filename = "graphs/enrichment-pt2-wlx-down-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topDowncrit2[[1]], showTerms = 6,
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Downregulated in Patient 2 critical (IGHV4-34)"))


rm(list=setdiff(ls(), "BCR"))

#For Patient 3:
pt3 <- BCR[,BCR$patient == "Patient 3" & BCR$v_gene == "IGHV4-39"]
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
  top_n(50, wt=-padj) %>%
  head(10) %>%
  head(10)

topDownmod <- ranked.genes %>%
  dplyr::filter(group == "mild186_Patient3") %>%
  filter(rank < 0) %>%
  top_n(50, wt=-padj)

topUpcrit <- ranked.genes %>% 
  dplyr::filter(group == "critical213_Patient3") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj) %>%
  head(10) %>% 
  head(10)

topDowncrit <- ranked.genes %>% 
  dplyr::filter(group == "critical213_Patient3") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-padj)

cat(topUpmod$feature, sep = ", ")
cat(topUpcrit$feature, sep = ", ")

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

cat(vapply(str_split(enrich.topUpmod3[["GO_Biological_Process_2021"]]$Term[1:20] ,"[(GO:*)]"), "[", "", 1), sep = ", ")

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
pt4 <- BCR[,BCR$patient == "Patient 4" & BCR$v_gene == "IGHV3-23"]
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
  top_n(50, wt=-padj) %>%
  head(10)

topDownmod <- ranked.genes %>%
  dplyr::filter(group == "mild227_Patient4") %>%
  filter(rank < 0) %>%
  top_n(50, wt=-padj)

topUpcrit <- ranked.genes %>% 
  dplyr::filter(group == "critical238_Patient4") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj) %>%
  head(10)

topDowncrit <- ranked.genes %>% 
  dplyr::filter(group == "critical238_Patient4") %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-padj)

cat(topUpmod$feature, sep = ", ")

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

cat(vapply(str_split(enrich.topUpmod4[["GO_Biological_Process_2021"]]$Term[1:20] ,"[(GO:*)]"), "[", "", 1), sep = ", ")


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

rm(list=setdiff(ls(), "BCR"))

#For Patient 5:
pt5 <- BCR[,BCR$patient == "Patient 5" & BCR$v_gene == "IGHV3-23"]
pt5$sample <- droplevels(pt5$sample)
DefaultAssay(pt5) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto
#, "moderate138_Patient5"

wlx.mrk.pt5 <- wilcoxauc(pt5, 'sample',
                         c("critical119_Patient5", "severe123_Patient5"),
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
cat(topUpmod$feature[1:10], sep = ", ")
cat(topUpsev$feature[1:10], sep = ", ")

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

cat(vapply(str_split(enrich.topUpmod5[["GO_Biological_Process_2021"]]$Term[1:20] ,"[(GO:*)]"), "[", "", 1), sep = ", ")
cat(vapply(str_split(enrich.topUpsev5[["GO_Biological_Process_2021"]]$Term[1:20] ,"[(GO:*)]"), "[", "", 1), sep = ", ")

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

rm(list=setdiff(ls(), "BCR"))

#For Patient 6:
pt6 <- BCR[,BCR$patient == "Patient 6" & BCR$v_gene == "IGHV4-59"]
pt6$sample <- droplevels(pt6$sample)
DefaultAssay(pt6) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto: "severe122_Patient6", 

wlx.mrk.pt6 <- wilcoxauc(pt6, 'sample',
                         c("critical120_Patient6",  "moderate124_Patient6"),
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

cat(topUpcrit$feature[1:10], sep = ", ")
cat(topUpmod$feature[1:10], sep = ", ")

# save found markers as exportable table
write.table(topPathways,
            file = "ranked-genes-rrho-pt6.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

top <- topPathways %>% group_by(group) 

cat(topUpsev$feature, sep = ", ")
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

cat(vapply(str_split(enrich.topUpmod6[["GO_Biological_Process_2021"]]$Term[1:20] ,"[(GO:*)]"), "[", "", 1), sep = ", ")
cat(vapply(str_split(enrich.topUpsev6[["GO_Biological_Process_2021"]]$Term[1:20] ,"[(GO:*)]"), "[", "", 1), sep = ", ")


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



################################ISGs Heatmaps###################################

#####ISGs in DEGs, for Patients 1-4:
#Importing interferon-stimulated gene set:
ISGs <- read.csv("manuscript/ISGs_GeneSets_SU.csv")
ISGs <- ISGs$ISGs_geneset_227
ISGs <- gsub(x = ISGs, pattern = "HLA-", replacement = "HLA.")

#Cleaning gene set:
ISGs <- ISGs[!grepl(ISGs, pattern = "LOC") & !grepl(ISGs, pattern = "XENTR") & !grepl(ISGs, pattern = "GSON")]
ISGs <- ISGs[ISGs != ""]


#Subsetting BCR data:
obj.b <- BCR[,(BCR$patient != "Patient 5" & BCR$patient != "Patient 6") & BCR$severity != "healthy"]
obj.b$sample <- droplevels(obj.b$sample)
obj.b$zeverity <- "critical"
obj.b$zeverity[obj.b$severity == "mild" | obj.b$severity == "moderate"] <- "moderate"
DefaultAssay(obj.b) <- "RNA"


#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj.b, group_by = "zeverity", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(desc(abs(logFC)))

# rank genes using rrho algorithm
MNP.genes_B <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                # sort genes by decreasing rank
#head(n = 10)

# choose top upregulated genes
topUp_B_mod <- MNP.genes_B %>% 
  dplyr::filter(group == c('moderate')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>% 
  top_n(50, wt=-padj)

topUp_B_crit <- MNP.genes_B %>% 
  dplyr::filter(group == c('critical')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>%
  top_n(50, wt=-padj)


## moderate vs critical
topPathways_all <- bind_rows(topUp_B_mod, topUp_B_crit)


#Subsetting to only DEGs that are DEGs:
obj.b <- obj.b[ISGs[ISGs %in% topPathways_all$feature],]
obj.b <- NormalizeData(obj.b)

#Changing the data to "pseudo-bulk", using ISG set:
obj.b.tmp <- ScaleData(obj.b)
avg.exp.mat.b <- AverageExpression(obj.b.tmp, group.by = "sample", slot = 'scale.data')
avg.exp.rna.b <- avg.exp.mat.b$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.b.tmp@meta.data %>% 
  group_by(sample, zeverity) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data.B <- data.frame(t(avg.exp.rna.b), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data.B <- meta.data.B[,c(colnames(meta.data.B)[colnames(meta.data.B) %in% rownames(obj.b.tmp@assays$RNA@scale.data)], "sample", "zeverity")]


#Preparing aesthetic parameters:
expr.cols <- colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow"))
col.anno <- HeatmapAnnotation(Severity = meta.data$zeverity, show_annotation_name = F,   
                              col = list(Severity = c("critical" = "red", "moderate" = "blue")))
column_labels <- c("moderate272_Patient1" = "Pt1", "critical293_Patient1" = "Pt1",
                   "moderate303_Patient2" = "Pt2", "critical308_Patient2" = "Pt2",
                   "mild186_Patient3" = "Pt3", "critical213_Patient3" = "Pt3",
                   "mild227_Patient4" = "Pt4", "critical238_Patient4" = "Pt4")

#Preparing separate matrices for each cell type:
rownames(meta.data.B) <- meta.data.B$sample
meta.data.B$zeverity <- factor(meta.data.B$zeverity, levels = c("moderate", "critical"))
hm.b <- Heatmap(t(meta.data.B[,colnames(meta.data.B) %in% diff.genes$feature]), 
                column_split = meta.data.B$zeverity, name = "Expression", 
                    column_labels = column_labels, row_title = "B cells", column_title = " ",  
                    col = expr.cols, top_annotation = col.anno, row_names_gp = gpar(fontsize = 5), 
                    show_row_dend = F, show_column_dend = F, cluster_columns = F)

#For changing cells into squares:
# width = nrow(meta.data.B)*unit(5, "mm"), 
# height = ncol(meta.data.B)*unit(5, "mm"),

#Subsetting TCR data:
#TCR <- readRDS("06-TCR-combined.rds")
obj.t <- TCR[,(TCR$patient != "Patient5" & TCR$patient != "Patient6") & TCR$severity != "healthy"]


obj.cd4 <- obj.t[,obj.t$azimuthNames == "CD4 TCM"]
obj.cd4$sample <- droplevels(obj.cd4$sample)
obj.cd4$zeverity <- "critical"
obj.cd4$zeverity[obj.cd4$severity == "mild" | obj.cd4$severity == "moderate"] <- "moderate"
DefaultAssay(obj.cd4) <- "RNA"


#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj.cd4, group_by = "zeverity", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(desc(abs(logFC)))

# rank genes using rrho algorithm
MNP.genes_CD4 <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                # sort genes by decreasing rank
#head(n = 10)

# choose top upregulated genes
topUp_CD4_mod <- MNP.genes_CD4 %>% 
  dplyr::filter(group == c('moderate')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>% 
  top_n(50, wt=-padj)

topUp_CD4_crit <- MNP.genes_CD4 %>% 
  dplyr::filter(group == c('critical')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>%
  top_n(50, wt=-padj)


## moderate vs critical
topPathways_all <- bind_rows(topUp_CD4_mod, topUp_CD4_crit)
 

#Subsetting to only DEGs that are ISGs:
obj.cd4 <- obj.cd4[ISGs[ISGs %in% topPathways_all$feature],]
obj.cd4 <- NormalizeData(obj.cd4)

#Changing the data to "pseudo-bulk", using ISG set:
obj.cd4.tmp <- ScaleData(obj.cd4)
avg.exp.mat.cd4 <- AverageExpression(obj.cd4.tmp, group.by = "sample", slot = 'scale.data')
avg.exp.rna.cd4 <- avg.exp.mat.cd4$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.cd4.tmp@meta.data %>% 
  group_by(sample, zeverity) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data.CD4 <- data.frame(t(avg.exp.rna.cd4), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data.CD4 <- meta.data.CD4[,c(colnames(meta.data.CD4)[colnames(meta.data.CD4) %in% rownames(obj.cd4.tmp@assays$RNA@scale.data)], "sample", "zeverity")]

rownames(meta.data.CD4) <- meta.data.CD4$sample
meta.data.CD4$zeverity <- factor(meta.data.CD4$zeverity, levels = c("moderate", "critical"))

hm.cd4 <- Heatmap(t(meta.data.CD4[,colnames(meta.data.CD4) %in% topPathways_all$feature]), column_split = meta.data.CD4$zeverity,
                    column_labels = column_labels, row_title = "CD4 TCM",
                    col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 5),
                    show_row_dend = F, show_column_dend = F, cluster_columns = F)



obj.cd8 <- obj.t[,obj.t$azimuthNames == "CD8 TEM"]
obj.cd8$sample <- droplevels(obj.cd8$sample)
obj.cd8$zeverity <- "critical"
obj.cd8$zeverity[obj.cd8$severity == "mild" | obj.cd8$severity == "moderate"] <- "moderate"
DefaultAssay(obj.cd8) <- "RNA"


#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj.cd8, group_by = "zeverity", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(desc(abs(logFC)))

# rank genes using rrho algorithm
MNP.genes_CD8 <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                # sort genes by decreasing rank


# choose top upregulated genes
topUp_CD8_mod <- MNP.genes_CD8 %>% 
  dplyr::filter(group == c('moderate')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>% 
  top_n(50, wt=-padj)

topUp_CD8_crit <- MNP.genes_CD8 %>% 
  dplyr::filter(group == c('critical')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>%
  top_n(50, wt=-padj)


## moderate vs critical
topPathways_all <- bind_rows(topUp_CD8_mod, topUp_CD8_crit)
 

#Subsetting to only DEGs that are ISGs:
obj.cd8 <- obj.cd8[ISGs[ISGs %in% topPathways_all$feature],]
obj.cd8 <- NormalizeData(obj.cd8)

#Changing the data to "pseudo-bulk", using ISG set:
obj.cd8.tmp <- ScaleData(obj.cd8)
avg.exp.mat.cd8 <- AverageExpression(obj.cd8.tmp, group.by = "sample", slot = 'scale.data')
avg.exp.rna.cd8 <- avg.exp.mat.cd8$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.cd8.tmp@meta.data %>% 
  group_by(sample, zeverity) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data.CD8 <- data.frame(t(avg.exp.rna.cd8), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data.CD8 <- meta.data.CD8[,c(colnames(meta.data.CD8)[colnames(meta.data.CD8) %in% rownames(obj.cd8.tmp@assays$RNA@scale.data)], "sample", "zeverity")]

rownames(meta.data.CD8) <- meta.data.CD8$sample
meta.data.CD4$zeverity <- factor(meta.data.CD4$zeverity, levels = c("moderate", "critical"))
hm.cd8 <- Heatmap(t(meta.data.CD8[,colnames(meta.data.CD8) %in% topPathways_all$feature]), column_split = meta.data.CD8$zeverity,
                    column_labels = column_labels, row_title = "CD8 TEM",
                    col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 5),
                    show_row_dend = F, show_column_dend = F, cluster_columns = F)





#Adding all heatmaps to one, and saving it:
hm.list <- hm.b %v% hm.cd4 %v% hm.cd8

tiff("./graphs/hm-critical-moderate-b-t-cells-isg-in-deg-pt1-4.tiff",
     res = 300, width = 4, height = 8, units = "in")
draw(hm.list)
dev.off()


#####ISGs in DEGs, for Patients 5 and 6:
#Importing interferon-stimulated gene set:
ISGs <- read.csv("manuscript/ISGs_GeneSets_SU.csv")
ISGs <- ISGs$ISGs_geneset_227
ISGs <- gsub(x = ISGs, pattern = "HLA-", replacement = "HLA.")

#Cleaning gene set:
ISGs <- ISGs[!grepl(ISGs, pattern = "LOC") & !grepl(ISGs, pattern = "XENTR") & !grepl(ISGs, pattern = "GSON")]
ISGs <- ISGs[ISGs != ""]

obj.b <- BCR[,(BCR$patient == "Patient 5" | BCR$patient == "Patient 6") & BCR$severity != "severe"]
obj.b$sample <- droplevels(obj.b$sample)
DefaultAssay(obj.b) <- "RNA"


#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj.b, group_by = "severity", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(desc(abs(logFC)))


# rank genes using rrho algorithm
MNP.genes_B <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                # sort genes by decreasing rank


# choose top upregulated genes
topUp_B_mod <- MNP.genes_B %>% 
  dplyr::filter(group == c('moderate')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>% 
  top_n(50, wt=-padj)

topUp_B_crit <- MNP.genes_B %>% 
  dplyr::filter(group == c('critical')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>%
  top_n(50, wt=-padj)


## moderate vs critical
topPathways_all <- bind_rows(topUp_B_mod, topUp_B_crit)


#Subsetting to only DEGs that are DEGs:
obj.b <- obj.b[ISGs[ISGs %in% topPathways_all$feature],]
obj.b <- NormalizeData(obj.b)

#Changing the data to "pseudo-bulk", using ISG set:
obj.b.tmp <- ScaleData(obj.b)
avg.exp.mat.b <- AverageExpression(obj.b.tmp, group.by = "sample", slot = 'scale.data')
avg.exp.rna.b <- avg.exp.mat.b$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.b.tmp@meta.data %>% 
  group_by(sample, severity) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data.B <- data.frame(t(avg.exp.rna.b), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data.B <- meta.data.B[,c(colnames(meta.data.B)[colnames(meta.data.B) %in% rownames(obj.b.tmp@assays$RNA@scale.data)], "sample", "severity")]


#Preparing aesthetic parameters:
expr.cols <- colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow"))
col.anno <- HeatmapAnnotation(Severity = factor(meta.data.B$severity), show_annotation_name = F,   
                              col = list(Severity = c("moderate" = "blue", "critical" = "red")))
column_labels <- c("critical119_Patient5" = "Pt5", "moderate138_Patient5" = "Pt5",
                   "critical120_Patient6" = "Pt6", "moderate124_Patient6" = "Pt6")

#Preparing separate matrices for each cell type:
rownames(meta.data.B) <- meta.data.B$sample
meta.data.B$severity <- factor(meta.data.B$severity, levels = c("moderate", "critical"))
hm.b <- Heatmap(t(meta.data.B[,colnames(meta.data.B) %in% diff.genes$feature]), 
                column_split = meta.data.B$severity, name = "Expression", 
                column_labels = column_labels, row_title = "B cells", column_title = " ",  
                col = expr.cols, top_annotation = col.anno, row_names_gp = gpar(fontsize = 5), 
                show_row_dend = F, show_column_dend = F, cluster_columns = F)


#Subsetting TCR data:
#TCR <- readRDS("06-TCR-combined.rds")
obj.t <- TCR[,(TCR$patient == "Patient5" | TCR$patient == "Patient6") & TCR$severity != "severe"]


obj.cd4 <- obj.t[,obj.t$azimuthNames == "CD4 TCM"]
obj.cd4$sample <- droplevels(obj.cd4$sample)
DefaultAssay(obj.cd4) <- "RNA"


#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj.cd4, group_by = "severity", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(desc(abs(logFC)))

# rank genes using rrho algorithm
MNP.genes_CD4 <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                # sort genes by decreasing rank


# choose top upregulated genes
topUp_CD4_mod <- MNP.genes_CD4 %>% 
  dplyr::filter(group == c('moderate')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>% 
  top_n(50, wt=-padj)

topUp_CD4_crit <- MNP.genes_CD4 %>% 
  dplyr::filter(group == c('critical')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>%
  top_n(50, wt=-padj)


## moderate vs critical
topPathways_all <- bind_rows(topUp_CD4_mod, topUp_CD4_crit)


#Subsetting to only DEGs that are ISGs:
obj.cd4 <- obj.cd4[ISGs[ISGs %in% topPathways_all$feature],]
obj.cd4 <- NormalizeData(obj.cd4)

#Changing the data to "pseudo-bulk", using ISG set:
obj.cd4.tmp <- ScaleData(obj.cd4)
avg.exp.mat.cd4 <- AverageExpression(obj.cd4.tmp, group.by = "sample", slot = 'scale.data')
avg.exp.rna.cd4 <- avg.exp.mat.cd4$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.cd4.tmp@meta.data %>% 
  group_by(sample, severity) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data.CD4 <- data.frame(t(avg.exp.rna.cd4), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data.CD4 <- meta.data.CD4[,c(colnames(meta.data.CD4)[colnames(meta.data.CD4) %in% rownames(obj.cd4.tmp@assays$RNA@scale.data)], "sample", "severity")]

rownames(meta.data.CD4) <- meta.data.CD4$sample

meta.data.CD4$severity <- factor(meta.data.CD4$severity, levels = c("moderate", "critical"))

hm.cd4 <- Heatmap(t(meta.data.CD4[,colnames(meta.data.CD4) %in% topPathways_all$feature]), column_split = meta.data.CD4$severity,
                  column_labels = column_labels, row_title = "CD4 TCM",
                  col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 5),
                  show_row_dend = F, show_column_dend = F, cluster_columns = F)



obj.cd8 <- obj.t[,obj.t$azimuthNames == "CD8 TEM"]
obj.cd8$sample <- droplevels(obj.cd8$sample)
DefaultAssay(obj.cd8) <- "RNA"


#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj.cd8, group_by = "severity", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(desc(abs(logFC)))

# rank genes using rrho algorithm
MNP.genes_CD8 <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                # sort genes by decreasing rank


# choose top upregulated genes
topUp_CD8_mod <- MNP.genes_CD8 %>% 
  dplyr::filter(group == c('moderate')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>% 
  top_n(50, wt=-padj)

topUp_CD8_crit <- MNP.genes_CD8 %>% 
  dplyr::filter(group == c('critical')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>%
  top_n(50, wt=-padj)


## moderate vs critical
topPathways_all <- bind_rows(topUp_CD8_mod, topUp_CD8_crit)


#Subsetting to only DEGs that are ISGs:
obj.cd8 <- obj.cd8[ISGs[ISGs %in% topPathways_all$feature],]
obj.cd8 <- NormalizeData(obj.cd8)

#Changing the data to "pseudo-bulk", using ISG set:
obj.cd8.tmp <- ScaleData(obj.cd8)
avg.exp.mat.cd8 <- AverageExpression(obj.cd8.tmp, group.by = "sample", slot = 'scale.data')
avg.exp.rna.cd8 <- avg.exp.mat.cd8$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.cd8.tmp@meta.data %>% 
  group_by(sample, severity) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data.CD8 <- data.frame(t(avg.exp.rna.cd8), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data.CD8 <- meta.data.CD8[,c(colnames(meta.data.CD8)[colnames(meta.data.CD8) %in% rownames(obj.cd8.tmp@assays$RNA@scale.data)], "sample", "severity")]

rownames(meta.data.CD8) <- meta.data.CD8$sample
meta.data.CD8$severity <- factor(meta.data.CD8$severity, levels = c("moderate", "critical"))

hm.cd8 <- Heatmap(t(meta.data.CD8[,colnames(meta.data.CD8) %in% topPathways_all$feature]), column_split = meta.data.CD8$severity,
                  column_labels = column_labels, row_title = "CD8 TEM",
                  col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 5),
                  show_row_dend = F, show_column_dend = F, cluster_columns = F)


# width = ncol(t(meta.data.CD8))*unit(5, "mm"), 
# height = nrow(t(meta.data.CD8))*unit(5, "mm")
# 



#Adding all heatmaps to one, and saving it:
hm.list <- hm.b %v% hm.cd4 %v% hm.cd8

tiff("./graphs/hm-critical-moderate-b-t-cells-isg-in-deg-pt5-6.tiff",
     res = 300, width = 4, height = 8, units = "in")
draw(hm.list)
dev.off()





#####The same as above, for all samples:

obj.b <- BCR
DefaultAssay(obj.b) <- "RNA"
obj.b$illness <- "covid"
obj.b$illness[obj.b$severity == "healthy"] <- "control"

#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj.b, group_by = "illness", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(desc(abs(logFC)))


# rank genes using rrho algorithm
MNP.genes_B <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                # sort genes by decreasing rank


# choose top upregulated genes
topUp_B_cov <- MNP.genes_B %>% 
  dplyr::filter(group == c('covid')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>% 
  top_n(50, wt=-padj)

topUp_B_hc <- MNP.genes_B %>% 
  dplyr::filter(group == c('control')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>%
  top_n(50, wt=-padj)


## moderate vs critical
topPathways_all <- bind_rows(topUp_B_cov, topUp_B_hc)


#Subsetting to only DEGs that are DEGs:
obj.b <- obj.b[ISGs[ISGs %in% topPathways_all$feature],]
obj.b <- NormalizeData(obj.b)

#Changing the data to "pseudo-bulk", using ISG set:
obj.b.tmp <- ScaleData(obj.b)
avg.exp.mat.b <- AverageExpression(obj.b.tmp, group.by = "sample", slot = 'scale.data')
avg.exp.rna.b <- avg.exp.mat.b$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.b.tmp@meta.data %>% 
  group_by(sample, severity, illness) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data.B <- data.frame(t(avg.exp.rna.b), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data.B <- meta.data.B[,c(colnames(meta.data.B)[colnames(meta.data.B) %in% rownames(obj.b.tmp@assays$RNA@scale.data)], "sample", "severity", "illness")]


#Preparing aesthetic parameters:
expr.cols <- colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow"))
col.anno <- HeatmapAnnotation(Status = factor(meta.data.B$illness), Severity = factor(meta.data.B$severity), show_annotation_name = F,   
                              col = list(Status = c("control" = "green", "covid" = "#00FFFF"),
                                         Severity = c("mild" = "#000080", "moderate" = "blue",
                                                      "severe" = "pink", "critical" = "red", "healthy" = "#FAFDF3")))
column_labels <- c("healthy1_control1" = "HC1", "healthy2_control2" = "HC2", "healthy3_control3" = "HC3",
                   "moderate272_Patient1" = "Pt1_m", "critical293_Patient1" = "Pt1_c",
                   "moderate303_Patient2" = "Pt2_m", "critical308_Patient2" = "Pt2_c",
                   "mild186_Patient3" = "Pt3_m", "critical213_Patient3" = "Pt3_c",
                   "mild227_Patient4" = "Pt4_m", "critical238_Patient4" = "Pt4_c",
                   "critical119_Patient5" = "Pt5_c", "severe123_Patient5" = "Pt5_s", "moderate138_Patient5" = "Pt5_m",
                   "critical120_Patient6" = "Pt6_c", "severe122_Patient6" = "Pt6_s", "moderate124_Patient6" = "Pt6_m")


#Preparing separate matrices for each cell type:
hm.b <- Heatmap(t(meta.data.B[,colnames(meta.data.B) %in% diff.genes$feature]), 
                column_split = meta.data.B$illness, name = "Expression", 
                # column_order = c("healthy1_control1", "healthy2_control2", "healthy3_control3",
                #                  "moderate272_Patient1", "moderate303_Patient2", "mild186_Patient3", "mild227_Patient4",
                #                  "moderate138_Patient5", "moderate124_Patient6", "severe123_Patient5", "severe122_Patient6",
                #                  "critical293_Patient1", "critical308_Patient2", "critical213_Patient3", "critical238_Patient4",
                #                  "critical119_Patient5", "critical120_Patient6"),                        #If we need a certain column order
                column_labels = column_labels, row_title = "B cells", column_title = " ",  
                col = expr.cols, top_annotation = col.anno, row_names_gp = gpar(fontsize = 5), 
                show_row_dend = F, show_column_dend = F, cluster_columns = F)


#Subsetting TCR data:
TCR <- readRDS("06-TCR-combined.rds")
obj.t <- TCR[, TCR$sample != "healthy4_control4"]
obj.t$sample <- factor(obj.t$sample)
obj.t$sample <- droplevels(obj.t$sample)

obj.cd4 <- obj.t[,obj.t$azimuthNames == "CD4 TCM"]
DefaultAssay(obj.cd4) <- "RNA"
obj.cd4$illness <- "covid"
obj.cd4$illness[obj.cd4$severity == "healthy"] <- "control"

#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj.cd4, group_by = "illness", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(desc(abs(logFC)))

# rank genes using rrho algorithm
MNP.genes_CD4 <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                # sort genes by decreasing rank


# choose top upregulated genes
topUp_CD4_cov <- MNP.genes_CD4 %>% 
  dplyr::filter(group == c('covid')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>% 
  top_n(50, wt=-padj)

topUp_CD4_hc <- MNP.genes_CD4 %>% 
  dplyr::filter(group == c('control')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>%
  top_n(50, wt=-padj)


## moderate vs critical
topPathways_all <- bind_rows(topUp_CD4_cov, topUp_CD4_hc)


#Subsetting to only DEGs that are ISGs:
obj.cd4 <- obj.cd4[ISGs[ISGs %in% topPathways_all$feature],]
obj.cd4 <- NormalizeData(obj.cd4)

#Changing the data to "pseudo-bulk", using ISG set:
obj.cd4.tmp <- ScaleData(obj.cd4)
avg.exp.mat.cd4 <- AverageExpression(obj.cd4.tmp, group.by = "sample", slot = 'scale.data')
avg.exp.rna.cd4 <- avg.exp.mat.cd4$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.cd4.tmp@meta.data %>% 
  group_by(sample, severity, illness) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data.CD4 <- data.frame(t(avg.exp.rna.cd4), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data.CD4 <- meta.data.CD4[,c(colnames(meta.data.CD4)[colnames(meta.data.CD4) %in% rownames(obj.cd4.tmp@assays$RNA@scale.data)], "sample", "severity", "illness")]


hm.cd4 <- Heatmap(t(meta.data.CD4[,colnames(meta.data.CD4) %in% topPathways_all$feature]), column_split = meta.data.CD4$illness,
                  column_labels = column_labels, row_title = "CD4 TCM", 
                  col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 5), 
                  show_row_dend = F, show_column_dend = F, cluster_columns = F)



obj.cd8 <- obj.t[,obj.t$azimuthNames == "CD8 TEM"]
DefaultAssay(obj.cd8) <- "RNA"
obj.cd8$illness <- "covid"
obj.cd8$illness[obj.cd8$severity == "healthy"] <- "control"

#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj.cd8, group_by = "illness", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(desc(abs(logFC)))

# rank genes using rrho algorithm
MNP.genes_CD8 <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                # sort genes by decreasing rank


# choose top upregulated genes
topUp_CD8_cov <- MNP.genes_CD8 %>% 
  dplyr::filter(group == c('covid')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>% 
  top_n(50, wt=-padj)

topUp_CD8_hc <- MNP.genes_CD8 %>% 
  dplyr::filter(group == c('control')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>%
  top_n(50, wt=-padj)


## moderate vs critical
topPathways_all <- bind_rows(topUp_CD8_cov, topUp_CD8_hc)

#Subsetting to only DEGs that are ISGs:
obj.cd8 <- obj.cd8[ISGs[ISGs %in% topPathways_all$feature],]
obj.cd8 <- NormalizeData(obj.cd8)

#Changing the data to "pseudo-bulk", using ISG set:
obj.cd8.tmp <- ScaleData(obj.cd8)
avg.exp.mat.cd8 <- AverageExpression(obj.cd8.tmp, group.by = "sample", slot = 'scale.data')
avg.exp.rna.cd8 <- avg.exp.mat.cd8$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.cd8.tmp@meta.data %>% 
  group_by(sample, severity, illness) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data.CD8 <- data.frame(t(avg.exp.rna.cd8), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data.CD8 <- meta.data.CD8[,c(colnames(meta.data.CD8)[colnames(meta.data.CD8) %in% rownames(obj.cd8.tmp@assays$RNA@scale.data)], "sample", "severity", "illness")]


hm.cd8 <- Heatmap(t(meta.data.CD8[,colnames(meta.data.CD8) %in% topPathways_all$feature]), column_split = meta.data.CD8$illness,
                  column_labels = column_labels, row_title = "CD8 TEM", 
                  col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 5),
                  show_row_dend = F, show_column_dend = F, cluster_columns = F)






#Adding all heatmaps to one, and saving it:
hm.list <- hm.b %v% hm.cd4 %v% hm.cd8

tiff("./graphs/hm-all-b-t-cells-isg-in-deg.tiff",
     res = 300, width = 4, height = 8, units = "in")
draw(hm.list)
dev.off()










#####Dividing B cells:
#Importing interferon-stimulated gene set:
ISGs <- read.csv("manuscript/ISGs_GeneSets_SU.csv")
ISGs <- union(ISGs$ISGs_0034340_response_to_type1_interferon, 
              c(ISGs$ISGs_0034341_response_to_interferon_gamma, ISGs$ISGs_0035455_response_to_interferon_alpha, ISGs$ISGs_geneset_227))

#Cleaning gene set:
ISGs <- ISGs[!grepl(ISGs, pattern = "LOC") & !grepl(ISGs, pattern = "XENTR") & !grepl(ISGs, pattern = "GSON")]
ISGs <- ISGs[ISGs != ""]

obj <- BCR[,(BCR$patient == "Patient 5" | BCR$patient == "Patient 6") & BCR$severity != "severe"]
obj$sample <- droplevels(obj$sample)
obj$zeverity <- "critical"
obj$zeverity[obj$severity == "mild" | obj$severity == "moderate"] <- "moderate"
DefaultAssay(obj) <- "RNA"
obj <- obj[ISGs,]
obj <- NormalizeData(obj)

#Changing the data to "pseudo-bulk", using ISG set:
obj.tmp <- ScaleData(obj, features = ISGs)
avg.exp.mat <- AverageExpression(obj.tmp, features = ISGs, group.by = c("sample", "azimuthNames"),  slot = 'scale.data')
avg.exp.rna <- avg.exp.mat$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.tmp@meta.data %>% 
  group_by(sample, zeverity, azimuthNames) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data <- data.frame(t(avg.exp.rna), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data <- meta.data[,c(colnames(meta.data)[colnames(meta.data) %in% rownames(obj.tmp@assays$RNA@scale.data)], "sample", "azimuthNames", "zeverity")]

#Preparing aesthetic parameters:
expr.cols <- colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow"))
col.anno <- HeatmapAnnotation(Severity = factor(meta.data.B.int$zeverity), show_annotation_name = F,   
                              col = list(Severity = c("critical" = "red", "moderate" = "blue")))
column_labels <- c("critical119_Patient5" = "Pt5", "moderate138_Patient5" = "Pt5",
                   "critical120_Patient6" = "Pt6", "moderate124_Patient6" = "Pt6")

#Preparing separate matrices for each cell type:
meta.data.B.int <- meta.data[meta.data$azimuthNames == "B intermediate",]
rownames(meta.data.B.int) <- meta.data.B.int$sample
hm.b.int <- Heatmap(t(meta.data.B.int[,1:30]), column_split = meta.data.B.int$zeverity, name = "Expression",
                    column_labels = column_labels, row_title = "B intermediate", column_title = " ", 
                    col = expr.cols, top_annotation = col.anno, row_names_gp = gpar(fontsize = 5), 
                    show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F)


meta.data.B.mem <- meta.data[meta.data$azimuthNames == "B memory",]
rownames(meta.data.B.mem) <- meta.data.B.mem$sample
hm.b.mem <- Heatmap(t(meta.data.B.mem[,1:30]), column_split = meta.data.B.int$zeverity,
                    column_labels = column_labels, row_title = "B memory",
                    col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 5),
                    show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F)

meta.data.B.nai <- meta.data[meta.data$azimuthNames == "B naive",]
rownames(meta.data.B.nai) <- meta.data.B.nai$sample
hm.b.nai <- Heatmap(t(meta.data.B.nai[,1:30]), column_split = meta.data.B.int$zeverity,
                    column_labels = column_labels, row_title = "B naive",
                    col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 5),
                    show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F)

meta.data.Plasm <- meta.data[meta.data$azimuthNames == "Plasmablast",]
rownames(meta.data.Plasm) <- meta.data.Plasm$sample
hm.plasm <- Heatmap(t(meta.data.Plasm[,1:30]), column_split = meta.data.B.int$zeverity,
                    column_labels = column_labels, row_title = "Plasmablasts",
                    col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 5),
                    show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F)

#Adding all heatmaps to one, and saving it:
hm.list <- hm.b.int %v% hm.b.mem %v% hm.b.nai %v% hm.plasm

tiff("./graphs/hm-critical-moderate-cells-isgs-pt5-6.tiff",
     res = 300, width = 4, height = 8, units = "in")
draw(hm.list)
dev.off()



#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj, group_by = "zeverity", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(padj)

#Changing the data to "pseudo-bulk", using DEG set:
obj.tmp <- ScaleData(obj, features = diff.genes$feature)
avg.exp.mat <- AverageExpression(obj.tmp, features = diff.genes$feature, group.by = c("sample", "azimuthNames"),  slot = 'scale.data')
avg.exp.rna <- avg.exp.mat$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.tmp@meta.data %>% 
  group_by(sample, zeverity, azimuthNames) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data <- data.frame(t(avg.exp.rna), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data <- meta.data[,c(colnames(meta.data)[colnames(meta.data) %in% rownames(obj.tmp@assays$RNA@scale.data)], "sample", "azimuthNames", "zeverity")]

#Preparing aesthetic parameters:
expr.cols <- colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow"))
col.anno <- HeatmapAnnotation(Severity = factor(meta.data.B.int$zeverity), show_annotation_name = F,   
                              col = list(Severity = c("critical" = "red", "moderate" = "blue")))
column_labels <- c("critical119_Patient5" = "Pt5", "moderate138_Patient5" = "Pt5",
                   "critical120_Patient6" = "Pt6", "moderate124_Patient6" = "Pt6")

#Preparing separate matrices for each cell type:
meta.data.B.int <- meta.data[meta.data$azimuthNames == "B intermediate",]
rownames(meta.data.B.int) <- meta.data.B.int$sample
hm.b.int <- Heatmap(t(meta.data.B.int[,1:30]), column_split = meta.data.B.int$zeverity, name = "Expression",
                    column_labels = column_labels, row_title = "B intermediate", column_title = " ", 
                    col = expr.cols, top_annotation = col.anno, row_names_gp = gpar(fontsize = 5), 
                    show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F)


meta.data.B.mem <- meta.data[meta.data$azimuthNames == "B memory",]
rownames(meta.data.B.mem) <- meta.data.B.mem$sample
hm.b.mem <- Heatmap(t(meta.data.B.mem[,1:30]), column_split = meta.data.B.int$zeverity,
                    column_labels = column_labels, row_title = "B memory",
                    col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 5),
                    show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F)

meta.data.B.nai <- meta.data[meta.data$azimuthNames == "B naive",]
rownames(meta.data.B.nai) <- meta.data.B.nai$sample
hm.b.nai <- Heatmap(t(meta.data.B.nai[,1:30]), column_split = meta.data.B.int$zeverity,
                    column_labels = column_labels, row_title = "B naive",
                    col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 5),
                    show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F)

meta.data.Plasm <- meta.data[meta.data$azimuthNames == "Plasmablast",]
rownames(meta.data.Plasm) <- meta.data.Plasm$sample
hm.plasm <- Heatmap(t(meta.data.Plasm[,1:30]), column_split = meta.data.B.int$zeverity,
                    column_labels = column_labels, row_title = "Plasmablasts",
                    col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 5),
                    show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F)

#Adding all heatmaps to one, and saving it:
hm.list <- hm.b.int %v% hm.b.mem %v% hm.b.nai %v% hm.plasm

tiff("./graphs/hm-critical-moderate-cells-deg-pt5-6.tiff",
     res = 300, width = 4, height = 8, units = "in")
draw(hm.list)
dev.off()




