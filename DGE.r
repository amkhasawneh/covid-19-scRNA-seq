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

#Choosing the human Hallmark annotated gene set from the MsigDB:
m_df <- msigdbr(species = "Homo sapiens", category = "H")


#Preparing a list of the gene sets, for fgsea:
fgsea_sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
names(fgsea_sets) <- gsub(pattern = "HALLMARK_", replacement = "", x = names(fgsea_sets)) |> gsub(pattern = "_", replacement = " ")


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

rm(list=setdiff(ls(), c("BCR", "fgsea_sets", "m_df")))

#For Patient 1:
pt1 <- BCR[,BCR$patient == "Patient 1"]
pt1$sample <- droplevels(pt1$sample)
DefaultAssay(pt1) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto

wlx.mrk.pt1 <- wilcoxauc(pt1, 'sample', 
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

top <- topPathways %>% group_by(group) 

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

enrich.topUpmod1[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpmod1[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)
enrich.topUpmod1[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpmod1[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)


#Selecting only the feature and rank columns of DGE data for fgsea run:
MNP.genes_c1 <- ranked.genes %>%
  dplyr::filter(group == "critical293_Patient1") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)

MNP.genes_m1 <- ranked.genes %>%
  dplyr::filter(group == "moderate272_Patient1") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)


#Converting matrices to tables:
ranks_c <- deframe(MNP.genes_c1)
ranks_m <- deframe(MNP.genes_m1)

#Running fgsea based on gene set ranking (stats = ranks):
fgseaRes_c1 <- fgseaMultilevel(fgsea_sets, stats = ranks_c, nPermSimple = 1000, 
                               maxSize = 200) %>% arrange(padj)
fgseaRes_m1 <- fgseaMultilevel(fgsea_sets, stats = ranks_m, nPermSimple = 1000, 
                               maxSize = 200) %>% arrange(padj)



#Enrichment plots:
ggsave(filename = "graphs/enrichment-plot-pt1-critical-moderate.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_c1, desc(NES))[1]$pathway]],
                             ranks_c) + 
         labs(title = arrange(fgseaRes_c1, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))

ggsave(filename = "graphs/enrichment-plot-pt1-moderate-critical.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_m1, desc(NES))[1]$pathway]],
                             ranks_m) + 
         labs(title = arrange(fgseaRes_m1, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))



#Plotting pathways:
fgseaRes_c1$adjPvalue <- ifelse(fgseaRes_c1$padj <= 0.05, "significant", "non-significant")
cols <- c("significant" = "red") #"non-significant" = "grey", 
plot_c1 <- ggplot(fgseaRes_c1, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = "Patient 1 - critical vs moderate") + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='right', legend.key.size = unit(0.25, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) 
ggsave(filename = "graphs/gsea-pt1-critical-vs-moderate-top.tiff",
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_c1)

all1c <- data.frame(patient = "Patient 1", severity = "critical", progress = "progressing", outcome = "deceased", pathway = fgseaRes_c1$pathway, size = fgseaRes_c1$size,
                  NES = fgseaRes_c1$NES, padj = fgseaRes_c1$padj) %>% filter(padj < 0.05) %>% arrange(desc(NES))
all1m <- data.frame(patient = "Patient 1", severity = "moderate", progress = "progressing", outcome = "deceased", pathway = fgseaRes_m1$pathway, size = fgseaRes_m1$size,
                  NES = fgseaRes_m1$NES, padj = fgseaRes_m1$padj) %>% filter(padj < 0.05) %>% arrange(desc(NES))

allc <- NULL
allm <- NULL
allc <- rbind(allc, first(all1c[all1c$NES > 0,], 5), last(all1c[all1c$NES < 0,], 5))
allm <- rbind(allm, first(all1m[all1m$NES > 0,], 5), last(all1m[all1m$NES < 0,], 5))

rm(list=setdiff(ls(), c("BCR", "fgsea_sets", "m_df", "allc", "allm")))

#For Patient 2:
pt2 <- BCR[,BCR$patient == "Patient 2"]
pt2$sample <- droplevels(pt2$sample)
DefaultAssay(pt2) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto

wlx.mrk.pt2 <- wilcoxauc(pt2, 'sample',
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



enrich.topUpmod2[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpmod2[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)
enrich.topUpcrit2[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpcrit2[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)

#Selecting only the feature and rank columns of DGE data for fgsea run:
MNP.genes_c2 <- ranked.genes %>%
  dplyr::filter(group == "critical308_Patient2") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)

MNP.genes_m2 <- ranked.genes %>%
  dplyr::filter(group == "moderate303_Patient2") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)


#Converting matrices to tables:
ranks_c <- deframe(MNP.genes_c2)
ranks_m <- deframe(MNP.genes_m2)

#Running fgsea based on gene set ranking (stats = ranks):
fgseaRes_c2 <- fgseaMultilevel(fgsea_sets, stats = ranks_c, nPermSimple = 1000, 
                               maxSize = 200) %>% arrange(padj)
fgseaRes_m2 <- fgseaMultilevel(fgsea_sets, stats = ranks_m, nPermSimple = 1000, 
                               maxSize = 200) %>% arrange(padj)



#Enrichment plots:
ggsave(filename = "graphs/enrichment-plot-pt2-critical-moderate.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_c2, desc(NES))[1]$pathway]],
                             ranks_c) + 
         labs(title = arrange(fgseaRes_c2, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))

ggsave(filename = "graphs/enrichment-plot-pt2-moderate-critical.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_m2, desc(NES))[1]$pathway]],
                             ranks_m) + 
         labs(title = arrange(fgseaRes_m2, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))



#Plotting pathways:
fgseaRes_c2$adjPvalue <- ifelse(fgseaRes_c2$padj <= 0.05, "significant", "non-significant")
cols <- c("significant" = "red") #"non-significant" = "grey", 
plot_c2 <- ggplot(fgseaRes_c2, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = "Patient 2 - critical vs moderate") + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='right', legend.key.size = unit(0.25, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) 
ggsave(filename = "graphs/gsea-pt2-critical-vs-moderate-top.tiff",
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_c2)

all2c <- data.frame(patient = "Patient 2", severity = "critical", progress = "progressing", outcome = "deceased", pathway = fgseaRes_c2$pathway, size = fgseaRes_c2$size,
                    NES = fgseaRes_c2$NES, padj = fgseaRes_c2$padj) %>% filter(padj < 0.05) %>% arrange(desc(NES))
all2m <- data.frame(patient = "Patient 2", severity = "moderate", progress = "progressing", outcome = "deceased", pathway = fgseaRes_m2$pathway, size = fgseaRes_m2$size,
                    NES = fgseaRes_m2$NES, padj = fgseaRes_m2$padj) %>% filter(padj < 0.05) %>% arrange(desc(NES))


allc <- rbind(allc, first(all2c[all2c$NES > 0,], 5), last(all2c[all2c$NES < 0,], 5))
allm <- rbind(allm, first(all2m[all2m$NES > 0,], 5), last(all2m[all2m$NES < 0,], 5))

rm(list=setdiff(ls(), c("BCR", "fgsea_sets", "m_df", "allc", "allm")))

#For Patient 3:
pt3 <- BCR[,BCR$patient == "Patient 3"]
pt3$sample <- droplevels(pt3$sample)
DefaultAssay(pt3) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto

wlx.mrk.pt3 <- wilcoxauc(pt3, 'sample', 
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


enrich.topUpmod3[["GO_Biological_Process_2021"]]$Term <-  vapply(str_split(enrich.topUpmod3[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)
enrich.topUpcrit3[["GO_Biological_Process_2021"]]$Term <-  vapply(str_split(enrich.topUpcrit3[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)

#Selecting only the feature and rank columns of DGE data for fgsea run:
MNP.genes_c3 <- ranked.genes %>%
  dplyr::filter(group == "critical213_Patient3") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)

MNP.genes_m3 <- ranked.genes %>%
  dplyr::filter(group == "mild186_Patient3") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)


#Converting matrices to tables:
ranks_c <- deframe(MNP.genes_c3)
ranks_m <- deframe(MNP.genes_m3)

#Running fgsea based on gene set ranking (stats = ranks):
fgseaRes_c3 <- fgseaMultilevel(fgsea_sets, stats = ranks_c, nPermSimple = 1000, 
                               maxSize = 200) %>% arrange(padj)
fgseaRes_m3 <- fgseaMultilevel(fgsea_sets, stats = ranks_m, nPermSimple = 1000, 
                               maxSize = 200) %>% arrange(padj)


#Enrichment plots:
ggsave(filename = "graphs/enrichment-plot-pt3-critical-mild.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_c3, desc(NES))[1]$pathway]],
                             ranks_c) + 
         labs(title = arrange(fgseaRes_c3, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))

ggsave(filename = "graphs/enrichment-plot-pt3-mild-critical.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_m3, desc(NES))[1]$pathway]],
                             ranks_m) + 
         labs(title = arrange(fgseaRes_m3, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))



#Plotting pathways:
fgseaRes_c3$adjPvalue <- ifelse(fgseaRes_c3$padj <= 0.05, "significant", "non-significant")
cols <- c("significant" = "red") #"non-significant" = "grey", 
plot_c3 <- ggplot(fgseaRes_c3, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = "Patient 3 - critical vs mild") + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='right', legend.key.size = unit(0.25, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) 
ggsave(filename = "graphs/gsea-pt3-critical-vs-mild-top.tiff",
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_c3)


all3c <- data.frame(patient = "Patient 3", severity = "critical", progress = "progressing", outcome = "deceased", pathway = fgseaRes_c3$pathway, size = fgseaRes_c3$size,
                    NES = fgseaRes_c3$NES, padj = fgseaRes_c3$padj) %>% filter(padj < 0.05) %>% arrange(desc(NES))
all3m <- data.frame(patient = "Patient 3", severity = "mild", progress = "progressing", outcome = "deceased", pathway = fgseaRes_m3$pathway, size = fgseaRes_m3$size,
                    NES = fgseaRes_m3$NES, padj = fgseaRes_m3$padj) %>% filter(padj < 0.05)%>% arrange(desc(NES))


allc <- rbind(allc, first(all3c[all3c$NES > 0,], 5), last(all3c[all3c$NES < 0,], 5))
allm <- rbind(allm, first(all3m[all3m$NES > 0,], 5), last(all3m[all3m$NES < 0,], 5))

rm(list=setdiff(ls(), c("BCR", "fgsea_sets", "m_df", "allc", "allm")))



#For Patient 4:
pt4 <- BCR[,BCR$patient == "Patient 4"]
pt4$sample <- droplevels(pt4$sample)
DefaultAssay(pt4) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto

wlx.mrk.pt4 <- wilcoxauc(pt4, 'sample',
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

enrich.topUpmod4[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpmod4[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)
enrich.topUpcrit4[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpcrit4[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)


#Selecting only the feature and rank columns of DGE data for fgsea run:
MNP.genes_c <- ranked.genes %>%
  dplyr::filter(group == "critical238_Patient4") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)

MNP.genes_m <- ranked.genes %>%
  dplyr::filter(group == "mild227_Patient4") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)


#Converting matrices to tables:
ranks_c <- deframe(MNP.genes_c)
ranks_m <- deframe(MNP.genes_m)

#Running fgsea based on gene set ranking (stats = ranks):
fgseaRes_c4 <- fgseaMultilevel(fgsea_sets, stats = ranks_c, nPermSimple = 1000, 
                                maxSize = 200) %>% arrange(padj)
fgseaRes_m4 <- fgseaMultilevel(fgsea_sets, stats = ranks_m, nPermSimple = 1000, 
                                maxSize = 200) %>% arrange(padj)



#Enrichment plots:
ggsave(filename = "graphs/enrichment-plot-pt4-critical-mild.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_c4, desc(NES))[1]$pathway]],
                             ranks_c) + 
         labs(title = arrange(fgseaRes_c4, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))

ggsave(filename = "graphs/enrichment-plot-pt5-mild-critical.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_m4, desc(NES))[1]$pathway]],
                             ranks_m) + 
         labs(title = arrange(fgseaRes_m4, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))



#Plotting pathways:
fgseaRes_c4$adjPvalue <- ifelse(fgseaRes_c4$padj <= 0.05, "significant", "non-significant")
cols <- c("significant" = "red") #"non-significant" = "grey", 
plot_c4 <- ggplot(fgseaRes_c4, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = "Patient 4 - critical vs mild") + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='right', legend.key.size = unit(0.25, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) 
ggsave(filename = "graphs/gsea-pt4-critical-vs-mild-top.tiff",
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_c4)


all4c <- data.frame(patient = "Patient 4", severity = "critical", progress = "progressing", outcome = "survived", pathway = fgseaRes_c4$pathway, size = fgseaRes_c4$size,
                    NES = fgseaRes_c4$NES, padj = fgseaRes_c4$padj) %>% filter(padj < 0.05) %>% arrange(desc(NES))
all4m <- data.frame(patient = "Patient 4", severity = "mild", progress = "progressing", outcome = "survived", pathway = fgseaRes_m4$pathway, size = fgseaRes_m4$size,
                    NES = fgseaRes_m4$NES, padj = fgseaRes_m4$padj) %>% filter(padj < 0.05) %>% arrange(desc(NES))


allc <- rbind(allc, first(all4c[all4c$NES > 0,], 5), last(all4c[all4c$NES < 0,], 5))
allm <- rbind(allm, first(all4m[all4m$NES > 0,], 5), last(all4m[all4m$NES < 0,], 5))

rm(list=setdiff(ls(), c("BCR", "fgsea_sets", "m_df", "allc", "allm")))

#For Patient 5:
pt5 <- BCR[,BCR$patient == "Patient 5"]
pt5$sample <- droplevels(pt5$sample)
DefaultAssay(pt5) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto
#, "moderate138_Patient5"

wlx.mrk.pt5.cs <- wilcoxauc(pt5, 'sample',
                         c("critical119_Patient5", "severe123_Patient5"),
                               seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(padj)

wlx.mrk.pt5.cm <- wilcoxauc(pt5, 'sample',
                         c("critical119_Patient5", "moderate138_Patient5"),
                               seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(padj)

wlx.mrk.pt5.ms <- wilcoxauc(pt5, 'sample',
                         c("moderate138_Patient5", "severe123_Patient5"),
                               seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(padj)

head(wlx.mrk.pt5)

# we have all the genes for each cluster
dplyr::count(wlx.mrk.pt5, group)

# summarize the top abundant marker features for each group
top_markers(wlx.mrk.pt5)

write.table(top_markers(wlx.mrk.pt5), file = "top-wilcox-markers-pt5.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes.cs <- wlx.mrk.pt5.cs %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank

ranked.genes.cm <- wlx.mrk.pt5.cm %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank

ranked.genes.ms <- wlx.mrk.pt5.ms %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank


# save found markers as exportable table
write.table(ranked.genes.cm,
            file = "ranked-genes-pt5-cm.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUpcritsev <- ranked.genes.cs %>% 
  dplyr::filter(group == "critical119_Patient5") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)
topUpcritmod <- ranked.genes.cm %>% 
  dplyr::filter(group == "critical119_Patient5") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topUpsevmod <- ranked.genes.ms %>% 
  dplyr::filter(group == "severe123_Patient5") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)
topUpsevcrit <- ranked.genes.cs %>% 
  dplyr::filter(group == "severe123_Patient5") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topUpmodcrit <- ranked.genes.cm %>% 
  dplyr::filter(group == "moderate138_Patient5") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)
topUpmodsev <- ranked.genes.ms %>% 
  dplyr::filter(group == "moderate138_Patient5") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)


topPathways <- bind_rows(topUpcritmod, topUpcritsev, topUpsevcrit, topUpsevmod, topUpmodcrit, topUpmodsev)

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
enrich.topUpcrit5.cs <- enrichr(genes = topUpcritsev$feature, 
                               databases = "GO_Biological_Process_2021")
enrich.topUpcrit5.cm <- enrichr(genes = topUpcritmod$feature, 
                               databases = "GO_Biological_Process_2021")
enrich.topUpsev5.sc <- enrichr(genes = topUpsevcrit$feature, 
                               databases = "GO_Biological_Process_2021")
enrich.topUpsev5.sm <- enrichr(genes = topUpsevmod$feature, 
                               databases = "GO_Biological_Process_2021")
enrich.topUpmod5.mc <- enrichr(genes = topUpmodcrit$feature, 
                               databases = "GO_Biological_Process_2021")
enrich.topUpmod5.ms <- enrichr(genes = topUpmodsev$feature, 
                               databases = "GO_Biological_Process_2021")

enrich.topUpmod5.mc[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpmod5.mc[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)
enrich.topUpmod5.ms[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpmod5.ms[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)
enrich.topUpsev5.cs[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpsev5.cs[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)
enrich.topUpsev5.ms[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpsev5.ms[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)
enrich.topUpcrit5.cs[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpcrit5.cs[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)
enrich.topUpcrit5.cm[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpcrit5.cm[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)


ggsave(filename = "graphs/enrichment-pt5-wlx-crit-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit5.cm[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 5 critical vs moderate"))
ggsave(filename = "graphs/enrichment-pt5-wlx-crit-sev.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit5.cs[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 5 critical vs severe"))

ggsave(filename = "graphs/enrichment-pt5-wlx-sev-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpsev5.ms[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 5 severe vs moderate"))
ggsave(filename = "graphs/enrichment-pt5-wlx-sev-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpsev5.cs[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 5 severe vs critical"))

ggsave(filename = "graphs/enrichment-pt5-wlx-mod-sev.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod5.ms[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 5 moderate vs severe"))
ggsave(filename = "graphs/enrichment-pt5-wlx-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod5.mc[[1]], 
           numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
           xlab = NULL,
           ylab = NULL,
           title = "Upregulated in Patient 5 moderate vs critical"))

#Selecting only the feature and rank columns of DGE data for fgsea run:
MNP.genes_cs <- ranked.genes.cs %>%
  dplyr::filter(group == "critical119_Patient5") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)
MNP.genes_cm <- ranked.genes.cm %>%
  dplyr::filter(group == "critical119_Patient5") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)

MNP.genes_sc <- ranked.genes.cs %>%
  dplyr::filter(group == "severe123_Patient5") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)
MNP.genes_sm <- ranked.genes.ms %>%
  dplyr::filter(group == "severe123_Patient5") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)

MNP.genes_mc <- ranked.genes.cm %>%
  dplyr::filter(group == "moderate138_Patient5") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)
MNP.genes_ms <- ranked.genes.ms %>%
  dplyr::filter(group == "moderate138_Patient5") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)


#Converting matrices to tables:
ranks_cs <- deframe(MNP.genes_cs)
ranks_cm <- deframe(MNP.genes_cm)
ranks_sc <- deframe(MNP.genes_sc)
ranks_sm <- deframe(MNP.genes_sm)
ranks_mc <- deframe(MNP.genes_mc)
ranks_ms <- deframe(MNP.genes_ms)

#Running fgsea based on gene set ranking (stats = ranks):
fgseaRes_cs5 <- fgseaMultilevel(fgsea_sets, stats = ranks_cs, nPermSimple = 1000, 
                               maxSize = 200) %>% arrange(padj)
fgseaRes_cm5 <- fgseaMultilevel(fgsea_sets, stats = ranks_cm, nPermSimple = 1000, 
                               maxSize = 200) %>% arrange(padj)
fgseaRes_sc5 <- fgseaMultilevel(fgsea_sets, stats = ranks_sc, nPermSimple = 1000, 
                               maxSize = 200) %>% arrange(padj)
fgseaRes_sm5 <- fgseaMultilevel(fgsea_sets, stats = ranks_sm, nPermSimple = 1000, 
                               maxSize = 200) %>% arrange(padj)
fgseaRes_mc5 <- fgseaMultilevel(fgsea_sets, stats = ranks_mc, nPermSimple = 1000, 
                               maxSize = 200) %>% arrange(padj)
fgseaRes_ms5 <- fgseaMultilevel(fgsea_sets, stats = ranks_ms, nPermSimple = 1000, 
                               maxSize = 200) %>% arrange(padj)


#Enrichment plots:
ggsave(filename = "graphs/enrichment-plot-pt5-critical-severe.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_cs, desc(NES))[1]$pathway]],
                             ranks_cs) + 
         labs(title = arrange(fgseaRes_cs, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))
ggsave(filename = "graphs/enrichment-plot-pt5-critical-moderate.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_cm, desc(NES))[1]$pathway]],
                             ranks_cm) + 
         labs(title = arrange(fgseaRes_cm, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))

ggsave(filename = "graphs/enrichment-plot-pt5-severe-critical.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_sc, desc(NES))[1]$pathway]],
                             ranks_sc) + 
         labs(title = arrange(fgseaRes_sc, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))
ggsave(filename = "graphs/enrichment-plot-pt5-severe-moderate.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_sm, desc(NES))[1]$pathway]],
                             ranks_sm) + 
         labs(title = arrange(fgseaRes_sm, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))

ggsave(filename = "graphs/enrichment-plot-pt5-moderate-critical.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_mc, desc(NES))[1]$pathway]],
                             ranks_mc) + 
         labs(title = arrange(fgseaRes_mc, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))
ggsave(filename = "graphs/enrichment-plot-pt5-moderate-severe.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_ms, desc(NES))[1]$pathway]],
                             ranks_ms) + 
         labs(title = arrange(fgseaRes_ms, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))



#Plotting pathways:
fgseaRes_cs5$adjPvalue <- ifelse(fgseaRes_cs5$padj <= 0.05, "significant", "non-significant")
cols <- c("significant" = "red") #"non-significant" = "grey", 
plot_cs <- ggplot(fgseaRes_cs5, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = "Patient 5 - critical vs severe") + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='right', legend.key.size = unit(0.25, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) 
ggsave(filename = "graphs/gsea-pt5-critical-vs-severe-top.tiff",
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_cs)


fgseaRes_cm5$adjPvalue <- ifelse(fgseaRes_cm5$padj <= 0.05, "significant", "non-significant")
cols <- c("significant" = "red") #"non-significant" = "grey", 
plot_cm <- ggplot(fgseaRes_cm5, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = "Patient 5 - critical vs moderate") + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='right', legend.key.size = unit(0.25, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) 
ggsave(filename = "graphs/gsea-pt5-critical-vs-moderate-top.tiff",
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_cm)

fgseaRes_sc5$adjPvalue <- ifelse(fgseaRes_sc5$padj <= 0.05, "significant", "non-significant")
plot_sc <- ggplot(fgseaRes_sc5, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = "Patient 5 - severe vs critical") + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='right', legend.key.size = unit(0.25, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) 
ggsave(filename = "graphs/gsea-pt5-severe-vs-critical-top.tiff",
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_sc)

fgseaRes_sm5$adjPvalue <- ifelse(fgseaRes_sm5$padj <= 0.05, "significant", "non-significant")
plot_sm <- ggplot(fgseaRes_sm5, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = "Patient 5 - severe vs moderate") + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='right', legend.key.size = unit(0.25, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) 
ggsave(filename = "graphs/gsea-pt5-severe-vs-moderate-top.tiff",
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_sm)



fgseaRes_mc5$adjPvalue <- ifelse(fgseaRes_mc5$padj <= 0.05, "significant", "non-significant")
plot_mc <- ggplot(fgseaRes_mc5, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = "Patient 5 - moderate vs critical") + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='right', legend.key.size = unit(0.25, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) 
ggsave(filename = "graphs/gsea-pt5-moderatee-vs-critical-top.tiff",
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_mc)

fgseaRes_ms5$adjPvalue <- ifelse(fgseaRes_ms5$padj <= 0.05, "significant", "non-significant")
plot_ms <- ggplot(fgseaRes_ms, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = "Patient 5 - moderate vs severe") + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='right', legend.key.size = unit(0.25, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) 
ggsave(filename = "graphs/gsea-pt5-moderatee-vs-severe-top.tiff",
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_ms)




all5c <- data.frame(patient = "Patient 5", severity = "critical", progress = "recovering", outcome = "survived", pathway = fgseaRes_cm5$pathway, size = fgseaRes_cm5$size,
                    NES = fgseaRes_cm5$NES, padj = fgseaRes_cm5$padj) %>% filter(padj < 0.05) %>% arrange(desc(NES))
all5m <- data.frame(patient = "Patient 5", severity = "moderate", progress = "recovering", outcome = "survived", pathway = fgseaRes_mc5$pathway, size = fgseaRes_mc5$size,
                    NES = fgseaRes_mc5$NES, padj = fgseaRes_mc5$padj) %>% filter(padj < 0.05) %>% arrange(desc(NES))


allc <- rbind(allc, first(all5c[all5c$NES > 0,], 5), last(all5c[all5c$NES < 0,], 5))
allm <- rbind(allm, first(all5m[all5m$NES > 0,], 5), last(all5m[all5m$NES < 0,], 5))

rm(list=setdiff(ls(), c("BCR", "fgsea_sets", "m_df", "allc", "allm")))

#For Patient 6:
pt6 <- BCR[,BCR$patient == "Patient 6"]
pt6$sample <- droplevels(pt6$sample)
DefaultAssay(pt6) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto: "severe122_Patient6", 

wlx.mrk.pt6.cm <- wilcoxauc(pt6, 'sample',
                         c("critical120_Patient6",  "moderate124_Patient6"),
                         seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(padj)
wlx.mrk.pt6.cs <- wilcoxauc(pt6, 'sample',
                         c("critical120_Patient6",  "severe122_Patient6"),
                         seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(padj)
wlx.mrk.pt6.sm <- wilcoxauc(pt6, 'sample',
                         c("severe122_Patient6",  "moderate124_Patient6"),
                         seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(padj)

head(wlx.mrk.pt6)

# we have all the genes for each cluster
dplyr::count(wlx.mrk.pt6, group)

# summarize the top abundant marker features for each group
top_markers(wlx.mrk.pt6)

write.table(top_markers(wlx.mrk.pt6), file = "top-wilcox-markers-pt6.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes.cs <- wlx.mrk.pt6.cs %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank

ranked.genes.cm <- wlx.mrk.pt6.cm %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank

ranked.genes.sm <- wlx.mrk.pt6.sm %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank


# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-pt6.tsv", 
            sep="\t", 
            append = FALSE, 
            quote=FALSE, 
            row.names = TRUE, 
            col.names = NA)

# choose top upregulated genes
topUpcrit.cs <- ranked.genes.cs %>% 
  dplyr::filter(group == "critical120_Patient6") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)
topUpcrit.cm <- ranked.genes.cm %>% 
  dplyr::filter(group == "critical120_Patient6") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)


topUpsev.cs <- ranked.genes.cs %>% 
  dplyr::filter(group == "severe122_Patient6") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)
topUpsev.sm <- ranked.genes.sm %>% 
  dplyr::filter(group == "severe122_Patient6") %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)


topUpmod.cm <- ranked.genes.cm %>% 
  dplyr::filter(group == "moderate124_Patient6") %>%
  filter(rank > 0) %>% 
  top_n(60, wt=-padj)
topUpmod.sm <- ranked.genes.sm %>% 
  dplyr::filter(group == "moderate124_Patient6") %>%
  filter(rank > 0) %>% 
  top_n(60, wt=-padj)


topPathways <- bind_rows(topUpcrit.cs, topUpcrit.cm, topUpsev.cs, topUpsev.sm, topUpmod.cm, topUpmod.sm)

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
enrich.topUpcrit6.cs <- enrichr(genes = topUpcrit.cs$feature, 
                             databases = "GO_Biological_Process_2021")
enrich.topUpcrit6.cm <- enrichr(genes = topUpcrit.cm$feature, 
                             databases = "GO_Biological_Process_2021")
enrich.topUpsev6.cs <- enrichr(genes = topUpsev.cs$feature, 
                            databases = "GO_Biological_Process_2021")
enrich.topUpsev6.sm <- enrichr(genes = topUpsev.sm$feature, 
                            databases = "GO_Biological_Process_2021")
enrich.topUpmod6.cm <- enrichr(genes = topUpmod.cm$feature, 
                            databases = "GO_Biological_Process_2021")
enrich.topUpmod6.sm <- enrichr(genes = topUpmod.sm$feature, 
                            databases = "GO_Biological_Process_2021")

enrich.topUpmod6.cm[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpmod6.cm[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)
enrich.topUpmod6.sm[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpmod6.sm[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)
enrich.topUpsev6.cs[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpsev6.cs[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)
enrich.topUpsev6.sm[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpsev6.sm[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)
enrich.topUpcrit6.cs[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpcrit6.cs[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)
enrich.topUpcrit6.cm[["GO_Biological_Process_2021"]]$Term <- vapply(str_split(enrich.topUpcrit6.cm[["GO_Biological_Process_2021"]]$Term, "[(GO:*)]"), "[", "", 1)


ggsave(filename = "graphs/enrichment-pt6-wlx-crit-sev.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit6.cs[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 6 critical vs severe"))
ggsave(filename = "graphs/enrichment-pt6-wlx-crit-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit6.cm[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 6 critical vs moderate"))

ggsave(filename = "graphs/enrichment-pt6-wlx-sev-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpsev6.cs[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 6 severe vs critical"))
ggsave(filename = "graphs/enrichment-pt6-wlx-sev-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpsev6.sm[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 6 severe vs moderate"))

ggsave(filename = "graphs/enrichment-pt6-wlx-mod-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod6.cm[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 6 moderate vs critical"))
ggsave(filename = "graphs/enrichment-pt6-wlx-mod-sev.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod6.sm[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 6 moderate vs severe"))

#Separating moderate genes from critical genes:
MNP.genes <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%
  arrange(-rank)  

#Selecting only the feature and rank columns of DGE data for fgsea run:
MNP.genes_cs <- ranked.genes.cs %>%
  dplyr::filter(group == "critical120_Patient6") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)
MNP.genes_cm <- ranked.genes.cm %>%
  dplyr::filter(group == "critical120_Patient6") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)

MNP.genes_sc <- ranked.genes.cs %>%
  dplyr::filter(group == "severe122_Patient6") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)
MNP.genes_sm <- ranked.genes.sm %>%
  dplyr::filter(group == "severe122_Patient6") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)

MNP.genes_mc <- ranked.genes.cm %>%
  dplyr::filter(group == "moderate124_Patient6") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)
MNP.genes_ms <- ranked.genes.sm %>%
  dplyr::filter(group == "moderate124_Patient6") %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)


#Converting matrices to tables:
ranks_cs <- deframe(MNP.genes_cs)
ranks_cm <- deframe(MNP.genes_cm)
ranks_sc <- deframe(MNP.genes_sc)
ranks_sm <- deframe(MNP.genes_sm)
ranks_mc <- deframe(MNP.genes_mc)
ranks_ms <- deframe(MNP.genes_ms)

#Running fgsea based on gene set ranking (stats = ranks):
fgseaRes_cs6 <- fgseaMultilevel(fgsea_sets, stats = ranks_cs, nPermSimple = 1000, 
                              maxSize = 200) %>% arrange(padj)
fgseaRes_cm6 <- fgseaMultilevel(fgsea_sets, stats = ranks_cm, nPermSimple = 1000, 
                              maxSize = 200) %>% arrange(padj)
fgseaRes_sc6 <- fgseaMultilevel(fgsea_sets, stats = ranks_sc, nPermSimple = 1000, 
                              maxSize = 200) %>% arrange(padj)
fgseaRes_sm6 <- fgseaMultilevel(fgsea_sets, stats = ranks_sm, nPermSimple = 1000, 
                              maxSize = 200) %>% arrange(padj)
fgseaRes_mc6 <- fgseaMultilevel(fgsea_sets, stats = ranks_mc, nPermSimple = 1000, 
                              maxSize = 200) %>% arrange(padj)
fgseaRes_ms6 <- fgseaMultilevel(fgsea_sets, stats = ranks_ms, nPermSimple = 1000, 
                              maxSize = 200) %>% arrange(padj)


#Enrichment plots:
ggsave(filename = "graphs/enrichment-plot-pt6-critical-severe.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_cs, desc(NES))[1]$pathway]],
                             ranks_cs) + 
         labs(title = arrange(fgseaRes_cs, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))
ggsave(filename = "graphs/enrichment-plot-pt6-critical-moderate.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_cm, desc(NES))[1]$pathway]],
                             ranks_cm) + 
         labs(title = arrange(fgseaRes_cm, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))

ggsave(filename = "graphs/enrichment-plot-pt6-severe-critical.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_sc, desc(NES))[1]$pathway]],
                             ranks_sc) + 
         labs(title = arrange(fgseaRes_sc, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))
ggsave(filename = "graphs/enrichment-plot-pt6-severe-moderate.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_sm, desc(NES))[1]$pathway]],
                             ranks_sm) + 
         labs(title = arrange(fgseaRes_sm, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))

ggsave(filename = "graphs/enrichment-plot-pt6-moderate-critical.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_mc, desc(NES))[1]$pathway]],
                             ranks_mc) + 
         labs(title = arrange(fgseaRes_mc, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))
ggsave(filename = "graphs/enrichment-plot-pt6-moderate-severe.tiff", 
       dpi = 300, width = 10, height = 5, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_ms, desc(NES))[1]$pathway]],
                             ranks_ms) + 
         labs(title = arrange(fgseaRes_ms, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))



#Plotting pathways:
fgseaRes_cs6$adjPvalue <- ifelse(fgseaRes_cs6$padj <= 0.05, "significant", "non-significant")
cols <- c("significant" = "red") #"non-significant" = "grey", 
plot_cs <- ggplot(fgseaRes_cs6, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = "Patient 6 - critical vs severe") + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='right', legend.key.size = unit(0.25, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) 
ggsave(filename = "graphs/gsea-pt6-critical-vs-severe-top.tiff",
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_cs)


fgseaRes_cm6$adjPvalue <- ifelse(fgseaRes_cm6$padj <= 0.05, "significant", "non-significant")
plot_cm <- ggplot(fgseaRes_cm6, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = "Patient 6 - critical vs moderate") + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='right', legend.key.size = unit(0.25, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) 
ggsave(filename = "graphs/gsea-pt6-critical-vs-moderate-top.tiff",
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_cm)

fgseaRes_sc6$adjPvalue <- ifelse(fgseaRes_sc6$padj <= 0.05, "significant", "non-significant")
plot_sc <- ggplot(fgseaRes_sc6, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = "Patient 6 - severe vs critical") + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='right', legend.key.size = unit(0.25, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) 
ggsave(filename = "graphs/gsea-pt6-severe-vs-critical-top.tiff",
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_sc)

fgseaRes_sm6$adjPvalue <- ifelse(fgseaRes_sm6$padj <= 0.05, "significant", "non-significant")
plot_sm <- ggplot(fgseaRes_sm6, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = "Patient 6 - severe vs moderate") + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='right', legend.key.size = unit(0.25, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) 
ggsave(filename = "graphs/gsea-pt6-severe-vs-moderate-top.tiff",
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_sm)



fgseaRes_mc6$adjPvalue <- ifelse(fgseaRes_mc6$padj <= 0.05, "significant", "non-significant")
plot_mc <- ggplot(fgseaRes_mc6, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = "Patient 6 - moderate vs critical") + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='right', legend.key.size = unit(0.25, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) 
ggsave(filename = "graphs/gsea-pt6-moderatee-vs-critical-top.tiff",
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_mc)

fgseaRes_ms6$adjPvalue <- ifelse(fgseaRes_ms6$padj <= 0.05, "significant", "non-significant")
plot_ms <- ggplot(fgseaRes_ms6, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = "Patient 6 - moderate vs severe") + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='right', legend.key.size = unit(0.25, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) 
ggsave(filename = "graphs/gsea-pt6-moderatee-vs-severe-top.tiff",
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_ms)





all6c <- data.frame(patient = "Patient 6", severity = "critical", progress = "recovering", outcome = "survived", pathway = fgseaRes_cm6$pathway, size = fgseaRes_cm6$size,
                    NES = fgseaRes_cm6$NES, padj = fgseaRes_cm6$padj) %>% arrange(desc(NES)) %>% filter(padj < 0.05)
all6m <- data.frame(patient = "Patient 6", severity = "moderate", progress = "recovering", outcome = "survived", pathway = fgseaRes_mc6$pathway, size = fgseaRes_mc6$size,
                    NES = fgseaRes_mc6$NES, padj = fgseaRes_mc6$padj) %>% arrange(desc(NES)) %>% filter(padj < 0.05)


allc <- rbind(allc, first(all6c[all6c$NES > 0,], 5), last(all6c[all6c$NES < 0,], 5))
allm <- rbind(allm, first(all6m[all6m$NES > 0,], 5), last(all6m[all6m$NES < 0,], 5))

write.table(x = allc, file = "allc.tsv", row.names = F, sep = "\t")
write.table(x = allm, file = "allm.tsv", row.names = F, sep = "\t")

rm(list=setdiff(ls(), c("BCR", "fgsea_sets", "m_df", "allc", "allm")))
gc()


#Plotting all:

ggsave(filename = "all-pathways-moderate.tiff",
       path = "graphs/",
       dpi = "print", width = 15, height = 10,
       plot = ggplot(allm,
                     aes(x = NES, y = reorder(pathway,NES), color = progress, size = size)) +
         geom_point(alpha=0.5) +
         #facet_wrap(~patient, ncol = 3) +
         #facet_grid(vars(severity), vars(TRBV)) +
         #facet_grid(~TRBV + severity) +
         #labs(title = "Enriched Pathways", vjust = 1) +
         theme(plot.subtitle=element_text(size=12, face="italic", color="black")) +
         theme(axis.text=element_text(size=8),
               axis.title=element_text(size=10)) +
         ylab("pathways") +
         geom_vline(xintercept = 0, linetype="dotted") +
         scale_size(range = c(3, 8), name="genes") +
         #scale_fill_viridis(discrete=TRUE, guide="none", option="A") +
         theme_bw()
)

ggsave(filename = "all-pathways-critical.tiff",
       path = "graphs/",
       dpi = "print", width = 15, height = 10,
       plot = ggplot(allc,
                     aes(x = NES, y = reorder(pathway,NES), color = progress, size = size)) +
         geom_point(alpha=0.5) +
         #facet_wrap(~patient, ncol = 3) +
         #facet_grid(vars(severity), vars(TRBV)) +
         #facet_grid(~TRBV + severity) +
         #labs(title = "Enriched Pathways", vjust = 1) +
         theme(plot.subtitle=element_text(size=12, face="italic", color="black")) +
         theme(axis.text=element_text(size=8),
               axis.title=element_text(size=10)) +
         ylab("pathways") +
         geom_vline(xintercept = 0, linetype="dotted") +
         scale_size(range = c(3, 8), name="genes") +
         #scale_fill_viridis(discrete=TRUE, guide="none", option="A") +
         theme_bw()
)

#Adding a column for shared/unique pathways:
allm$shared <- allm$progress
allm[allm$pathway %in% intersect(allm[allm$progress == "progressing",]$pathway, allm[allm$progress == "recovering",]$pathway),]$shared <- "shared"

allc$shared <- allc$progress
allc[allc$pathway %in% intersect(allc[allc$progress == "progressing",]$pathway, allc[allc$progress == "recovering",]$pathway),]$shared <- "shared"

allm$shared <- factor(allm$shared, levels = c("shared","progressing", "recovering"))
allc$shared <- factor(allc$shared, levels = c("shared","progressing", "recovering"))

#Setting facet labels
new_labels <- c("shared" = "Common Pathways", "progressing" = "Pathways in Progressing Patients", "recovering" = "Pathways in Recovering Patients")

#Plotting bubble plots, divided by progress, including the shared pathways:
shared <- ggplot(transform(allm[allm$shared == "shared",]),
                 aes(x = NES, y = reorder(pathway,NES), color = progress, size = size)) +
  geom_point(alpha=0.7) +
  facet_grid(vars(shared), labeller = labeller(shared = new_labels), drop = TRUE, space = "free") +
  ylab("Pathways") +
  geom_vline(xintercept = 0, linetype="dotted") +
  scale_size(range = c(3, 8), name="genes") +
  scale_x_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3,4), limits = c(-4,4)) +
  theme_bw() +
  theme(plot.subtitle=element_text(size=12, face="italic", color="black")) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

progressing <- ggplot(transform(allm[allm$shared == "progressing",]),
                      aes(x = NES, y = reorder(pathway,NES), size = size)) +
  geom_point(alpha=0.7, color='#F8766D') +  
  facet_grid(vars(shared), labeller = labeller(shared = new_labels), drop = TRUE, space = "free") +
  theme(strip.text.x = element_blank()) +
  ylab("Pathways") +
  geom_vline(xintercept = 0, linetype="dotted") +
  scale_size(range = c(3, 8), name="genes") +
  scale_x_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3,4), limits = c(-4,4)) +
  theme_bw() + 
  theme(axis.text=element_text(size=12), 
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.title=element_text(size=14),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  

recovering <- ggplot(transform(allm[allm$shared == "recovering",]),
                      aes(x = NES, y = reorder(pathway,NES), size = size)) +
  geom_point(alpha=0.7, color='#00BFC4') +  
  facet_grid(vars(shared), labeller = labeller(shared = new_labels), drop = TRUE, space = "free") +
  theme(strip.text.x = element_blank()) +
  ylab("Pathways") +
  geom_vline(xintercept = 0, linetype="dotted") +
  scale_size(range = c(3, 8), name="genes") +
  scale_x_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3,4), limits = c(-4,4)) +
  theme_bw() +
  theme(axis.text=element_text(size=12), 
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.title=element_text(size=14))

bubble <- ggarrange(shared + rremove("xlab") + rremove("ylab"), progressing + rremove("xlab") + rremove("ylab"), recovering + rremove("ylab"),
                    ncol=1, nrow=3, common.legend = TRUE, legend="right",
                    labels = NULL, heights = c(0.5, 0.25, 0.25),
                    align = "v", 
                    font.label = list(size = 14, color = "black", face = "bold", family = NULL, position = "top"))

bubble <- annotate_figure(bubble, left = textGrob("Pathways", rot = 90, vjust = 1, gp = gpar(fontsize = 14)))

tiff(filename = "graphs/pathways-bubble.tiff",
       width = 20, height = 12, units = "in", res = 300)
bubble
dev.off()


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
ISGs <- ISGs$ISGs_geneset_227

#Cleaning gene set:
ISGs <- ISGs[!grepl(ISGs, pattern = "LOC") & !grepl(ISGs, pattern = "XENTR") & !grepl(ISGs, pattern = "GSON")]
ISGs <- ISGs[ISGs != ""]


obj.b <- BCR[,(BCR$patient != "Patient 5" & BCR$patient != "Patient 6") & BCR$severity != "healthy"]
obj.b$sample <- droplevels(obj.b$sample)
obj.b$zeverity <- "critical"
obj.b$zeverity[obj.b$severity == "mild" | obj.b$severity == "moderate"] <- "moderate"
DefaultAssay(obj.b) <- "RNA"

#For B intermediate:
obj.b.int <- obj.b[,obj.b$azimuthNames == "B intermediate"]

#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj.b.int, group_by = "zeverity", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(desc(abs(logFC)))

# rank genes using rrho algorithm
MNP.genes_B.int <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                # sort genes by decreasing rank
#head(n = 10)

# choose top upregulated genes
topUp_B_mod.int <- MNP.genes_B.int %>% 
  dplyr::filter(group == c('moderate')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>% 
  top_n(50, wt=-padj)

topUp_B_crit.int <- MNP.genes_B.int %>% 
  dplyr::filter(group == c('critical')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>%
  top_n(50, wt=-padj)


## moderate vs critical
topPathways_all.int <- bind_rows(topUp_B_mod.int, topUp_B_crit.int)


#Subsetting to only DEGs that are DEGs:
obj.b.int <- obj.b.int[ISGs[ISGs %in% topPathways_all$feature],]
obj.b.int <- NormalizeData(obj.b.int)

#Changing the data to "pseudo-bulk", using ISG set:
obj.b.tmp <- ScaleData(obj.b.int)
avg.exp.mat.b.int <- AverageExpression(obj.b.tmp, group.by = "sample", slot = 'scale.data')
avg.exp.rna.b.int <- avg.exp.mat.b.int$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.b.tmp@meta.data %>% 
  group_by(sample, zeverity) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data.B.int <- data.frame(t(avg.exp.rna.b.int), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data.B.int <- meta.data.B.int[,c(colnames(meta.data.B.int)[colnames(meta.data.B.int) %in% rownames(obj.b.tmp@assays$RNA@scale.data)], "sample", "zeverity")]


#Preparing aesthetic parameters:
expr.cols <- colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow"))
col.anno <- HeatmapAnnotation(Severity = meta.data.B.int$zeverity, show_annotation_name = F,   
                              col = list(Severity = c("critical" = "red", "moderate" = "blue")))
column_labels <- c("moderate272_Patient1" = "Pt1", "critical293_Patient1" = "Pt1",
                   "moderate303_Patient2" = "Pt2", "critical308_Patient2" = "Pt2",
                   "mild186_Patient3" = "Pt3", "critical213_Patient3" = "Pt3",
                   "mild227_Patient4" = "Pt4", "critical238_Patient4" = "Pt4")

#Preparing separate matrices for each cell type:
rownames(meta.data.B.int) <- meta.data.B.int$sample
meta.data.B.int$zeverity <- factor(meta.data.B.int$zeverity, levels = c("moderate", "critical"))
hm.b.int <- Heatmap(t(meta.data.B.int[,colnames(meta.data.B.int) %in% diff.genes$feature]), 
                column_split = meta.data.B.int$zeverity, name = "Expression", 
                column_labels = column_labels, row_title = "B intermediate", column_title = " ",  
                col = expr.cols, top_annotation = col.anno, row_names_gp = gpar(fontsize = 7), 
                show_row_dend = F, show_column_dend = F, cluster_columns = F)


#For B memory:
obj.b.mem <- obj.b[,obj.b$azimuthNames == "B memory"]

#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj.b.mem, group_by = "zeverity", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(desc(abs(logFC)))

# rank genes using rrho algorithm
MNP.genes_B.mem <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                # sort genes by decreasing rank
#head(n = 10)

# choose top upregulated genes
topUp_B_mod.mem <- MNP.genes_B.mem %>% 
  dplyr::filter(group == c('moderate')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>% 
  top_n(50, wt=-padj)

topUp_B_crit.mem <- MNP.genes_B.mem %>% 
  dplyr::filter(group == c('critical')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>%
  top_n(50, wt=-padj)


## moderate vs critical
topPathways_all.mem <- bind_rows(topUp_B_mod.mem, topUp_B_crit.mem)


#Subsetting to only DEGs that are DEGs:
obj.b.mem <- obj.b.mem[ISGs[ISGs %in% topPathways_all$feature],]
obj.b.mem  <- NormalizeData(obj.b.mem)

#Changing the data to "pseudo-bulk", using ISG set:
obj.b.tmp <- ScaleData(obj.b.mem)
avg.exp.mat.b.mem <- AverageExpression(obj.b.tmp, group.by = "sample", slot = 'scale.data')
avg.exp.rna.b.mem <- avg.exp.mat.b.mem$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.b.tmp@meta.data %>% 
  group_by(sample, zeverity) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data.B.mem <- data.frame(t(avg.exp.rna.b.mem), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data.B.mem <- meta.data.B.mem[,c(colnames(meta.data.B.mem)[colnames(meta.data.B.mem) %in% rownames(obj.b.tmp@assays$RNA@scale.data)], "sample", "zeverity")]



rownames(meta.data.B.mem) <- meta.data.B.mem$sample
hm.b.mem <- Heatmap(t(meta.data.B.mem[,colnames(meta.data.B.mem) %in% diff.genes$feature]),
                    column_split = meta.data.B.mem$zeverity,
                    column_labels = column_labels, row_title = "B memory",
                    col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 7),
                    show_row_dend = F, show_column_dend = F, cluster_columns = F)


#####For B naive:

obj.b.nai <- obj.b[,obj.b$azimuthNames == "B naive"]

#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj.b.nai, group_by = "zeverity", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(desc(abs(logFC)))

# rank genes using rrho algorithm
MNP.genes_B.nai <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                # sort genes by decreasing rank
#head(n = 10)

# choose top upregulated genes
topUp_B_mod.nai <- MNP.genes_B.nai %>% 
  dplyr::filter(group == c('moderate')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>% 
  top_n(50, wt=-padj)

topUp_B_crit.nai <- MNP.genes_B.nai %>% 
  dplyr::filter(group == c('critical')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>%
  top_n(50, wt=-padj)


## moderate vs critical
topPathways_all.nai <- bind_rows(topUp_B_mod.nai, topUp_B_crit.nai)


#Subsetting to only DEGs that are DEGs:
obj.b.nai <- obj.b.nai[ISGs[ISGs %in% topPathways_all$feature],]
obj.b.nai  <- NormalizeData(obj.b.nai)

#Changing the data to "pseudo-bulk", using ISG set:
obj.b.tmp <- ScaleData(obj.b.nai)
avg.exp.mat.b.nai <- AverageExpression(obj.b.tmp, group.by = "sample", slot = 'scale.data')
avg.exp.rna.b.nai <- avg.exp.mat.b.nai$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.b.tmp@meta.data %>% 
  group_by(sample, zeverity) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data.B.nai <- data.frame(t(avg.exp.rna.b.nai), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data.B.nai <- meta.data.B.nai[,c(colnames(meta.data.B.nai)[colnames(meta.data.B.nai) %in% rownames(obj.b.tmp@assays$RNA@scale.data)], "sample", "zeverity")]



rownames(meta.data.B.nai) <- meta.data.B.nai$sample
hm.b.nai <- Heatmap(t(meta.data.B.nai[,colnames(meta.data.B.mem) %in% diff.genes$feature]), 
                    column_split = meta.data.B.nai$zeverity,
                    column_labels = column_labels, row_title = "B naive",
                    col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 7),
                    show_row_dend = F, show_column_dend = F, cluster_columns = F)

#####For plasmablasts:

obj.b.pla <- obj.b[,obj.b$azimuthNames == "Plasmablast"]

#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj.b.pla, group_by = "zeverity", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(desc(abs(logFC)))

# rank genes using rrho algorithm
MNP.genes_B.pla <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                # sort genes by decreasing rank
#head(n = 10)

# choose top upregulated genes
topUp_B_mod.pla <- MNP.genes_B.pla %>% 
  dplyr::filter(group == c('moderate')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>% 
  top_n(50, wt=-padj)

topUp_B_crit.pla <- MNP.genes_B.pla %>% 
  dplyr::filter(group == c('critical')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>%
  top_n(50, wt=-padj)


## moderate vs critical
topPathways_all.pla <- bind_rows(topUp_B_mod.pla, topUp_B_crit.pla)


#Subsetting to only DEGs that are DEGs:
obj.b.pla <- obj.b.pla[ISGs[ISGs %in% topPathways_all$feature],]
obj.b.pla <- NormalizeData(obj.b.pla)

#Changing the data to "pseudo-bulk", using ISG set:
obj.b.tmp <- ScaleData(obj.b.pla)
avg.exp.mat.b.pla <- AverageExpression(obj.b.tmp, group.by = "sample", slot = 'scale.data')
avg.exp.rna.b.pla <- avg.exp.mat.b.pla$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.b.tmp@meta.data %>% 
  group_by(sample, zeverity) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data.B.pla <- data.frame(t(avg.exp.rna.b.pla), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data.B.pla <- meta.data.B.pla[,c(colnames(meta.data.B.pla)[colnames(meta.data.B.pla) %in% rownames(obj.b.tmp@assays$RNA@scale.data)], "sample", "zeverity")]


rownames(meta.data.B.pla) <- meta.data.B.pla$sample
hm.plasm <- Heatmap(t(meta.data.B.pla[,colnames(meta.data.B.mem) %in% diff.genes$feature]), column_split = meta.data.B.pla$zeverity,
                    column_labels = column_labels, row_title = "Plasmablasts",
                    col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 7),
                    show_row_dend = F, show_column_dend = F, cluster_columns = F)

#Adding all heatmaps to one, and saving it:
hm.list <- hm.b.int %v% hm.b.mem %v% hm.b.nai %v% hm.plasm

tiff("./graphs/hm-critical-moderate-b-cells-isg-in-deg-pt1-4.tiff",
     res = 300, width = 4, height = 8, units = "in")
draw(hm.list)
dev.off()


#####For Patients 5 and 6


obj.b <- BCR[,(BCR$patient == "Patient 5" | BCR$patient == "Patient 6") & BCR$severity != "severe"]
obj.b$sample <- droplevels(obj.b$sample)
DefaultAssay(obj.b) <- "RNA"



#For B intermediate:
obj.b.int <- obj.b[,obj.b$azimuthNames == "B intermediate"]

#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj.b.int, group_by = "severity", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(desc(abs(logFC)))

# rank genes using rrho algorithm
MNP.genes_B.int <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                # sort genes by decreasing rank

# choose top upregulated genes
topUp_B_mod.int <- MNP.genes_B.int %>% 
  dplyr::filter(group == c('moderate')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>% 
  top_n(50, wt=-padj)

topUp_B_crit.int <- MNP.genes_B.int %>% 
  dplyr::filter(group == c('critical')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>%
  top_n(50, wt=-padj)


## moderate vs critical
topPathways_all.int <- bind_rows(topUp_B_mod.int, topUp_B_crit.int)


#Subsetting to only DEGs that are DEGs:
obj.b.int <- obj.b.int[ISGs[ISGs %in% topPathways_all$feature],]
obj.b.int <- NormalizeData(obj.b.int)

#Changing the data to "pseudo-bulk", using ISG set:
obj.b.tmp <- ScaleData(obj.b.int)
avg.exp.mat.b.int <- AverageExpression(obj.b.tmp, group.by = "sample", slot = 'scale.data')
avg.exp.rna.b.int <- avg.exp.mat.b.int$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.b.tmp@meta.data %>% 
  group_by(sample, severity) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data.B.int <- data.frame(t(avg.exp.rna.b.int), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data.B.int <- meta.data.B.int[,c(colnames(meta.data.B.int)[colnames(meta.data.B.int) %in% rownames(obj.b.tmp@assays$RNA@scale.data)], "sample", "severity")]


#Preparing aesthetic parameters:
expr.cols <- colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow"))
col.anno <- HeatmapAnnotation(Severity = factor(meta.data.B.int$severity), show_annotation_name = F,   
                              col = list(Severity = c("critical" = "red", "moderate" = "blue")))
column_labels <- c("critical119_Patient5" = "Pt5", "moderate138_Patient5" = "Pt5",
                   "critical120_Patient6" = "Pt6", "moderate124_Patient6" = "Pt6")

#Preparing separate matrices for each cell type:
rownames(meta.data.B.int) <- meta.data.B.int$sample
meta.data.B.int$severity <- factor(meta.data.B.int$severity, levels = c("moderate", "critical"))
hm.b.int <- Heatmap(t(meta.data.B.int[,colnames(meta.data.B.int) %in% diff.genes$feature]), 
                    column_split = meta.data.B.int$severity, name = "Expression", 
                    column_labels = column_labels, row_title = "B intermediate", column_title = " ",  
                    col = expr.cols, top_annotation = col.anno, row_names_gp = gpar(fontsize = 7), 
                    show_row_dend = F, show_column_dend = F, cluster_columns = F)


#For B memory:
obj.b.mem <- obj.b[,obj.b$azimuthNames == "B memory"]

#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj.b.mem, group_by = "severity", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(desc(abs(logFC)))

# rank genes using rrho algorithm
MNP.genes_B.mem <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                # sort genes by decreasing rank

# choose top upregulated genes
topUp_B_mod.mem <- MNP.genes_B.mem %>% 
  dplyr::filter(group == c('moderate')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>% 
  top_n(50, wt=-padj)

topUp_B_crit.mem <- MNP.genes_B.mem %>% 
  dplyr::filter(group == c('critical')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>%
  top_n(50, wt=-padj)


## moderate vs critical
topPathways_all.mem <- bind_rows(topUp_B_mod.mem, topUp_B_crit.mem)


#Subsetting to only DEGs that are DEGs:
obj.b.mem <- obj.b.mem[ISGs[ISGs %in% topPathways_all.mem$feature],]
obj.b.mem  <- NormalizeData(obj.b.mem)

#Changing the data to "pseudo-bulk", using ISG set:
obj.b.tmp <- ScaleData(obj.b.mem, features = ISGs[ISGs %in% topPathways_all.mem$feature])
avg.exp.mat.b.mem <- AverageExpression(obj.b.tmp, group.by = "sample", slot = 'scale.data')
avg.exp.rna.b.mem <- avg.exp.mat.b.mem$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.b.tmp@meta.data %>% 
  group_by(sample, severity) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data.B.mem <- data.frame(t(avg.exp.rna.b.mem[rownames(avg.exp.rna.b.mem) %in% diff.genes,]), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data.B.mem <- meta.data.B.mem[,c(colnames(meta.data.B.mem)[colnames(meta.data.B.mem) %in% rownames(obj.b.tmp@assays$RNA@scale.data)], "sample", "severity")]



rownames(meta.data.B.mem) <- meta.data.B.mem$sample
hm.b.mem <- Heatmap(t(meta.data.B.mem[,colnames(meta.data.B.mem) %in% diff.genes$feature]),
                    column_split = meta.data.B.mem$severity,
                    column_labels = column_labels, row_title = "B memory",
                    col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 7),
                    show_row_dend = F, show_column_dend = F, cluster_columns = F)


#####For B naive:

obj.b.nai <- obj.b[,obj.b$azimuthNames == "B naive"]

#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj.b.nai, group_by = "severity", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(desc(abs(logFC)))

# rank genes using rrho algorithm
MNP.genes_B.nai <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                # sort genes by decreasing rank
#head(n = 10)

# choose top upregulated genes
topUp_B_mod.nai <- MNP.genes_B.nai %>% 
  dplyr::filter(group == c('moderate')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>% 
  top_n(50, wt=-padj)

topUp_B_crit.nai <- MNP.genes_B.nai %>% 
  dplyr::filter(group == c('critical')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>%
  top_n(50, wt=-padj)


## moderate vs critical
topPathways_all.nai <- bind_rows(topUp_B_mod.nai, topUp_B_crit.nai)


#Subsetting to only DEGs that are DEGs:
obj.b.nai <- obj.b.nai[ISGs[ISGs %in% topPathways_all$feature],]
obj.b.nai  <- NormalizeData(obj.b.nai)

#Changing the data to "pseudo-bulk", using ISG set:
obj.b.tmp <- ScaleData(obj.b.nai)
avg.exp.mat.b.nai <- AverageExpression(obj.b.tmp, group.by = "sample", slot = 'scale.data')
avg.exp.rna.b.nai <- avg.exp.mat.b.nai$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.b.tmp@meta.data %>% 
  group_by(sample, severity) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data.B.nai <- data.frame(t(avg.exp.rna.b.nai[rownames(avg.exp.rna.b.nai) %in% diff.genes,]), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data.B.nai <- meta.data.B.nai[,c(colnames(meta.data.B.nai)[colnames(meta.data.B.nai) %in% rownames(obj.b.tmp@assays$RNA@scale.data)], "sample", "severity")]


#STAT1
rownames(meta.data.B.nai) <- meta.data.B.nai$sample
hm.b.nai <- Heatmap(t(meta.data.B.nai[,colnames(meta.data.B.nai) %in% diff.genes$feature]), 
                    column_split = meta.data.B.nai$severity, 
                    column_labels = column_labels, row_title = "B naive",
                    col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 7),
                    show_row_dend = F, show_column_dend = F, cluster_columns = F)

#####For plasmablasts:

obj.b.pla <- obj.b[,obj.b$azimuthNames == "Plasmablast"]

#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(obj.b.pla, group_by = "severity", seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(desc(abs(logFC)))

# rank genes using rrho algorithm
MNP.genes_B.pla <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                # sort genes by decreasing rank

# choose top upregulated genes
topUp_B_mod.pla <- MNP.genes_B.pla %>% 
  dplyr::filter(group == c('moderate')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>% 
  top_n(50, wt=-padj)

topUp_B_crit.pla <- MNP.genes_B.pla %>% 
  dplyr::filter(group == c('critical')) %>%
  filter(rank > 0) %>% 
  subset(padj < 0.05) %>%
  top_n(50, wt=-padj)


## moderate vs critical
topPathways_all.pla <- bind_rows(topUp_B_mod.pla, topUp_B_crit.pla)


#Subsetting to only DEGs that are DEGs:
obj.b.pla <- obj.b.pla[ISGs[ISGs %in% topPathways_all$feature],]
obj.b.pla <- NormalizeData(obj.b.pla)

#Changing the data to "pseudo-bulk", using ISG set:
obj.b.tmp <- ScaleData(obj.b.pla)
avg.exp.mat.b.pla <- AverageExpression(obj.b.tmp, group.by = "sample", slot = 'scale.data')
avg.exp.rna.b.pla <- avg.exp.mat.b.pla$RNA

#Extracting some meta data for "clustering" the heatmap:
obj.b.tmp@meta.data %>% 
  group_by(sample, severity) %>% count() -> thing

#Adding average/bulk expression data to meta data:
meta.data.B.pla <- data.frame(t(avg.exp.rna.b.pla), thing)

#Keeping only genes that exist in the scaled data in the Seurat object:
meta.data.B.pla <- meta.data.B.pla[,c(colnames(meta.data.B.pla)[colnames(meta.data.B.pla) %in% rownames(obj.b.tmp@assays$RNA@scale.data)], "sample", "severity")]


rownames(meta.data.B.pla) <- meta.data.B.pla$sample
hm.plasm <- Heatmap(t(meta.data.B.pla[,colnames(meta.data.B.pla) %in% diff.genes$feature]), column_split = meta.data.B.pla$severity,
                    column_labels = column_labels, row_title = "Plasmablasts",
                    col = expr.cols, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 7),
                    show_row_dend = F, show_column_dend = F, cluster_columns = F)





#Adding all heatmaps to one, and saving it: hm.b.mem %v% hm.b.nai %v%
hm.list <- hm.b.int %v% hm.plasm

tiff("./graphs/hm-critical-moderate-b-cells-isg-in-deg-pt5-6.tiff",
     res = 300, width = 4, height = 8, units = "in")
draw(hm.list)
dev.off()





################################Outcome Comparison#############################

#Comparing each of the first 3 to the last 3:

#Subsetting to moderate/mild and critical patients:
mod <- BCR[,BCR$severity == "moderate" | BCR$severity == "mild"]
mod$sample <- droplevels(mod$sample)
crit <- BCR[,BCR$severity == "critical"]
crit$sample <- droplevels(crit$sample)

for (s in c(mod, crit)) {
  for (i in c("Patient 1", "Patient 2", "Patient 3")) {
    for (j in c("Patient 4", "Patient 5", "Patient 6")) {
    #Wilcoxon test:
    wlx.mrk <- wilcoxauc(s, 'patient', c(i, j),
                         seurat_assay='RNA', assay = "data") %>%
      subset(padj < 0.05) %>% arrange(padj)
    
    #Ranking genes using rrho algorithm
    ranked.genes <- wlx.mrk %>%
      mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
      arrange(-rank)      
    
    
    #Selecting only the feature and rank columns of DGE data for fgsea run:
    MNP.genes.ded <- ranked.genes %>%
      dplyr::filter(group == i) %>%
      arrange(-rank) %>% 
      dplyr::select(feature, rank)
    
    MNP.genes.srv <- ranked.genes %>%
      dplyr::filter(group == j) %>%
      arrange(-rank) %>% 
      dplyr::select(feature, rank)
    
    
    #Converting matrices to tables:
    ranks_d <- deframe(MNP.genes.ded)
    ranks_s <- deframe(MNP.genes.srv)
    
    #Running fgsea based on gene set ranking (stats = ranks):
    fgseaRes_d <- fgseaMultilevel(fgsea_sets, stats = ranks_d, nPermSimple = 1000, 
                                   maxSize = 200) %>% arrange(padj)
    fgseaRes_s <- fgseaMultilevel(fgsea_sets, stats = ranks_s, nPermSimple = 1000, 
                                   maxSize = 200) %>% arrange(padj)
    
    
    
    #Plotting pathways:
    fgseaRes_d$adjPvalue <- ifelse(fgseaRes_d$padj <= 0.05, "significant", "non-significant")
    cols <- c("significant" = "red") #"non-significant" = "grey", 
    plot_d <- ggplot(fgseaRes_d, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
      geom_col() +
      scale_fill_manual(values = cols) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title = paste0(i, " vs ", j, " - ", ifelse(levels(factor(s$severity))[1] == "critical", "critical", "moderate"))) + 
      theme(axis.text.x = element_text(vjust = 0.5, size = 10),
            axis.text.y = element_text(size = 10, vjust = 0.5),
            plot.title = element_text(hjust = 0.5, size = 14),
            axis.title = element_text(size = 10),
            legend.position='right', legend.key.size = unit(0.25, 'cm'),
            legend.key.height = unit(0.25, 'cm'),
            legend.key.width = unit(0.25, 'cm'), 
            legend.title = element_text(size=10),
            legend.text = element_text(size=10)) 
    ggsave(path =  "graphs/", filename = paste0(i, "-vs-", j, "-pathways-", ifelse(levels(factor(s$severity))[1] == "critical", "critical", "moderate"), ".tiff"),
           width = 10, height = 10, dpi = 300, units = "in", plot = plot_d)
    
    #Saving a table of the results:
    write.table(x = fgseaRes_d %>% mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ", ", character(1L))),
                file = paste0(i, "-vs-", j, "-pathways-", ifelse(levels(factor(s$severity))[1] == "critical", "critical", "moderate"), ".tsv"),
                sep = "\t", col.names = NA)
   
    
  }

 }
}




