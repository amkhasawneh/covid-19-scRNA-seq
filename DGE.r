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

#Plasmablasts:healthy2_control2 has 2 cells
plasma <- BCR[,BCR$azimuthNames == "Plasmablast"]
plasma$sample <- droplevels(plasma$sample)
Idents(plasma) <- "sample"
plasma <- ScaleData(plasma)
mrk <- FindAllMarkers(plasma, min.pct = 0.25, logfc.threshold = 0.58) %>%
  subset(p_val_adj < 0.05) %>%
  arrange(p_val_adj)
ggsave(filename = "graphs/DEGs-plasmablasts.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(plasma, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                    genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-plasmablasts.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
                         title = "Plasmablasts",
                         y = "Count", orderBy = "P.value"))

#B intermediate: 
bint <- BCR[,BCR$azimuthNames == "B intermediate"]
bint$sample <- droplevels(bint$sample)
Idents(bint) <- "sample"
bint <- ScaleData(bint)
mrk <- FindAllMarkers(bint, min.pct = 0.25, logfc.threshold = 0.58) %>%
  subset(p_val_adj < 0.05) %>%
  arrange(p_val_adj)
ggsave(filename = "graphs/DEGs-B intermediate.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(plasma, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                    genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-B intermediate.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
                         title = "B intermediate",
                         y = "Count", orderBy = "P.value"))


Idents(BCR) <- "sample"
#Patient 1:
pt1 <- BCR[,BCR$patient == "Patient 1"]
pt1$sample <- droplevels(pt1$sample)
pt1 <- ScaleData(pt1)
mrk <- FindMarkers(pt1[,pt1$v_gene == "IGHV1-18"], ident.1 = "critical293_Patient1",
                   min.pct = 0.25, logfc.threshold = 0.58) 
ggsave(filename = "graphs/DEGs-pt1-IGHV1-18-all.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(pt1, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                    genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-pt1-IGHV1-18-all.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
                         title = "Patient 1 - IGHV1-18",
                         y = "Count", orderBy = "P.value"))

mrk <- FindMarkers(pt1[,pt1$v_gene == "IGHV1-18" & pt1$azimuthNames == "B intermediate"], ident.1 = "critical293_Patient1",
                   min.pct = 0.25, logfc.threshold = 0.58) 
ggsave(filename = "graphs/DEGs-pt1-IGHV1-18-bint.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(pt1, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                       genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-pt1-IGHV1-18-bint.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
                         title = "Patient 1 - IGHV1-18 - B intermediate",
                         y = "Count", orderBy = "P.value"))

#Patient 2:
pt2 <- BCR[,BCR$patient == "Patient 2"]
pt2$sample <- droplevels(pt2$sample)
pt2 <- ScaleData(pt2)
mrk <- FindMarkers(pt2[,pt2$v_gene == "IGHV4-59"], ident.1 = "critical308_Patient2",
                   min.pct = 0.25, logfc.threshold = 0.58) 
ggsave(filename = "graphs/DEGs-pt2-IGHV4-59-all.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(pt2, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                    genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-pt2-IGHV4-59-all.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
                         title = "Patient 2 - IGHV4-59",
                         y = "Count", orderBy = "P.value"))

mrk <- FindMarkers(pt2[,pt2$v_gene == "IGHV4-59" & pt2$azimuthNames == "Plasmablast"], ident.1 = "critical308_Patient2",
                   min.pct = 0.25, logfc.threshold = 0.58) 
ggsave(filename = "graphs/DEGs-pt2-IGHV4-59-plasmablasts.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(pt2, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                    genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-pt2-IGHV4-59-plasmablasts.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
                         title = "Patient 2 - IGHV4-59 - Plasmablasts",
                         y = "Count", orderBy = "P.value"))

mrk <- FindMarkers(pt2[,pt2$v_gene == "IGHV4-34"], ident.1 = "critical308_Patient2",
                   min.pct = 0.25, logfc.threshold = 0.58) 
ggsave(filename = "graphs/DEGs-pt2-IGHV4-34-all.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(pt2, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                    genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-pt2-IGHV4-34-all.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
                         title = "Patient 2 - IGHV4-34",
                         y = "Count", orderBy = "P.value"))

mrk <- FindMarkers(pt2[,pt2$v_gene == "IGHV4-34" & pt2$azimuthNames == "Plasmablast"], ident.1 = "critical308_Patient2",
                   min.pct = 0.25, logfc.threshold = 0.58) 
ggsave(filename = "graphs/DEGs-pt2-IGHV4-34-plasmablasts.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(pt2, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                    genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-pt2-IGHV4-34-plasmablasts.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
                         title = "Patient 2 - IGHV4-34 - Plasmablasts",
                         y = "Count", orderBy = "P.value"))
#Patient 3:
pt3 <- BCR[,BCR$patient == "Patient 3"]
pt3$sample <- droplevels(pt3$sample)
pt3 <- ScaleData(pt3)
mrk <- FindMarkers(pt3[,pt3$v_gene == "IGHV4-39"], ident.1 = "critical213_Patient3",
                   min.pct = 0.25, logfc.threshold = 0.58) %>% subset(p_val_adj < 0.05) %>% arrange(p_val_adj)
#Cell group 2 has fewer than 3 cells

mrk <- FindMarkers(pt3[,pt3$v_gene == "IGHV3-48"], ident.1 = "critical213_Patient3",
                   min.pct = 0.25, logfc.threshold = 0.58) %>% subset(p_val_adj < 0.05) %>% arrange(p_val_adj)
ggsave(filename = "graphs/DEGs-pt3-IGHV3-48-all.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(pt3, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                    genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-pt3-IGHV3-48-all.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
                         title = "Patient 3 - IGHV3-48",
                         y = "Count", orderBy = "P.value"))




#Patient 4:
pt4 <- BCR[,BCR$patient == "Patient 4"]
pt4$sample <- droplevels(pt4$sample)
pt4 <- ScaleData(pt4)
mrk <- FindMarkers(pt4[,pt4$v_gene == "IGHV3-33"], ident.1 = "critical238_Patient4",
                   min.pct = 0.25, logfc.threshold = 0.58)
ggsave(filename = "graphs/DEGs-pt4-IGHV3-33-all.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(pt4, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                    genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-pt4-IGHV3-33-all.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
                         title = "Patient 4 - IGHV3-33",
                         y = "Count", orderBy = "P.value"))


mrk <- FindMarkers(pt4[,pt4$v_gene == "IGHV3-23"], ident.1 = "critical238_Patient4",
                   min.pct = 0.25, logfc.threshold = 0.58)
ggsave(filename = "graphs/DEGs-pt4-IGHV3-23-all.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(pt4, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                    genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-pt4-IGHV3-23-all.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
                         title = "Patient 4 - IGHV3-23",
                         y = "Count", orderBy = "P.value"))




#Patient 5:
pt5 <- BCR[,BCR$patient == "Patient 5"]
pt5$sample <- droplevels(pt5$sample)
pt5 <- ScaleData(pt5)
mrk <- FindMarkers(pt5[,pt5$v_gene == "IGHV3-23"], ident.1 = "critical119_Patient5",
                   min.pct = 0.25, logfc.threshold = 0.58)
ggsave(filename = "graphs/DEGs-pt5-IGHV3-23-all.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(pt5, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                    genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-pt5-IGHV3-23-all.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
                         title = "Patient 5 - IGHV3-23",
                         y = "Count", orderBy = "P.value"))



#Patient 6:
pt6 <- BCR[,BCR$patient == "Patient 6"]
pt6$sample <- droplevels(pt6$sample)
pt6 <- ScaleData(pt6)
mrk <- FindMarkers(pt6[,pt6$v_gene == "IGHV3-23"], ident.1 = "critical120_Patient6",
                   min.pct = 0.25, logfc.threshold = 0.58)
ggsave(filename = "graphs/DEGs-pt6-IGHV3-23-all.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(pt6, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                    genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-pt6-IGHV3-23-all.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
                         title = "Patient 6 - IGHV3-23",
                         y = "Count", orderBy = "P.value"))


Idents(BCR) <- "sample"

DoHeatmap(BCR, features = rownames(FindAllMarkers(BCR, min.pct = 0.25, logfc.threshold = 0.58, only.pos = T) %>%
  subset(p_val_adj < 0.05) %>%
  arrange(p_val_adj)), group.by = "sample")

#IGHV1-18: moderate124_Patient6 has 2 cells, severe122_Patient6 0
vg <- BCR[,BCR$v_gene == "IGHV1-18"]
vg$sample <-  droplevels(vg$sample)
vg <- ScaleData(vg)
mrk <- FindAllMarkers(vg, min.pct = 0.25, logfc.threshold = 0.58) %>%
  subset(p_val_adj < 0.05) %>%
  arrange(p_val_adj)
ggsave(filename = "graphs/DEGs-IGHV1-18.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(vg, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                       genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-IGHV1-18.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
           title = "IGHV1-18",
           y = "Count", orderBy = "P.value"))

#IGHV4-39: mild186_Patient3 has 1 cell
vg <- BCR[,BCR$v_gene == "IGHV4-39"]
vg$sample <-  droplevels(vg$sample)
vg <- ScaleData(vg)
mrk <- FindAllMarkers(vg, min.pct = 0.25, logfc.threshold = 0.58) %>%
  subset(p_val_adj < 0.05) %>%
  arrange(p_val_adj)
ggsave(filename = "graphs/DEGs-IGHV4-39.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(vg, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                       genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-IGHV4-39.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
           title = "IGHV4-39",
           y = "Count", orderBy = "P.value"))

#IGHV4-34: severe123_Patient5  has 2 cells
vg <- BCR[,BCR$v_gene == "IGHV4-34"]
vg$sample <-  droplevels(vg$sample)
vg <- ScaleData(vg)
mrk <- FindAllMarkers(vg, min.pct = 0.25, logfc.threshold = 0.58) %>%
  subset(p_val_adj < 0.05) %>%
  arrange(p_val_adj)
ggsave(filename = "graphs/DEGs-IGHV4-34.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(vg, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                       genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-IGHV4-34.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
           title = "IGHV4-34",
           y = "Count", orderBy = "P.value"))

#IGHV3-23: mild186_Patient3 has 1 cell
vg <- BCR[,BCR$v_gene == "IGHV3-23"]
vg$sample <-  droplevels(vg$sample)
vg <- ScaleData(vg)
mrk <- FindAllMarkers(vg, min.pct = 0.25, logfc.threshold = 0.58, only.pos = T) %>%
  subset(p_val_adj < 0.05) %>%
  arrange(p_val_adj)
ggsave(filename = "graphs/DEGs-IGHV3-23.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(vg, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021", #None significant
                       genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-IGHV3-23.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
           title = "IGHV3-23",
           y = "Count", orderBy = "P.value"))

#IGHV3-33: 
vg <- BCR[,BCR$v_gene == "IGHV3-33"]
vg$sample <-  droplevels(vg$sample)
vg <- ScaleData(vg)
mrk <- FindAllMarkers(vg, min.pct = 0.25, logfc.threshold = 0.58, only.pos = T) %>%
  subset(p_val_adj < 0.05) %>%
  arrange(p_val_adj)
ggsave(filename = "graphs/DEGs-IGHV3-33.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(vg, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                       genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-IGHV3-33.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
           title = "IGHV3-33",
           y = "Count", orderBy = "P.value"))

#IGHV3-30: moderate272_Patient1 has 2 cells, moderate124_Patient6 2
vg <- BCR[,BCR$v_gene == "IGHV3-30"]
vg$sample <-  droplevels(vg$sample)
vg <- ScaleData(vg)
mrk <- FindAllMarkers(vg, min.pct = 0.25, logfc.threshold = 0.58, only.pos = T) %>%
  subset(p_val_adj < 0.05) %>%
  arrange(p_val_adj)
ggsave(filename = "graphs/DEGs-IGHV3-30.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(vg, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                       genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-IGHV3-30.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
           title = "IGHV3-33",
           y = "Count", orderBy = "P.value"))

#IGHV3-48: moderate138_Patient5 has 2 cells
vg <- BCR[,BCR$v_gene == "IGHV3-48"]
vg$sample <-  droplevels(vg$sample)
vg <- ScaleData(vg)
mrk <- FindAllMarkers(vg, min.pct = 0.25, logfc.threshold = 0.58, only.pos = T) %>%
  subset(p_val_adj < 0.05) %>%
  arrange(p_val_adj)
ggsave(filename = "graphs/DEGs-IGHV3-48.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(vg, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021", #None statistically significant
                       genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-IGHV3-48.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
           title = "IGHV3-48",
           y = "Count", orderBy = "P.value"))

#IGHV4-59: 
vg <- BCR[,BCR$v_gene == "IGHV4-59"]
vg$sample <-  droplevels(vg$sample)
vg <- ScaleData(vg)
mrk <- FindAllMarkers(vg, min.pct = 0.25, logfc.threshold = 0.58) %>%
  subset(p_val_adj < 0.05) %>%
  arrange(p_val_adj)
ggsave(filename = "graphs/DEGs-IGHV4-59.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(vg, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                    genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-IGHV4-59.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
                         title = "IGHV4-59",
                         y = "Count", orderBy = "P.value"))

#IGHV1-69D: moderate303_Patient2 has 2 cells
vg <- BCR[,BCR$v_gene == "IGHV1-69D"]
vg$sample <-  droplevels(vg$sample)
vg <- ScaleData(vg)
mrk <- FindAllMarkers(vg, min.pct = 0.25, logfc.threshold = 0.58) %>%
  subset(p_val_adj < 0.05) %>%
  arrange(p_val_adj)
ggsave(filename = "graphs/DEGs-IGHV1-69D.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(vg, features = rownames(mrk), group.by = "sample"))
enriched <- enrichr(databases = "GO_Biological_Process_2021",
                    genes = rownames(mrk))
ggsave(filename = "graphs/enrichment-IGHV1-69D.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched[[1]], numChar = 70, 
                         title = "IGHV1-69D",
                         y = "Count", orderBy = "P.value"))


#Patient1 - IGHV3-33:
BCR$compare <- 0
BCR$compare[BCR$v_gene == "IGHV3-33" & BCR$patient == "Patient 1"] <- 1
Idents(BCR) <- "compare"
degs333 <- FindMarkers(BCR, logfc.threshold = 0.58, ident.1 = 1, min.pct = 0.25) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj) 
write.table(x = degs118, file = "DEGs-pt1-IGHV1-18.tsv", sep = "\t", col.names = NA)
ggsave(filename = "graphs/DEGs-pt1-IGHV3-33.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(BCR, group.by = "sample", size = 6, 
                        features = rownames(degs333)))
enriched333 <- enrichr(databases = "GO_Biological_Process_2021",
                    genes = rownames(degs333))
ggsave(filename = "graphs/enrichment-pt1-IGHV3-33.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched333[[1]], showTerms = 9, numChar = 70, 
                         title = "Patient 1 - IGHV3-33",
                         y = "Count", orderBy = "P.value"))

#Patient1 - IGHV1-18:
BCR$compare <- 0
BCR$compare[BCR$v_gene == "IGHV1-18" & BCR$patient == "Patient 1"] <- 1
Idents(BCR) <- "compare"
degs118 <- FindMarkers(BCR, logfc.threshold = 0.58, ident.1 = 1, min.pct = 0.25) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj) 
write.table(x = degs118, file = "DEGs-pt1-IGHV1-18.tsv", sep = "\t", col.names = NA)
ggsave(filename = "graphs/DEGs-pt1-IGHV1-18.jpeg", dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(BCR, group.by = "sample", size = 6, 
                        features = rownames(degs118)))
enriched118 <- enrichr(databases = "GO_Biological_Process_2021",
                    genes = rownames(degs118))
ggsave(filename = "graphs/enrichment-pt1-IGHV1-18.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched118[[1]], showTerms = 23, numChar = 70, 
                         title = "Patient 1 - IGHV1-18",
                         y = "Count", orderBy = "P.value"))


#Patient2 - IGHV4-34
BCR$compare <- 0
BCR$compare[BCR$v_gene == "IGHV4-34" & BCR$patient == "Patient 2"] <- 1
Idents(BCR) <- "compare"
degs434 <- FindMarkers(BCR, logfc.threshold = 0.58, ident.1 = 1, min.pct = 0.25) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)
degs434 <- degs434[!grepl("IGH", rownames(degs434)) & !grepl("IGL", rownames(degs434)) & !grepl("IGK", rownames(degs434)),]
write.table(x = degs434, file = "DEGs-pt2-IGHV4-34.tsv", sep = "\t", col.names = NA)
ggsave(filename = "graphs/DEGs-pt2-IGHV4-34.jpeg",  dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(BCR, group.by = "sample", size = 6, 
                        features = rownames(degs434[1:34,])))
enriched434 <- enrichr(databases = "GO_Biological_Process_2021",
                       genes = rownames(degs434))
ggsave(filename = "graphs/enrichment-pt2-IGHV4-34.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched434[[1]], showTerms = 30, numChar = 70, 
                         title = "Patient 2 - IGHV4-34",
                         y = "Count", orderBy = "P.value"))

#Patient2 - IGHV4-59
BCR$compare <- 0
BCR$compare[BCR$v_gene == "IGHV4-59" & BCR$patient == "Patient 2"] <- 1
Idents(BCR) <- "compare"
degs459 <- FindMarkers(BCR, logfc.threshold = 0.58, ident.1 = 1, min.pct = 0.25) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)
degs459 <- degs459[!grepl("IGH", rownames(degs459)) & !grepl("IGL", rownames(degs459)) & !grepl("IGK", rownames(degs459)),]
write.table(x = degs459, file = "DEGs-pt2-IGHV4-59.tsv", sep = "\t", col.names = NA)
ggsave(filename = "graphs/DEGs-pt2-IGHV4-59.jpeg",  dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(BCR, group.by = "sample", size = 6, 
                        features = rownames(degs459[1:34,])))
enriched459 <- enrichr(databases = "GO_Biological_Process_2021",
                       genes = rownames(degs459))
ggsave(filename = "graphs/enrichment-pt2-IGHV4-59.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched459[[1]], showTerms = 30, numChar = 70, 
                         title = "Patient 2 - IGHV4-59",
                         y = "Count", orderBy = "P.value"))

#Patient3 - IGHV4-39
BCR$compare <- 0
BCR$compare[BCR$v_gene == "IGHV4-39" & BCR$patient == "Patient 3"] <- 1
Idents(BCR) <- "compare"
degs439 <- FindMarkers(BCR, logfc.threshold = 0.58, ident.1 = 1, min.pct = 0.25) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)
degs439 <- degs439[!grepl("IGH", rownames(degs439)) & !grepl("IGL", rownames(degs439)) & !grepl("IGK", rownames(degs439)),]
write.table(x = degs439, file = "DEGs-pt3-IGHV4-39.tsv", sep = "\t", col.names = NA)
ggsave(filename = "graphs/DEGs-pt3-IGHV4-39.jpeg",  dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(BCR, group.by = "sample", size = 6, 
                        features = rownames(degs439[1:34,])))
enriched439 <- enrichr(databases = "GO_Biological_Process_2021",
                       genes = rownames(degs439))
ggsave(filename = "graphs/enrichment-pt3-IGHV4-39.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched439[[1]], showTerms = 30, numChar = 70, 
                         title = "Patient 3 - IGHV4-39",
                         y = "Count", orderBy = "P.value"))

#Patient3 - IGHV3-48
BCR$compare <- 0
BCR$compare[BCR$v_gene == "IGHV3-48" & BCR$patient == "Patient 3"] <- 1
Idents(BCR) <- "compare"
degs348 <- FindMarkers(BCR, logfc.threshold = 0.58, ident.1 = 1, min.pct = 0.25) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)
degs348 <- degs348[!grepl("IGH", rownames(degs348)) & !grepl("IGL", rownames(degs348)) & !grepl("IGK", rownames(degs348)),]
write.table(x = degs348, file = "DEGs-pt3-IGHV3-48.tsv", sep = "\t", col.names = NA)
ggsave(filename = "graphs/DEGs-pt3-IGHV3-48.jpeg",  dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(BCR, group.by = "sample", size = 6, 
                        features = rownames(degs348[1:34,])))
enriched348 <- enrichr(databases = "GO_Biological_Process_2021",
                       genes = rownames(degs348))
ggsave(filename = "graphs/enrichment-pt3-IGHV3-48.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched434[[1]], showTerms = 30, numChar = 70, 
                         title = "Patient 3 - IGHV3-48",
                         y = "Count", orderBy = "P.value"))

#Patient4 - IGHV3-23
BCR$compare <- 0
BCR$compare[BCR$v_gene == "IGHV3-23" & BCR$patient == "Patient 4"] <- 1
Idents(BCR) <- "compare"
degs3234 <- FindMarkers(BCR, logfc.threshold = 0.58, ident.1 = 1, min.pct = 0.25) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)
degs3234 <- degs3234[!grepl("IGH", rownames(degs3234)) & !grepl("IGL", rownames(degs323)) & !grepl("IGK", rownames(degs3234)),]
write.table(x = degs323, file = "DEGs-pt4-IGHV3-23.tsv", sep = "\t", col.names = NA)
ggsave(filename = "graphs/DEGs-pt4-IGHV3-23.jpeg",  dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(BCR, group.by = "sample", size = 6, 
                        features = rownames(degs3234[1:34,])))
enriched3234 <- enrichr(databases = "GO_Biological_Process_2021",
                       genes = rownames(degs3234))
ggsave(filename = "graphs/enrichment-pt4-IGHV3-23.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched3234[[1]], showTerms = 30, numChar = 70, 
                         title = "Patient 4 - IGHV3-23",
                         y = "Count", orderBy = "P.value"))

#Patient5 - IGHV3-23
BCR$compare <- 0
BCR$compare[BCR$v_gene == "IGHV3-23" & BCR$patient == "Patient 5"] <- 1
Idents(BCR) <- "compare"
degs3235 <- FindMarkers(BCR, logfc.threshold = 0.58, ident.1 = 1, min.pct = 0.25) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)
degs3235 <- degs3235[!grepl("IGH", rownames(degs3235)) & !grepl("IGL", rownames(degs3235)) & !grepl("IGK", rownames(degs3235)),]
write.table(x = degs3235, file = "DEGs-pt5-IGHV3-23.tsv", sep = "\t", col.names = NA)
ggsave(filename = "graphs/DEGs-pt5-IGHV3-23.jpeg",  dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(BCR, group.by = "sample", size = 6, 
                        features = rownames(degs3235[1:34,])))
enriched3235 <- enrichr(databases = "GO_Biological_Process_2021",
                        genes = rownames(degs3235))
ggsave(filename = "graphs/enrichment-pt5-IGHV3-23.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched3235[[1]], showTerms = 30, numChar = 70, 
                         title = "Patient 5 - IGHV3-23",
                         y = "Count", orderBy = "P.value"))

#Patient6 - IGHV3-23
BCR$compare <- 0
BCR$compare[BCR$v_gene == "IGHV3-23" & BCR$patient == "Patient 6"] <- 1
Idents(BCR) <- "compare"
degs3236 <- FindMarkers(BCR, logfc.threshold = 0.58, ident.1 = 1, min.pct = 0.25) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)
degs3236 <- degs3236[!grepl("IGH", rownames(degs3236)) & !grepl("IGL", rownames(degs3236)) & !grepl("IGK", rownames(degs3236)),]
write.table(x = degs3236, file = "DEGs-pt6-IGHV3-23.tsv", sep = "\t", col.names = NA)
ggsave(filename = "graphs/DEGs-pt6-IGHV3-23.jpeg",  dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(BCR, group.by = "sample", size = 6, 
                        features = rownames(degs3236[1:34,])))
enriched3236 <- enrichr(databases = "GO_Biological_Process_2021",
                        genes = rownames(degs3236))
ggsave(filename = "graphs/enrichment-pt6-IGHV3-23.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enriched3236[[1]], showTerms = 30, numChar = 70, 
                         title = "Patient 6 - IGHV3-23",
                         y = "Count", orderBy = "P.value"))


#Comparing Patient1 to the other two dead ones:
Idents(BCR) <- "patient"
degsdd <- FindMarkers(BCR[,BCR$outcome == "Deceased"], logfc.threshold = 0.58, min.pct = 0.25,
                      ident.1 = "Patient 1") %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)

DoHeatmap(BCR, features = rownames(degsdd[1:40,]), group.by = "sample")


Idents(BCR) <- "sample"
for (i in levels(factor(BCR$azimuthNames))) {
  for (j in levels(factor(BCR$sample))) {
    try({
      degs <- FindMarkers(BCR[,BCR$azimuthNames == i], logfc.threshold = 0.58, ident.1 = j, min.pct = 0.25) %>%
        subset(p_val_adj < 0.05) %>% arrange(p_val_adj) 
      write.table(x = degs, file = paste0("DEGs-", i, "-", j, ".tsv"), sep = "\t", col.names = NA)
    })
    
    
    
    
    
  }
  
}

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

# summarize the top abundant marker features for each group
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




#For IGHV1-18 in Patient 1:
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

# summarize the top abundant marker features for each group
top_markers(wlx.mrk.pt1)

write.table(top_markers(wlx.mrk.pt1), file = "top-wilcox-markers-v-1-18-pt1.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk.pt1 %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
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
            file = "ranked-genes-rrho-v-1-18-pt1.tsv", 
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

ggsave(filename = "graphs/wilcox-rrho-v-1-18-pt1.jpeg", dpi = "print",
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


################################V Gene Wilcoxon Analysis########################
#For IGHV1-18 in Patient 1:
pt1 <- BCR[,BCR$patient == "Patient 1" & BCR$v_gene == "IGHV1-18" & BCR$azimuthNames == "B intermediate"]
pt1$sample <- droplevels(pt1$sample)
DefaultAssay(pt1) <- "RNA"

# perform a fast Wilcoxon rank sum test with presto

wlx.mrk.pt1 <- wilcoxauc(pt1, 'sample', c("moderate272_Patient1", "critical293_Patient1"),
                         seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(padj)

head(wlx.mrk.pt1)

# we have all the genes for each cluster
dplyr::count(wlx.mrk.pt1, group)

# summarize the top abundant marker features for each group
top_markers(wlx.mrk.pt1)

write.table(top_markers(wlx.mrk.pt1), file = "top-wilcox-markers-b-int-v-1-18-pt1.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk.pt1 %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank
#head(n = 10)

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-b-int-v-1-18-pt1.tsv", 
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
            file = "ranked-genes-rrho-b-int-v-1-18-pt1.tsv", 
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

ggsave(filename = "graphs/wilcox-rrho-b-int-v-1-18-pt1.jpeg", dpi = "print",
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


ggsave(filename = "graphs/enrichment-b-int-v-1-18-pt1-wlx-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod1[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 1 moderate (IGHV1-18 - B intermediate)"))

ggsave(filename = "graphs/enrichment-b-int-v-1-18-pt1-wlx-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit1[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 1 critical (IGHV1-18 - B intermediate)"))



#For Patient 2:
pt2 <- BCR[,BCR$patient == "Patient 2" & BCR$v_gene == "IGHV4-34" & BCR$azimuthNames == "Plasmablast"]
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

write.table(top_markers(wlx.mrk.pt2), file = "top-wilcox-markers-plasm-v-4-34-pt2.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk.pt2 %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-plasm-v-4-34-pt2.tsv", 
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
            file = "ranked-genes-rrho-plasm-v-4-34-pt2.tsv", 
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

ggsave(filename = "graphs/wilcox-rrho-plasm-v-4-34-pt2.jpeg", dpi = "print",
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


ggsave(filename = "graphs/enrichment-plasm-v-4-34-pt2-wlx-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit2[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 2 critical (IGHV4-34 - Plasmablasts)"))

ggsave(filename = "graphs/enrichment-plasm-v-4-34-pt2-wlx-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod2[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 2 moderate (IGHV4-34 - Plasmablasts)"))


rm(list=setdiff(ls(), "BCR"))


#For Patient 3 IGHV4-39:
pt3 <- BCR[,BCR$patient == "Patient 3" & BCR$v_gene == "IGHV4-39" & BCR$azimuthNames == "Plasmablast"]
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

write.table(top_markers(wlx.mrk.pt3), file = "top-wilcox-markers-plasm-v-4-39-pt3.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk.pt3 %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-plasm-v-4-39-pt3.tsv", 
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
            file = "ranked-genes-rrho-plasm-v-4-39-pt3.tsv", 
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

ggsave(filename = "graphs/wilcox-rrho-plasm-v-4-39-pt3.jpeg", dpi = "print",
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


ggsave(filename = "graphs/enrichment-plasm-v-4-39-pt3-wlx-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit3[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 3 critical (IGHV4-39 - Plasmablasts)"))

ggsave(filename = "graphs/enrichment-plasm-v-4-39-pt3-wlx-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod3[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 3 moderate (IGHV4-39 - Plasmablasts)"))


rm(list=setdiff(ls(), "BCR"))

#For Patient 4 IGHV3-23:
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

write.table(top_markers(wlx.mrk.pt4), file = "top-wilcox-markers-v-3-23-pt4.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk.pt4 %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
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
            file = "ranked-genes-rrho-v-3-23-pt4.tsv", 
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

ggsave(filename = "graphs/wilcox-rrho-v-3-23-pt4.jpeg", dpi = "print",
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


ggsave(filename = "graphs/enrichment-v-3-23-pt4-wlx-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit4[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 4 critical (IGHV3-23)"))

ggsave(filename = "graphs/enrichment-v-3-23-pt4-wlx-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod4[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 4 moderate (IGHV3-23)"))


rm(list=setdiff(ls(), "BCR"))


#For Patient 4 IGHV3-23 in B intermediate cells:
pt4 <- BCR[,BCR$patient == "Patient 4" & BCR$v_gene == "IGHV3-23" & BCR$azimuthNames == "B intermediate"]
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

write.table(top_markers(wlx.mrk.pt4), file = "top-wilcox-markers-b-int-v-3-23-pt4.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk.pt4 %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-b-int-v-3-23-pt4.tsv", 
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
            file = "ranked-genes-rrho-v-b-int-3-23-pt4.tsv", 
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

pt4 <- ScaleData(object = pt4, 
                 features = as.character(unique(top$feature), assay = "RNA"))

ggsave(filename = "graphs/wilcox-rrho-b-int-v-3-23-pt4.jpeg", dpi = "print",
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


ggsave(filename = "graphs/enrichment-v-3-23-pt4-wlx-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit4[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 4 critical (IGHV3-23)"))

ggsave(filename = "graphs/enrichment-v-3-23-pt4-wlx-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod4[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 4 moderate (IGHV3-23)"))


rm(list=setdiff(ls(), "BCR"))


#For Patient 5 and IGHV3-23:
pt5 <- BCR[,BCR$patient == "Patient 5" & BCR$v_gene == "IGHV3-23"]
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

write.table(top_markers(wlx.mrk.pt5), file = "top-wilcox-markers-v-3-23-pt5.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk.pt5 %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-v-3-23-pt5.tsv", 
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
            file = "ranked-genes-rrho-v-3-23-pt5.tsv", 
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

ggsave(filename = "graphs/wilcox-rrho-v-3-23-pt5.jpeg", dpi = "print",
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


ggsave(filename = "graphs/enrichment-v-3-23-pt5-wlx-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit5[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 5 critical (IGHV3-23)"))

ggsave(filename = "graphs/enrichment-v-3-23-pt5-wlx-sev.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpsev5[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 5 severe (IGHV3-23)"))

ggsave(filename = "graphs/enrichment-v-3-23-pt5-wlx-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod5[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 5 moderate (IGHV3-23)"))



#For Patient 5 and IGHV3-23 plasmablasts:
pt5 <- BCR[,BCR$patient == "Patient 5" & BCR$v_gene == "IGHV3-23" & BCR$azimuthNames == "Plasmablast"]
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

write.table(top_markers(wlx.mrk.pt5), file = "top-wilcox-markers-plasm-v-3-23-pt5.tsv", sep = "\t", col.names = NA)

# rank genes using rrho algorithm
ranked.genes <- wlx.mrk.pt5 %>%
  mutate(rank = -log10(padj) * sign(logFC)) %>%   # rank genes by strength of significance, keeping the direction of the fold change
  arrange(-rank)                                  # sort genes by decreasing rank

# save found markers as exportable table
write.table(ranked.genes,
            file = "ranked-genes-plasm-v-3-23-pt5.tsv", 
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
            file = "ranked-genes-rrho-plasm-v-3-23-pt5.tsv", 
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

ggsave(filename = "graphs/wilcox-rrho-plasm-v-3-23-pt5.jpeg", dpi = "print",
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


ggsave(filename = "graphs/enrichment-v-3-23-pt5-wlx-crit.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpcrit5[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 5 critical (IGHV3-23)"))

ggsave(filename = "graphs/enrichment-v-3-23-pt5-wlx-sev.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpsev5[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 5 severe (IGHV3-23)"))

ggsave(filename = "graphs/enrichment-v-3-23-pt5-wlx-mod.jpeg", dpi = "print",
       height = 10, width = 10, units = "in",
       plot = plotEnrich(enrich.topUpmod5[[1]], 
                         numChar = 80, y = "Count", orderBy = "Adjusted.P.value",
                         xlab = NULL,
                         ylab = NULL,
                         title = "Upregulated in Patient 5 moderate (IGHV3-23)"))


rm(list=setdiff(ls(), "BCR"))
gc()

#For Patient 6 IGHV3-23:
pt6 <- BCR[,BCR$patient == "Patient 6" & BCR$v_gene == "IGHV3-23"]
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


