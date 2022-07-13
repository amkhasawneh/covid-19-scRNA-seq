################################Loading#########################################
library(escape) #devtools::install_github("ncborcherding/escape@dev")
library(dittoSeq) #BiocManager::install("dittoSeq")
library(SingleCellExperiment)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(enrichR)

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
    group_by(v_gene, j_gene, sample) %>% na.omit() %>% dplyr::count() %>% arrange(desc(n)) %>%
    as.data.frame() -> matrix
  hiclo <- rbind(hiclo, matrix[1,])
  
}



degs[!grepl("IGH", rownames(degs)) & !grepl("IGL", rownames(degs)) & !grepl("IGK", rownames(degs)),]

BCR@meta.data %>% group_by(v_gene, j_gene) %>% count() %>% arrange(desc(n))




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
#IGHV1-18:
vg <- BCR[,BCR$v_gene == "IGHV1-18"]
vg$sample <-  droplevels(vg$sample)
vg <- ScaleData(vg)
mrk <- FindAllMarkers(vg, min.pct = 0.25, logfc.threshold = 0.58) %>%
  subset(p_val_adj < 0.05) %>%
  arrange(p_val_adj)
DoHeatmap(vg, features = rownames(mrk), group.by = "sample")

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


