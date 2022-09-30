#Laoding:
library(presto)
library(ComplexHeatmap)
library(tidyverse)
BCR <- readRDS("05-BCR-combined.rds")

###For all:

for (i in levels(factor(BCR$patient[BCR$severity != "healthy"]))) {
  obj <- BCR[,BCR$patient == i & BCR$severity != "severe"]
  DefaultAssay(obj) <- "RNA"
  obj$severity[obj$severity == "mild"] <- "moderate"
  obj$sample <- droplevels(obj$sample)
  
  diff.genes <- wilcoxauc(obj, group_by = "severity", seurat_assay = "RNA", assay = "data") %>%
    subset(padj < 0.05) %>% arrange(desc(abs(logFC)))
  
  MNP.genes <- diff.genes %>%
    mutate(rank = -log10(pval) * sign(logFC)) %>%
    arrange(-rank)
  topUpmod <- MNP.genes %>% 
    dplyr::filter(group == c('moderate')) %>%
    filter(rank > 0) %>% 
    subset(padj < 0.05) %>% 
    top_n(50, wt=-padj)
  topUpcrit <- MNP.genes %>% 
    dplyr::filter(group == c('critical')) %>%
    filter(rank > 0) %>% 
    subset(padj < 0.05) %>%
    top_n(50, wt=-padj)
  topPathways <- bind_rows(topUpmod, topUpcrit) %>%
    filter(abs(avgExpr) > 0.5)
  #obj <- obj[ISGs[ISGs %in% topPathways$feature],]
  obj <- NormalizeData(obj)
  obj <- ScaleData(obj, features = topPathways$feature)
  expr <- obj@assays$RNA@scale.data
  #rownames(expr) <- gsub(x = rownames(expr), pattern = "HLA-", replacement = "HLA.")
  #diff.genes$feature <- gsub(x = diff.genes$feature, pattern = "HLA-", replacement = "HLA.")
  #rownames(expr) <- gsub(x = rownames(expr), pattern = "HLA-", replacement = "HLA.")
  #diff.genes$feature <- gsub(x = diff.genes$feature, pattern = "HLA-", replacement = "HLA.")
  meta <- obj@meta.data %>% group_by(patient, severity)
  meta <- data.frame(t(expr), meta)
  meta <- meta[,c(colnames(meta)[colnames(meta) %in% topPathways$feature], "patient", "severity")]

  expr.cols <- colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow"))
  col.anno.m <- HeatmapAnnotation(Severity = meta$severity[meta$severity == "moderate"], show_annotation_name = F, 
                                col = list(Severity = c("moderate" = "blue")))
  
  col.anno.c <- HeatmapAnnotation(Severity = meta$severity[meta$severity == "critical"], show_annotation_name = F, 
                                col = list(Severity = c("critical" = "red")))
  
  
  column_labels <- levels(meta$severity) 
  
  hm.m <- Heatmap(t(meta[meta$severity == "moderate",colnames(meta) %in% diff.genes$feature]), 
                    name = "Expression", row_names_gp = gpar(fontsize = 7),
                    column_labels = column_labels, 
                    col = expr.cols, show_column_names = F, top_annotation = col.anno.m,
                    show_row_dend = F, show_column_dend = F, cluster_columns = T) 
  
  
  hm.c <- Heatmap(t(meta[meta$severity == "critical",colnames(meta) %in% diff.genes$feature]), 
                    name = "Expression", row_names_gp = gpar(fontsize = 7),
                    column_labels = column_labels, 
                    col = expr.cols, show_column_names = F, top_annotation = col.anno.c,
                    show_row_dend = F, show_column_dend = F, cluster_columns = T) 
  
  
  
 
  
  
  tiff(paste0("./graphs/hm-critical-moderate-", i,".tiff"),
       res = 300, width = 5, height = 10, units = "in")
  if(levels(factor(meta$patient)) == "Patient 5" | levels(factor(meta$patient)) == "Patient 6") {
    draw(hm.c + hm.m)
  } else {
           draw(hm.m + hm.c)      
         }
     
  
  dev.off()
  
  
}
