BCR <- readRDS("05-BCR-combined.rds")

#Importing interferon-stimulated gene set:
ISGs <- read.csv("manuscript/ISGs_GeneSets_SU.csv")
ISGs <- ISGs$ISGs_geneset_227

#Cleaning gene set:
ISGs <- ISGs[!grepl(ISGs, pattern = "LOC") & !grepl(ISGs, pattern = "XENTR") & !grepl(ISGs, pattern = "GSON")]
ISGs <- ISGs[ISGs != ""]

for (i in levels(factor(BCR$patient[BCR$severity != "healthy"]))) {
 
   obj <- BCR[,BCR$patient == i & BCR$severity != "severe"]
   obj$severity[obj$severity == "mild"] <- "moderate"
   obj$sample <- droplevels(obj$sample)
   expr.cols <- colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow"))
   col.anno <- HeatmapAnnotation(Severity = meta$severity, show_annotation_name = F,   
                                 col = list(Severity = c("critical" = "red", "moderate" = "blue")))
   column_labels <- obj$severity
   for (j in levels(factor(obj$azimuthNames))) {
     diff.genes <- wilcoxauc(obj,  groupby = "severity", seurat_assay = "RNA", assay = "data") %>%
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
   topPathways <- bind_rows(topUpmod, topUpcrit)
   obj <- obj[ISGs[ISGs %in% topPathways$feature],]
   obj <- NormalizeData(obj)
   obj <- ScaleData(obj)
   expr <- obj@assays$RNA@scale.data
   meta <- obj@metadata %>% group_by(sample, severity) %>% count()
   meta <- data.frame(t(expr), meta)
   meta <- meta[,c(colnames(meta)[colnames(meta) %in% rownames(obj@assays$RNA@scale.data)], "sample", "severity")]
   hm <- Heatmap(t(meta[,colnames(meta) %in% diff.genes$feature]), 
                   column_split = meta$severity, name = "Expression", 
                   column_labels = column_labels, row_title = j, column_title = " ",  
                   col = expr.cols, top_annotation = col.anno, row_names_gp = gpar(fontsize = 7), 
                   show_row_dend = F, show_column_dend = F, cluster_columns = F)
   }
   
   hm.list <- hm.list %v% hm
   tiff(paste0("./graphs/hm-critical-moderate-b-cells-isg-in-deg-", i, ".tiff"),
        res = 300, width = 4, height = 8, units = "in")
   draw(hm.list)
   dev.off()
 
}
