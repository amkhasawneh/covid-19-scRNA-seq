################################Loading#########################################
library(Seurat)
library(presto)
library(msigdbr)
library(fgsea)
library(tidyverse)


BCR <- readRDS("05-BCR-combined.rds")


#Choosing the human Hallmark annotated gene set from the MsigDB:
m_df <- msigdbr(species = "Homo sapiens", category = "H")


#Preparing a list of the gene sets, for fgsea:
fgsea_sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

################################Gene Info - Wilcoxon############################



#Printing out the most abundant V-J combinations for each sample:
hiclo <- NULL
for (i in levels(factor(BCR@meta.data$patient[BCR$severity != "healthy"]))) {
  BCR@meta.data[BCR$patient == i & !is.na(BCR$v_gene),] %>%
    group_by(v_gene, patient)  %>% dplyr::count() %>% na.omit() %>%
    arrange(desc(n)) %>% as.data.frame() -> matrix
  hiclo <- rbind(hiclo, matrix[1,])
  
}


for (i in levels(factor(BCR$patient[BCR$severity != "healthy"]) %>% droplevels())) {
  #Subsetting the data, based on patient and top V gene:
  obj <- BCR[,BCR$patient == i & BCR$v_gene == hiclo$v_gene[hiclo$patient == i]]
  obj$sample <- droplevels(obj$sample)
  DefaultAssay(obj) <- "RNA"
  
  if (levels(factor(obj$patient)) == "Patient 5" | levels(factor(obj$patient)) == "Patient 6" ) {
    #Running the "quick" Wilcoxon test to extract DEGs:
    diff.genes <- wilcoxauc(obj, "severity", c("critical", "severe"),
                            seurat_assay='RNA', assay = "data") %>%
      subset(padj < 0.05) %>% arrange(padj)
    
    #Separating moderate genes from critical genes:
    MNP.genes <- diff.genes %>%
      mutate(rank = -log10(pval) * sign(logFC)) %>%
      arrange(-rank)  
    
    #Selecting only the feature and rank columns of DGE data for fgsea run:
    MNP.genes_c <- MNP.genes %>%
      dplyr::filter(group == levels(as.factor(diff.genes$group))[1]) %>%
      arrange(-rank) %>% 
      dplyr::select(feature, rank)
    
    MNP.genes_s <- MNP.genes %>%
      dplyr::filter(group == levels(as.factor(diff.genes$group))[2]) %>%
      arrange(-rank) %>% 
      dplyr::select(feature, rank)
    
    
    #Converting matrices to tables:
    ranks_c <- deframe(MNP.genes_c)
    ranks_s <- deframe(MNP.genes_s)
    
    #Running fgsea based on gene set ranking (stats = ranks):
    fgseaRes_c <- fgseaMultilevel(fgsea_sets, stats = ranks_c, nPermSimple = 1000, 
                                  maxSize = 200) %>% arrange(padj)
    fgseaRes_s <- fgseaMultilevel(fgsea_sets, stats = ranks_s, nPermSimple = 1000, 
                                  maxSize = 200) %>% arrange(padj)
    
    
    
    #Tidying the data:
    fgseaResTidy_c <- fgseaRes_c %>%
      as_tibble() %>%
      arrange(desc(NES))
    
    fgseaResTidy_s <- fgseaRes_s %>%
      as_tibble() %>%
      arrange(desc(NES))
    
    
    # choose top upregulated genes
    topUpcrit <- MNP.genes %>% 
      dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
      filter(rank > 0) %>% 
      top_n(50, wt=-padj)
    
    topDowncrit <- MNP.genes %>% 
      dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
      filter(rank < 0) %>% 
      top_n(50, wt=-padj)
    
    topUpsev <- MNP.genes %>% 
      dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
      filter(rank > 0) %>% 
      top_n(60, wt=-padj)
    
    topDownsev <- MNP.genes %>%
      dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
      filter(rank < 0) %>%
      top_n(50, wt=-padj)
    
    
    
    
    
    topPathways <- bind_rows(topUpcrit, topDowncrit, topUpsev, topDownsev)
    
    top <- topPathways %>% group_by(group) 
    
    # create a scale.data slot for the selected genes in subset data
    obj <- obj[,obj$severity == "critical" | obj$severity == "severe"]
    obj$sample <- droplevels(obj$sample)
    
    alldata <- ScaleData(object = obj, 
                         features = as.character(unique(top$feature), assay = "RNA"))
    
    DefaultAssay(alldata) <- "RNA"
    
    ggsave(filename = paste0("graphs/wilcox-rrho-hm-", i, "-", levels(as.factor(MNP.genes$group))[1], "-",
                             levels(as.factor(MNP.genes$group))[2], ".tiff"), 
           dpi = 300, width = 15, height = 10, units = "in",
           plot = DoHeatmap(alldata,
                            features = unique(top$feature), group.by = "sample",
                            size = 4, raster=TRUE, label=FALSE) + 
             labs(title = paste0(i, " ", levels(as.factor(MNP.genes$group))[1], " vs ",
                                 levels(as.factor(MNP.genes$group))[2])) +
             theme(axis.text.y = element_text(size = 9),
                   plot.title = element_text(hjust = 0.5, size = 16))) 
    
    #Plotting pathways:
    fgseaResTidy_c$adjPvalue <- ifelse(fgseaResTidy_c$padj <= 0.05, "significant", "non-significant")
    fgseaResTidy_c$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_c$pathway)
    cols <- c("significant" = "red") #"non-significant" = "grey", 
    plot_c <- ggplot(fgseaResTidy_c, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
      geom_col() +
      scale_fill_manual(values = cols) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title = paste0(i, " ", levels(as.factor(diff.genes$group))[1])) + 
      theme(axis.text.x = element_text(vjust = 0.5, size = 10),
            axis.text.y = element_text(size = 10, vjust = 0.5),
            plot.title = element_text(hjust = 0.5, size = 14),
            axis.title = element_text(size = 10),
            legend.position='right', legend.key.size = unit(0.25, 'cm'),
            legend.key.height = unit(0.25, 'cm'),
            legend.key.width = unit(0.25, 'cm'), 
            legend.title = element_text(size=10),
            legend.text = element_text(size=10)) 
    ggsave(filename = paste0("graphs/gsea-",   i, "-", levels(as.factor(diff.genes$group))[1], "-vs-", levels(as.factor(diff.genes$group))[2],  ".tiff"),
           width = 10, height = 10, dpi = 300, units = "in", plot = plot_c)
  
    fgseaResTidy_s$adjPvalue <- ifelse(fgseaResTidy_s$padj <= 0.05, "significant", "non-significant")
    fgseaResTidy_s$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_s$pathway)
    cols <- c("significant" = "red") #"non-significant" = "grey", 
    plot_s <- ggplot(fgseaResTidy_s, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
      geom_col() +
      scale_fill_manual(values = cols) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title = paste0(i, " ", levels(as.factor(diff.genes$group))[2])) + 
      theme(axis.text.x = element_text(vjust = 0.5, size = 10),
            axis.text.y = element_text(size = 10, vjust = 0.5),
            plot.title = element_text(hjust = 0.5, size = 14),
            axis.title = element_text(size = 10),
            legend.position='none') 
    ggsave(filename = paste0("graphs/gsea-",  i, "-", levels(as.factor(diff.genes$group))[2], "-vs-", levels(as.factor(diff.genes$group))[1],  ".tiff"),
           width = 10, height = 10, dpi = 300, units = "in", plot = plot_s)
    
    obj <- BCR[,BCR$patient == i & BCR$v_gene == hiclo$v_gene[hiclo$patient == i]]
    obj$sample <- droplevels(obj$sample)
    DefaultAssay(obj) <- "RNA"
    
    #Running the "quick" Wilcoxon test to extract DEGs:
    diff.genes <- wilcoxauc(obj, "severity", c("critical", "moderate"),
                            seurat_assay='RNA', assay = "data") %>%
      subset(padj < 0.05) %>% arrange(padj)
    
    #Separating moderate genes from critical genes:
    MNP.genes <- diff.genes %>%
      mutate(rank = -log10(pval) * sign(logFC)) %>%
      arrange(-rank)  
    
    #Selecting only the feature and rank columns of DGE data for fgsea run:
    MNP.genes_c <- MNP.genes %>%
      dplyr::filter(group == levels(as.factor(diff.genes$group))[1]) %>%
      arrange(-rank) %>% 
      dplyr::select(feature, rank)
    
    MNP.genes_m <- MNP.genes %>%
      dplyr::filter(group == levels(as.factor(diff.genes$group))[2]) %>%
      arrange(-rank) %>% 
      dplyr::select(feature, rank)
    
    #Converting matrices to tables:
    ranks_c <- deframe(MNP.genes_c)
    ranks_m <- deframe(MNP.genes_m)
        
    #Running fgsea based on gene set ranking (stats = ranks):
    fgseaRes_c <- fgseaMultilevel(fgsea_sets, stats = ranks_c, nPermSimple = 1000, 
                                  maxSize = 200) %>% arrange(padj)
    fgseaRes_m <- fgseaMultilevel(fgsea_sets, stats = ranks_m, nPermSimple = 1000, 
                                  maxSize = 200) %>% arrange(padj)
    
    
    
    #Tidying the data:
    fgseaResTidy_c <- fgseaRes_c %>%
      as_tibble() %>%
      arrange(desc(NES))
    
    fgseaResTidy_m <- fgseaRes_m %>%
      as_tibble() %>%
      arrange(desc(NES))
    
    
    # choose top upregulated genes
     topUpcrit <- MNP.genes %>% 
      dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
      filter(rank > 0) %>% 
      top_n(50, wt=-padj)
    
    topDowncrit <- MNP.genes %>% 
      dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
      filter(rank < 0) %>% 
      top_n(50, wt=-padj)
    
    topUpmod <- MNP.genes %>% 
      dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
      filter(rank > 0) %>% 
      top_n(60, wt=-padj)
    
    topDownmod <- MNP.genes %>%
      dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
      filter(rank < 0) %>%
      top_n(50, wt=-padj)
    
   
    
    
    
    topPathways <- bind_rows(topUpcrit, topDowncrit, topUpmod, topDownmod)
    
    top <- topPathways %>% group_by(group) 
    
    # create a scale.data slot for the selected genes in subset data
    obj <- obj[,obj$severity == "critical" | obj$severity == "moderate"]
    obj$sample <- droplevels(obj$sample)
    
    alldata <- ScaleData(object = obj, 
                         features = as.character(unique(top$feature), assay = "RNA"))
    
    DefaultAssay(alldata) <- "RNA"
    
    ggsave(filename = paste0("graphs/wilcox-rrho-hm-", i, "-", levels(as.factor(MNP.genes$group))[1], "-",
                             levels(as.factor(MNP.genes$group))[2], ".tiff"), 
           dpi = 300, width = 15, height = 10, units = "in",
           plot = DoHeatmap(alldata,
                            features = unique(top$feature), group.by = "sample",
                            size = 4, raster=TRUE, label=FALSE) + 
             labs(title = paste0(i, " ", levels(as.factor(MNP.genes$group))[1], " vs ",
                                 levels(as.factor(MNP.genes$group))[2])) +
             theme(axis.text.y = element_text(size = 9),
                   plot.title = element_text(hjust = 0.5, size = 16)))                  
    
    #Plotting pathways:
    fgseaResTidy_c$adjPvalue <- ifelse(fgseaResTidy_c$padj <= 0.05, "significant", "non-significant")
    fgseaResTidy_c$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_c$pathway)
    cols <- c("significant" = "red") #"non-significant" = "grey", 
    plot_c <- ggplot(fgseaResTidy_c, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
      geom_col() +
      scale_fill_manual(values = cols) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title = paste0(i, " ", levels(as.factor(diff.genes$group))[1])) + 
      theme(axis.text.x = element_text(vjust = 0.5, size = 10),
            axis.text.y = element_text(size = 10, vjust = 0.5),
            plot.title = element_text(hjust = 0.5, size = 14),
            axis.title = element_text(size = 10),
            legend.position='right', legend.key.size = unit(0.25, 'cm'),
            legend.key.height = unit(0.25, 'cm'),
            legend.key.width = unit(0.25, 'cm'), 
            legend.title = element_text(size=10),
            legend.text = element_text(size=10)) 
    ggsave(filename = paste0("graphs/gsea-",  i, "-", levels(as.factor(diff.genes$group))[1], "-vs-", levels(as.factor(diff.genes$group))[2],  ".tiff"),
           width = 10, height = 10, dpi = 300, units = "in", plot = plot_c)
    
    fgseaResTidy_m$adjPvalue <- ifelse(fgseaResTidy_m$padj <= 0.05, "significant", "non-significant")
    fgseaResTidy_m$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_m$pathway)
    cols <- c("significant" = "red") #"non-significant" = "grey", 
    plot_m <- ggplot(fgseaResTidy_m, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
      geom_col() +
      scale_fill_manual(values = cols) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title = paste0(i, " ", levels(as.factor(diff.genes$group))[2])) + 
      theme(axis.text.x = element_text(vjust = 0.5, size = 10),
            axis.text.y = element_text(size = 10, vjust = 0.5),
            plot.title = element_text(hjust = 0.5, size = 14),
            axis.title = element_text(size = 10),
            legend.position='none') 
    ggsave(filename = paste0("graphs/gsea-",  i, "-", levels(as.factor(diff.genes$group))[2], "-vs-", levels(as.factor(diff.genes$group))[1], ".tiff"),
           width = 10, height = 10, dpi = 300, units = "in", plot = plot_m)
    
    obj <- BCR[,BCR$patient == i & BCR$v_gene == hiclo$v_gene[hiclo$patient == i]]
    obj$sample <- droplevels(obj$sample)
    DefaultAssay(obj) <- "RNA"
   
    
    #Running the "quick" Wilcoxon test to extract DEGs:
    diff.genes <- wilcoxauc(obj, "severity", c("severe", "moderate"),
                            seurat_assay='RNA', assay = "data") %>%
      subset(padj < 0.05) %>% arrange(padj)
    
    #Separating moderate genes from critical genes:
    MNP.genes <- diff.genes %>%
      mutate(rank = -log10(pval) * sign(logFC)) %>%
      arrange(-rank)  
    
    #Selecting only the feature and rank columns of DGE data for fgsea run:
    MNP.genes_s <- MNP.genes %>%
      dplyr::filter(group == levels(as.factor(diff.genes$group))[2]) %>%
      arrange(-rank) %>% 
      dplyr::select(feature, rank)
    
    MNP.genes_m <- MNP.genes %>%
      dplyr::filter(group == levels(as.factor(diff.genes$group))[1]) %>%
      arrange(-rank) %>% 
      dplyr::select(feature, rank)
    
    #Converting matrices to tables:
    ranks_s <- deframe(MNP.genes_s)
    ranks_m <- deframe(MNP.genes_m)
    
    #Running fgsea based on gene set ranking (stats = ranks):
    fgseaRes_s <- fgseaMultilevel(fgsea_sets, stats = ranks_s, nPermSimple = 1000, 
                                  maxSize = 200) %>% arrange(padj)
    fgseaRes_m <- fgseaMultilevel(fgsea_sets, stats = ranks_m, nPermSimple = 1000, 
                                  maxSize = 200) %>% arrange(padj)
    
    #Tidying the data:
    fgseaResTidy_s <- fgseaRes_s %>%
      as_tibble() %>%
      arrange(desc(NES))
    
    fgseaResTidy_m <- fgseaRes_m %>%
      as_tibble() %>%
      arrange(desc(NES))

    
    # choose top upregulated genes
    topUpmod <- MNP.genes %>% 
      dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
      filter(rank > 0) %>% 
      top_n(60, wt=-padj)
    
    topDownmod <- MNP.genes %>%
      dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
      filter(rank < 0) %>%
      top_n(50, wt=-padj)
    
    topUpsev <- MNP.genes %>% 
      dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
      filter(rank > 0) %>% 
      top_n(50, wt=-padj)
    
    topDownsev <- MNP.genes %>% 
      dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
      filter(rank < 0) %>% 
      top_n(50, wt=-padj)
    
    
    
    topPathways <- bind_rows(topUpsev, topDownsev, topUpmod, topDownmod)
    
    top <- topPathways %>% group_by(group) 
    
    # create a scale.data slot for the selected genes in subset data
    obj <- obj[,obj$severity == "severe" | obj$severity == "moderate"]
    obj$sample <- droplevels(obj$sample)
    
    alldata <- ScaleData(object = obj, 
                         features = as.character(unique(top$feature), assay = "RNA"))
    
    DefaultAssay(alldata) <- "RNA"
    
    ggsave(filename = paste0("graphs/wilcox-rrho-hm-", i, "-", levels(as.factor(MNP.genes$group))[1], "-",
                             levels(as.factor(MNP.genes$group))[2], ".tiff"), 
           dpi = 300, width = 15, height = 10, units = "in",
           plot = DoHeatmap(alldata,
                            features = unique(top$feature), group.by = "sample",
                            size = 4, raster=TRUE, label=FALSE) + 
             labs(title = paste0(i, " ", levels(as.factor(MNP.genes$group))[1], " vs ",
                                levels(as.factor(MNP.genes$group))[2])) +
             theme(axis.text.y = element_text(size = 9),
                   plot.title = element_text(hjust = 0.5, size = 16)))
    
    #Plotting pathways:
    fgseaResTidy_s$adjPvalue <- ifelse(fgseaResTidy_s$padj <= 0.05, "significant", "non-significant")
    fgseaResTidy_s$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_s$pathway)
    cols <- c("significant" = "red") #"non-significant" = "grey", 
    plot_s <- ggplot(fgseaResTidy_s, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
      geom_col() +
      scale_fill_manual(values = cols) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title = paste0(i, " ", "severe")) + 
      theme(axis.text.x = element_text(vjust = 0.5, size = 10),
            axis.text.y = element_text(size = 10, vjust = 0.5),
            plot.title = element_text(hjust = 0.5, size = 14),
            axis.title = element_text(size = 10),
            legend.position='none') 
    ggsave(filename = paste0("graphs/gsea-",  i, "-", "severe-vs-moderate.tiff"),
           width = 10, height = 10, dpi = 300, units = "in", plot = plot_s)
    
    fgseaResTidy_m$adjPvalue <- ifelse(fgseaResTidy_m$padj <= 0.05, "significant", "non-significant")
    fgseaResTidy_m$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_m$pathway)
    cols <- c("significant" = "red") #"non-significant" = "grey", 
    plot_m <- ggplot(fgseaResTidy_m, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
      geom_col() +
      scale_fill_manual(values = cols) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title = paste0(i, " ", "moderate")) + 
      theme(axis.text.x = element_text(vjust = 0.5, size = 10),
            axis.text.y = element_text(size = 10, vjust = 0.5),
            plot.title = element_text(hjust = 0.5, size = 14),
            axis.title = element_text(size = 10),
            legend.position='none') 
    ggsave(filename = paste0("graphs/gsea-",  i, "-", "moderate-vs-severe.tiff"),
           width = 10, height = 10, dpi = 300, units = "in", plot = plot_m)
    
    
  } else {
    #Running the "quick" Wilcoxon test to extract DEGs:
    diff.genes <- wilcoxauc(obj, "severity", 
                            seurat_assay='RNA', assay = "data") %>%
      subset(padj < 0.05) %>% arrange(padj)
    
    #Separating moderate genes from critical genes:
    MNP.genes <- diff.genes %>%
      mutate(rank = -log10(pval) * sign(logFC)) %>%
      arrange(-rank)  
    
    #Selecting only the feature and rank columns of DGE data for fgsea run:
    MNP.genes_c <- MNP.genes %>%
      dplyr::filter(group == levels(as.factor(diff.genes$group))[1]) %>%
      arrange(-rank) %>% 
      dplyr::select(feature, rank)
    
    MNP.genes_m <- MNP.genes %>%
      dplyr::filter(group == levels(as.factor(diff.genes$group))[2]) %>%
      arrange(-rank) %>% 
      dplyr::select(feature, rank)
    
    #Converting matrices to tables:
    ranks_c <- deframe(MNP.genes_c)
    ranks_m <- deframe(MNP.genes_m)
    
    #Running fgsea based on gene set ranking (stats = ranks):
    fgseaRes_c <- fgseaMultilevel(fgsea_sets, stats = ranks_c, nPermSimple = 1000, 
                                  maxSize = 200) %>% arrange(padj)
    fgseaRes_m <- fgseaMultilevel(fgsea_sets, stats = ranks_m, nPermSimple = 1000, 
                                  maxSize = 200) %>% arrange(padj)
    
    #Tidying the data:
    fgseaResTidy_m <- fgseaRes_m %>%
      as_tibble() %>%
      arrange(desc(NES))
    
    fgseaResTidy_c <- fgseaRes_c %>%
      as_tibble() %>%
      arrange(desc(NES))
    
    
    
    # choose top upregulated genes
    topUpcrit <- MNP.genes %>% 
      dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
      filter(rank > 0) %>% 
      top_n(50, wt=-padj)
    
    topDowncrit <- MNP.genes %>% 
      dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
      filter(rank < 0) %>% 
      top_n(50, wt=-padj)
    
    topUpmod <- MNP.genes %>% 
      dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
      filter(rank > 0) %>% 
      top_n(60, wt=-padj)
    
    topDownmod <- MNP.genes %>%
      dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
      filter(rank < 0) %>%
      top_n(50, wt=-padj)
    
    topPathways <- bind_rows(topUpcrit, topDowncrit, topUpmod, topDownmod)
    
    top <- topPathways %>% group_by(group) 
    
    # create a scale.data slot for the selected genes in subset data
    alldata <- ScaleData(object = obj[,obj$severity == levels(as.factor(MNP.genes$group))[1] | obj$severity == levels(as.factor(MNP.genes$group))[2]], 
                         features = as.character(unique(top$feature), assay = "RNA"))
    
    DefaultAssay(alldata) <- "RNA"
    
    ggsave(filename = paste0("graphs/wilcox-rrho-hm-", i, ".tiff"), dpi = 300,
           width = 15, height = 10, units = "in",
           plot = DoHeatmap(alldata,
                            features = unique(top$feature), group.by = "sample",
                            size = 4, raster=TRUE, label=FALSE) + 
             labs(title = i) +
             theme(axis.text.y = element_text(size = 9),
                   plot.title = element_text(hjust = 0.5, size = 16)))
    
    #Plotting pathways:
    fgseaResTidy_c$adjPvalue <- ifelse(fgseaResTidy_c$padj <= 0.05, "significant", "non-significant")
    fgseaResTidy_c$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_c$pathway)
    cols <- c("significant" = "red") #"non-significant" = "grey", 
    plot_c <- ggplot(fgseaResTidy_c, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
      geom_col() +
      scale_fill_manual(values = cols) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title = paste0(i, " ", levels(as.factor(diff.genes$group))[1])) + 
      theme(axis.text.x = element_text(vjust = 0.5, size = 10),
            axis.text.y = element_text(size = 10, vjust = 0.5),
            plot.title = element_text(hjust = 0.5, size = 14),
            axis.title = element_text(size = 10),
            legend.position='right', legend.key.size = unit(0.25, 'cm'),
            legend.key.height = unit(0.25, 'cm'),
            legend.key.width = unit(0.25, 'cm'), 
            legend.title = element_text(size=10),
            legend.text = element_text(size=10)) 
    ggsave(filename = paste0("graphs/gsea-",  i, "-", levels(as.factor(diff.genes$group))[1], ".tiff"),
           width = 10, height = 10, dpi = 300, units = "in", plot = plot_c)
    
    fgseaResTidy_m$adjPvalue <- ifelse(fgseaResTidy_m$padj <= 0.05, "significant", "non-significant")
    fgseaResTidy_m$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_m$pathway)
    cols <- c("significant" = "red") #"non-significant" = "grey", 
    plot_m <- ggplot(fgseaResTidy_m, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
      geom_col() +
      scale_fill_manual(values = cols) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title = paste0(i, " ", levels(as.factor(diff.genes$group))[2])) + 
      theme(axis.text.x = element_text(vjust = 0.5, size = 10),
            axis.text.y = element_text(size = 10, vjust = 0.5),
            plot.title = element_text(hjust = 0.5, size = 14),
            axis.title = element_text(size = 10),
            legend.position='none') 
    ggsave(filename = paste0("graphs/gsea-",  i, "-", levels(as.factor(diff.genes$group))[2], ".tiff"),
           width = 10, height = 10, dpi = 300, units = "in", plot = plot_m)
    
  }
  
}


################################Using the most commong V genes to all###########

rm(list = setdiff(ls(), c("BCR", "m_df", "fgsea_sets")))
for (v in c("IGHV3-23", "IGHV4-34", "IGHV4-59")) {
  
for (i in levels(factor(BCR$patient[BCR$severity != "healthy"]) %>% droplevels())) {
  #Subsetting the data, based on patient and top V gene:
  obj <- BCR[,BCR$patient == i & BCR$v_gene == v]
  obj$sample <- droplevels(obj$sample)
  DefaultAssay(obj) <- "RNA"
  try(
    {
      if (levels(factor(obj$patient)) == "Patient 5" | levels(factor(obj$patient)) == "Patient 6" ) {
        #Running the "quick" Wilcoxon test to extract DEGs:
        diff.genes <- wilcoxauc(obj, "severity", c("critical", "severe"),
                                seurat_assay='RNA', assay = "data") %>%
          subset(padj < 0.05) %>% arrange(padj)
        
        #Separating moderate genes from critical genes:
        MNP.genes <- diff.genes %>%
          mutate(rank = -log10(pval) * sign(logFC)) %>%
          arrange(-rank)  
        
        #Selecting only the feature and rank columns of DGE data for fgsea run:
        MNP.genes_c <- MNP.genes %>%
          dplyr::filter(group == levels(as.factor(diff.genes$group))[1]) %>%
          arrange(-rank) %>% 
          dplyr::select(feature, rank)
        
        MNP.genes_s <- MNP.genes %>%
          dplyr::filter(group == levels(as.factor(diff.genes$group))[2]) %>%
          arrange(-rank) %>% 
          dplyr::select(feature, rank)
        
        
        #Converting matrices to tables:
        ranks_c <- deframe(MNP.genes_c)
        ranks_s <- deframe(MNP.genes_s)
        
        #Running fgsea based on gene set ranking (stats = ranks):
        fgseaRes_c <- fgseaMultilevel(fgsea_sets, stats = ranks_c, nPermSimple = 1000, 
                                      maxSize = 200) %>% arrange(padj)
        fgseaRes_s <- fgseaMultilevel(fgsea_sets, stats = ranks_s, nPermSimple = 1000, 
                                      maxSize = 200) %>% arrange(padj)
        
        #Enrichment plots:
        ggsave(filename = paste0("graphs/enrichment-plot-", i, "-", levels(as.factor(MNP.genes$group))[1], "-",
                                 levels(as.factor(MNP.genes$group))[2], "-", v, ".tiff"), 
               dpi = 300, width = 10, height = 5, units = "in",
               plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_c, desc(NES))[1]$pathway]],
                       ranks_c) + 
          labs(title = arrange(fgseaRes_c, desc(NES))[1]$pathway) +
          xlab("Rank") + ylab("Enrichment Score") +
          theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
          theme(plot.title = element_text(hjust = 0.5, size = 16)))
        
        ggsave(filename = paste0("graphs/enrichment-plot-", i, "-", levels(as.factor(MNP.genes$group))[2], "-",
                                 levels(as.factor(MNP.genes$group))[1], "-", v, ".tiff"), 
               dpi = 300, width = 10, height = 5, units = "in",
               plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_s, desc(NES))[1]$pathway]],
                       ranks_s) + 
          labs(title = arrange(fgseaRes_s, desc(NES))[1]$pathway) +
          xlab("Rank") + ylab("Enrichment Score") +
          theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
          theme(plot.title = element_text(hjust = 0.5, size = 16)))
        
        #Tidying the data:
        fgseaResTidy_c <- fgseaRes_c %>%
          as_tibble() %>%
          arrange(desc(NES) ) %>%
          mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ", ", character(1L)))
        
        write.table(fgseaResTidy_c, file = paste0(i, "-critical-severe-", v, ".tsv"), sep = "\t", col.names = NA)
        
        
        
        fgseaResTidy_s <- fgseaRes_s %>%
          as_tibble() %>%
          arrange(desc(NES)) %>%
          mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ", ", character(1L)))
        
        write.table(fgseaResTidy_s, file = paste0(i, "-severe-critical-", v, ".tsv"), sep = "\t", col.names = NA)
        
                
        
        # choose top upregulated genes
        topUpcrit <- MNP.genes %>% 
          dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
          filter(rank > 0) %>% 
          top_n(50, wt=-padj)
        
        topDowncrit <- MNP.genes %>% 
          dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
          filter(rank < 0) %>% 
          top_n(50, wt=-padj)
        
        topUpsev <- MNP.genes %>% 
          dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
          filter(rank > 0) %>% 
          top_n(60, wt=-padj)
        
        topDownsev <- MNP.genes %>%
          dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
          filter(rank < 0) %>%
          top_n(50, wt=-padj)
        
        
        
        
        
        topPathways <- bind_rows(topUpcrit, topDowncrit, topUpsev, topDownsev)
        
        top <- topPathways %>% group_by(group) 
        
        #Enrichment table graphs:
        tiff(filename = paste0("graphs/gsea-table-", i, "-", v, "-", levels(as.factor(MNP.genes$group))[2], "-", levels(as.factor(MNP.genes$group))[1], ".tiff"), res = 300,
             width = 10, height = 5, units = "in")
        plotGseaTable(fgsea_sets[arrange(fgseaRes_s, desc(NES))[1]$pathway], 
                      ranks_s, 
                      fgseaResTidy_s, 
                      gseaParam = 0.5)
        dev.off()
        
        tiff(filename = paste0("graphs/gsea-table-", i, "-", v, "-", levels(as.factor(MNP.genes$group))[1], "-", levels(as.factor(MNP.genes$group))[2], ".tiff"), res = 300,
             width = 10, height = 5, units = "in")
        plotGseaTable(fgsea_sets[arrange(fgseaRes_c, desc(NES))[1]$pathway], 
                      ranks_c, 
                      fgseaResTidy_c, 
                      gseaParam = 0.5) + 
          theme(text = element_text(size = 15))
        dev.off()
        
        # create a scale.data slot for the selected genes in subset data
        obj <- obj[,obj$severity == "critical" | obj$severity == "severe"]
        obj$sample <- droplevels(obj$sample)
        
        alldata <- ScaleData(object = obj, 
                             features = as.character(unique(top$feature), assay = "RNA"))
        
        DefaultAssay(alldata) <- "RNA"
        
        ggsave(filename = paste0("graphs/wilcox-rrho-hm-", i, "-", levels(as.factor(MNP.genes$group))[1], "-",
                                 levels(as.factor(MNP.genes$group))[2], "-", v, ".tiff"), 
               dpi = 300, width = 15, height = 10, units = "in",
               plot = DoHeatmap(alldata,
                                features = unique(top$feature), group.by = "sample",
                                size = 4, raster=TRUE, label=FALSE) + 
                 labs(title = paste0(i, " ", levels(as.factor(MNP.genes$group))[1], " vs ",
                                     levels(as.factor(MNP.genes$group))[2])) +
                 theme(axis.text.y = element_text(size = 9),
                       plot.title = element_text(hjust = 0.5, size = 16))) 
        
        #Plotting pathways:
        fgseaResTidy_c$adjPvalue <- ifelse(fgseaResTidy_c$padj <= 0.05, "significant", "non-significant")
        fgseaResTidy_c$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_c$pathway)
        cols <- c("significant" = "red") #"non-significant" = "grey", 
        plot_c <- ggplot(fgseaResTidy_c, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
          geom_col() +
          scale_fill_manual(values = cols) +
          coord_flip() +
          labs(x="Pathway", y="Normalized Enrichment Score",
               title = paste0(i, " ", levels(as.factor(diff.genes$group))[1])) + 
          theme(axis.text.x = element_text(vjust = 0.5, size = 10),
                axis.text.y = element_text(size = 10, vjust = 0.5),
                plot.title = element_text(hjust = 0.5, size = 14),
                axis.title = element_text(size = 10),
                legend.position='right', legend.key.size = unit(0.25, 'cm'),
                legend.key.height = unit(0.25, 'cm'),
                legend.key.width = unit(0.25, 'cm'), 
                legend.title = element_text(size=10),
                legend.text = element_text(size=10)) 
        ggsave(filename = paste0("graphs/gsea-",   i, "-", levels(as.factor(diff.genes$group))[1], "-vs-", levels(as.factor(diff.genes$group))[2], "-", v,  ".tiff"),
               width = 10, height = 10, dpi = 300, units = "in", plot = plot_c)
        
        fgseaResTidy_s$adjPvalue <- ifelse(fgseaResTidy_s$padj <= 0.05, "significant", "non-significant")
        fgseaResTidy_s$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_s$pathway)
        cols <- c("significant" = "red") #"non-significant" = "grey", 
        plot_s <- ggplot(fgseaResTidy_s, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
          geom_col() +
          scale_fill_manual(values = cols) +
          coord_flip() +
          labs(x="Pathway", y="Normalized Enrichment Score",
               title = paste0(i, " ", levels(as.factor(diff.genes$group))[2])) + 
          theme(axis.text.x = element_text(vjust = 0.5, size = 10),
                axis.text.y = element_text(size = 10, vjust = 0.5),
                plot.title = element_text(hjust = 0.5, size = 14),
                axis.title = element_text(size = 10),
                legend.position='none') 
        ggsave(filename = paste0("graphs/gsea-",  i, "-", levels(as.factor(diff.genes$group))[2], "-vs-", levels(as.factor(diff.genes$group))[1],  ".tiff"),
               width = 10, height = 10, dpi = 300, units = "in", plot = plot_s)
        
        obj <- BCR[,BCR$patient == i & BCR$v_gene == v]
        obj$sample <- droplevels(obj$sample)
        DefaultAssay(obj) <- "RNA"
        
        #Running the "quick" Wilcoxon test to extract DEGs:
        diff.genes <- wilcoxauc(obj, "severity", c("critical", "moderate"),
                                seurat_assay='RNA', assay = "data") %>%
          subset(padj < 0.05) %>% arrange(padj)
        
        #Separating moderate genes from critical genes:
        MNP.genes <- diff.genes %>%
          mutate(rank = -log10(pval) * sign(logFC)) %>%
          arrange(-rank)  
        
        #Selecting only the feature and rank columns of DGE data for fgsea run:
        MNP.genes_c <- MNP.genes %>%
          dplyr::filter(group == levels(as.factor(diff.genes$group))[1]) %>%
          arrange(-rank) %>% 
          dplyr::select(feature, rank)
        
        MNP.genes_m <- MNP.genes %>%
          dplyr::filter(group == levels(as.factor(diff.genes$group))[2]) %>%
          arrange(-rank) %>% 
          dplyr::select(feature, rank)
        
        #Converting matrices to tables:
        ranks_c <- deframe(MNP.genes_c)
        ranks_m <- deframe(MNP.genes_m)
        
        #Running fgsea based on gene set ranking (stats = ranks):
        fgseaRes_c <- fgseaMultilevel(fgsea_sets, stats = ranks_c, nPermSimple = 1000, 
                                      maxSize = 200) %>% arrange(padj)
        fgseaRes_m <- fgseaMultilevel(fgsea_sets, stats = ranks_m, nPermSimple = 1000, 
                                      maxSize = 200) %>% arrange(padj)
        
        #Enrichment plots:
        ggsave(filename = paste0("graphs/enrichment-plot-", i, "-", levels(as.factor(MNP.genes$group))[1], "-",
                                 levels(as.factor(MNP.genes$group))[2], "-", v, ".tiff"), 
               dpi = 300, width = 10, height = 5, units = "in",
               plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_c, desc(NES))[1]$pathway]],
                                     ranks_c) + 
                 labs(title = arrange(fgseaRes_c, desc(NES))[1]$pathway) +
                 xlab("Rank") + ylab("Enrichment Score") +
                 theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
                 theme(plot.title = element_text(hjust = 0.5, size = 16)))
        
        ggsave(filename = paste0("graphs/enrichment-plot-", i, "-", levels(as.factor(MNP.genes$group))[2], "-",
                                 levels(as.factor(MNP.genes$group))[1], "-", v, ".tiff"), 
               dpi = 300, width = 10, height = 5, units = "in",
               plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_m, desc(NES))[1]$pathway]],
                                     ranks_m) + 
                 labs(title = arrange(fgseaRes_m, desc(NES))[1]$pathway) +
                 xlab("Rank") + ylab("Enrichment Score") +
                 theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
                 theme(plot.title = element_text(hjust = 0.5, size = 16)))
        
        #Tidying the data:
        fgseaResTidy_c <- fgseaRes_c %>%
          as_tibble() %>%
          arrange(desc(NES)) %>%
          mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ", ", character(1L)))
        
        write.table(fgseaResTidy_c, file = paste0(i, "-critical-moderate-", v, ".tsv"), sep = "\t", col.names = NA)
        
                
        fgseaResTidy_m <- fgseaRes_m %>%
          as_tibble() %>%
          arrange(desc(NES)) %>%
          mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ", ", character(1L)))
        
        write.table(fgseaResTidy_m, file = paste0(i, "-moderate-critical-", v, ".tsv"), sep = "\t", col.names = NA)
        
        
        # choose top upregulated genes
        topUpcrit <- MNP.genes %>% 
          dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
          filter(rank > 0) %>% 
          top_n(50, wt=-padj)
        
        topDowncrit <- MNP.genes %>% 
          dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
          filter(rank < 0) %>% 
          top_n(50, wt=-padj)
        
        topUpmod <- MNP.genes %>% 
          dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
          filter(rank > 0) %>% 
          top_n(60, wt=-padj)
        
        topDownmod <- MNP.genes %>%
          dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
          filter(rank < 0) %>%
          top_n(50, wt=-padj)
        
        
        
        
        
        topPathways <- bind_rows(topUpcrit, topDowncrit, topUpmod, topDownmod)
        
        top <- topPathways %>% group_by(group) 
        
        #Enrichment table graphs:
        tiff(filename = paste0("graphs/gsea-table-", i, "-", v, "-", levels(as.factor(MNP.genes$group))[2], "-", levels(as.factor(MNP.genes$group))[1], ".tiff"), res = 300,
             width = 10, height = 5, units = "in")
        plotGseaTable(fgsea_sets[arrange(fgseaRes_m, desc(NES))[1]$pathway], 
                      ranks_m, 
                      fgseaResTidy_m, 
                      gseaParam = 0.5)
        dev.off()
        
        tiff(filename = paste0("graphs/gsea-table-", i, "-", v, "-", levels(as.factor(MNP.genes$group))[1], "-", levels(as.factor(MNP.genes$group))[2],  ".tiff"), res = 300,
             width = 10, height = 5, units = "in")
        plotGseaTable(fgsea_sets[arrange(fgseaRes_c, desc(NES))[1]$pathway], 
                      ranks_c, 
                      fgseaResTidy_c, 
                      gseaParam = 0.5)
        dev.off()
        
        # create a scale.data slot for the selected genes in subset data
        obj <- obj[,obj$severity == "critical" | obj$severity == "moderate"]
        obj$sample <- droplevels(obj$sample)
        
        alldata <- ScaleData(object = obj, 
                             features = as.character(unique(top$feature), assay = "RNA"))
        
        DefaultAssay(alldata) <- "RNA"
        
        ggsave(filename = paste0("graphs/wilcox-rrho-hm-", i, "-", levels(as.factor(MNP.genes$group))[1], "-",
                                 levels(as.factor(MNP.genes$group))[2], "-", v, ".tiff"), 
               dpi = 300, width = 15, height = 10, units = "in",
               plot = DoHeatmap(alldata,
                                features = unique(top$feature), group.by = "sample",
                                size = 4, raster=TRUE, label=FALSE) + 
                 labs(title = paste0(i, " ", levels(as.factor(MNP.genes$group))[1], " vs ",
                                     levels(as.factor(MNP.genes$group))[2])) +
                 theme(axis.text.y = element_text(size = 9),
                       plot.title = element_text(hjust = 0.5, size = 16)))                  
        
        #Plotting pathways:
        fgseaResTidy_c$adjPvalue <- ifelse(fgseaResTidy_c$padj <= 0.05, "significant", "non-significant")
        fgseaResTidy_c$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_c$pathway)
        cols <- c("significant" = "red") #"non-significant" = "grey", 
        plot_c <- ggplot(fgseaResTidy_c, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
          geom_col() +
          scale_fill_manual(values = cols) +
          coord_flip() +
          labs(x="Pathway", y="Normalized Enrichment Score",
               title = paste0(i, " ", levels(as.factor(diff.genes$group))[1])) + 
          theme(axis.text.x = element_text(vjust = 0.5, size = 10),
                axis.text.y = element_text(size = 10, vjust = 0.5),
                plot.title = element_text(hjust = 0.5, size = 14),
                axis.title = element_text(size = 10),
                legend.position='right', legend.key.size = unit(0.25, 'cm'),
                legend.key.height = unit(0.25, 'cm'),
                legend.key.width = unit(0.25, 'cm'), 
                legend.title = element_text(size=10),
                legend.text = element_text(size=10)) 
        ggsave(filename = paste0("graphs/gsea-",  i, "-", levels(as.factor(diff.genes$group))[1], "-vs-", levels(as.factor(diff.genes$group))[2], "-", v, ".tiff"),
               width = 10, height = 10, dpi = 300, units = "in", plot = plot_c)
        
        fgseaResTidy_m$adjPvalue <- ifelse(fgseaResTidy_m$padj <= 0.05, "significant", "non-significant")
        fgseaResTidy_m$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_m$pathway)
        cols <- c("significant" = "red") #"non-significant" = "grey", 
        plot_m <- ggplot(fgseaResTidy_m, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
          geom_col() +
          scale_fill_manual(values = cols) +
          coord_flip() +
          labs(x="Pathway", y="Normalized Enrichment Score",
               title = paste0(i, " ", levels(as.factor(diff.genes$group))[2])) + 
          theme(axis.text.x = element_text(vjust = 0.5, size = 10),
                axis.text.y = element_text(size = 10, vjust = 0.5),
                plot.title = element_text(hjust = 0.5, size = 14),
                axis.title = element_text(size = 10),
                legend.position='none') 
        ggsave(filename = paste0("graphs/gsea-",  i, "-", levels(as.factor(diff.genes$group))[2], "-vs-", levels(as.factor(diff.genes$group))[1], "-", v, ".tiff"),
               width = 10, height = 10, dpi = 300, units = "in", plot = plot_m)
        
        obj <- BCR[,BCR$patient == i & BCR$v_gene == v]
        obj$sample <- droplevels(obj$sample)
        DefaultAssay(obj) <- "RNA"
        
        
        #Running the "quick" Wilcoxon test to extract DEGs:
        diff.genes <- wilcoxauc(obj, "severity", c("severe", "moderate"),
                                seurat_assay='RNA', assay = "data") %>%
          subset(padj < 0.05) %>% arrange(padj)
        
        #Separating moderate genes from critical genes:
        MNP.genes <- diff.genes %>%
          mutate(rank = -log10(pval) * sign(logFC)) %>%
          arrange(-rank)  
        
        #Selecting only the feature and rank columns of DGE data for fgsea run:
        MNP.genes_s <- MNP.genes %>%
          dplyr::filter(group == levels(as.factor(diff.genes$group))[2]) %>%
          arrange(-rank) %>% 
          dplyr::select(feature, rank)
        
        MNP.genes_m <- MNP.genes %>%
          dplyr::filter(group == levels(as.factor(diff.genes$group))[1]) %>%
          arrange(-rank) %>% 
          dplyr::select(feature, rank)
        
        #Converting matrices to tables:
        ranks_s <- deframe(MNP.genes_s)
        ranks_m <- deframe(MNP.genes_m)
        
        #Running fgsea based on gene set ranking (stats = ranks):
        fgseaRes_s <- fgseaMultilevel(fgsea_sets, stats = ranks_s, nPermSimple = 1000, 
                                      maxSize = 200) %>% arrange(padj)
        fgseaRes_m <- fgseaMultilevel(fgsea_sets, stats = ranks_m, nPermSimple = 1000, 
                                      maxSize = 200) %>% arrange(padj)
        
        #Enrichment plots:
        ggsave(filename = paste0("graphs/enrichment-plot-", i, "-", levels(as.factor(MNP.genes$group))[1], "-",
                                 levels(as.factor(MNP.genes$group))[2], "-", v, ".tiff"), 
               dpi = 300, width = 10, height = 5, units = "in",
               plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_m, desc(NES))[1]$pathway]],
                                     ranks_m) + 
                 labs(title = arrange(fgseaRes_m, desc(NES))[1]$pathway) +
                 xlab("Rank") + ylab("Enrichment Score") +
                 theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
                 theme(plot.title = element_text(hjust = 0.5, size = 16)))
        
        ggsave(filename = paste0("graphs/enrichment-plot-", i, "-", levels(as.factor(MNP.genes$group))[2], "-",
                                 levels(as.factor(MNP.genes$group))[1], "-", v, ".tiff"), 
               dpi = 300, width = 10, height = 5, units = "in",
               plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_s, desc(NES))[1]$pathway]],
                                     ranks_s) + 
                 labs(title = arrange(fgseaRes_s, desc(NES))[1]$pathway) +
                 xlab("Rank") + ylab("Enrichment Score") +
                 theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
                 theme(plot.title = element_text(hjust = 0.5, size = 16)))
        
        #Tidying the data:
        fgseaResTidy_s <- fgseaRes_s %>%
          as_tibble() %>%
          arrange(desc(NES)) %>%
          mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ", ", character(1L)))
        
        write.table(fgseaResTidy_s, file = paste0(i, "-severe-moderate-", v, ".tsv"), sep = "\t", col.names = NA)
        
                
        fgseaResTidy_m <- fgseaRes_m %>%
          as_tibble() %>%
          arrange(desc(NES)) %>%
          mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ", ", character(1L)))
        
        write.table(fgseaResTidy_m, file = paste0(i, "-moderate-severe-", v, ".tsv"), sep = "\t", col.names = NA)
        
        # choose top upregulated genes
        topUpmod <- MNP.genes %>% 
          dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
          filter(rank > 0) %>% 
          top_n(60, wt=-padj)
        
        topDownmod <- MNP.genes %>%
          dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
          filter(rank < 0) %>%
          top_n(50, wt=-padj)
        
        topUpsev <- MNP.genes %>% 
          dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
          filter(rank > 0) %>% 
          top_n(50, wt=-padj)
        
        topDownsev <- MNP.genes %>% 
          dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
          filter(rank < 0) %>% 
          top_n(50, wt=-padj)
        
        
        
        topPathways <- bind_rows(topUpsev, topDownsev, topUpmod, topDownmod)
        
        top <- topPathways %>% group_by(group) 
        
        #Enrichment table graphs:
        tiff(filename = paste0("graphs/gsea-table-", i, "-", v, "-", levels(as.factor(MNP.genes$group))[1], "-", levels(as.factor(MNP.genes$group))[2], ".tiff"), res = 300,
             width = 10, height = 5, units = "in")
        plotGseaTable(fgsea_sets[arrange(fgseaRes_m, desc(NES))[1]$pathway], 
                      ranks_m, 
                      fgseaResTidy_m, 
                      gseaParam = 0.5)
        dev.off()
        
        tiff(filename = paste0("graphs/gsea-table-", i, "-", v, "-", levels(as.factor(MNP.genes$group))[2], "-", levels(as.factor(MNP.genes$group))[1], ".tiff"), res = 300,
             width = 10, height = 5, units = "in")
        plotGseaTable(fgsea_sets[arrange(fgseaRes_s, desc(NES))[1]$pathway], 
                      ranks_s, 
                      fgseaResTidy_s, 
                      gseaParam = 0.5)
        dev.off()
        
        # create a scale.data slot for the selected genes in subset data
        obj <- obj[,obj$severity == "severe" | obj$severity == "moderate"]
        obj$sample <- droplevels(obj$sample)
        
        alldata <- ScaleData(object = obj, 
                             features = as.character(unique(top$feature), assay = "RNA"))
        
        DefaultAssay(alldata) <- "RNA"
        
        ggsave(filename = paste0("graphs/wilcox-rrho-hm-", i, "-", levels(as.factor(MNP.genes$group))[1], "-",
                                 levels(as.factor(MNP.genes$group))[2], "-", v, ".tiff"), 
               dpi = 300, width = 15, height = 10, units = "in",
               plot = DoHeatmap(alldata,
                                features = unique(top$feature), group.by = "sample",
                                size = 4, raster=TRUE, label=FALSE) + 
                 labs(title = paste0(i, " ", levels(as.factor(MNP.genes$group))[1], " vs ",
                                     levels(as.factor(MNP.genes$group))[2])) +
                 theme(axis.text.y = element_text(size = 9),
                       plot.title = element_text(hjust = 0.5, size = 16)))
        
        #Plotting pathways:
        fgseaResTidy_s$adjPvalue <- ifelse(fgseaResTidy_s$padj <= 0.05, "significant", "non-significant")
        fgseaResTidy_s$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_s$pathway)
        cols <- c("significant" = "red") #"non-significant" = "grey", 
        plot_s <- ggplot(fgseaResTidy_s, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
          geom_col() +
          scale_fill_manual(values = cols) +
          coord_flip() +
          labs(x="Pathway", y="Normalized Enrichment Score",
               title = paste0(i, " ", "severe")) + 
          theme(axis.text.x = element_text(vjust = 0.5, size = 10),
                axis.text.y = element_text(size = 10, vjust = 0.5),
                plot.title = element_text(hjust = 0.5, size = 14),
                axis.title = element_text(size = 10),
                legend.position='none') 
        ggsave(filename = paste0("graphs/gsea-",  i, "-", "severe-vs-moderate", "-", v, ".tiff"),
               width = 10, height = 10, dpi = 300, units = "in", plot = plot_s)
        
        fgseaResTidy_m$adjPvalue <- ifelse(fgseaResTidy_m$padj <= 0.05, "significant", "non-significant")
        fgseaResTidy_m$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_m$pathway)
        cols <- c("significant" = "red") #"non-significant" = "grey", 
        plot_m <- ggplot(fgseaResTidy_m, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
          geom_col() +
          scale_fill_manual(values = cols) +
          coord_flip() +
          labs(x="Pathway", y="Normalized Enrichment Score",
               title = paste0(i, " ", "moderate")) + 
          theme(axis.text.x = element_text(vjust = 0.5, size = 10),
                axis.text.y = element_text(size = 10, vjust = 0.5),
                plot.title = element_text(hjust = 0.5, size = 14),
                axis.title = element_text(size = 10),
                legend.position='none') 
        ggsave(filename = paste0("graphs/gsea-",  i, "-", "moderate-vs-severe", "-", v, ".tiff"),
               width = 10, height = 10, dpi = 300, units = "in", plot = plot_m)
        
        
      } else {
        #Running the "quick" Wilcoxon test to extract DEGs:
        diff.genes <- wilcoxauc(obj, "severity", 
                                seurat_assay='RNA', assay = "data") %>%
          subset(padj < 0.05) %>% arrange(padj)
        
        #Separating moderate genes from critical genes:
        MNP.genes <- diff.genes %>%
          mutate(rank = -log10(pval) * sign(logFC)) %>%
          arrange(-rank)  
        
        #Selecting only the feature and rank columns of DGE data for fgsea run:
        MNP.genes_c <- MNP.genes %>%
          dplyr::filter(group == levels(as.factor(diff.genes$group))[1]) %>%
          arrange(-rank) %>% 
          dplyr::select(feature, rank)
        
        MNP.genes_m <- MNP.genes %>%
          dplyr::filter(group == levels(as.factor(diff.genes$group))[2]) %>%
          arrange(-rank) %>% 
          dplyr::select(feature, rank)
        
        #Converting matrices to tables:
        ranks_c <- deframe(MNP.genes_c)
        ranks_m <- deframe(MNP.genes_m)
        
        #Running fgsea based on gene set ranking (stats = ranks):
        fgseaRes_c <- fgseaMultilevel(fgsea_sets, stats = ranks_c, nPermSimple = 1000, 
                                      maxSize = 200) %>% arrange(padj)
        fgseaRes_m <- fgseaMultilevel(fgsea_sets, stats = ranks_m, nPermSimple = 1000, 
                                      maxSize = 200) %>% arrange(padj)
        
        #Enrichment plots:
        ggsave(filename = paste0("graphs/enrichment-plot-", i, "-", levels(as.factor(MNP.genes$group))[1], "-",
                                 levels(as.factor(MNP.genes$group))[2], "-", v, ".tiff"), 
               dpi = 300, width = 10, height = 5, units = "in",
               plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_c, desc(NES))[1]$pathway]],
                                     ranks_c) + 
                 labs(title = arrange(fgseaRes_c, desc(NES))[1]$pathway) +
                 xlab("Rank") + ylab("Enrichment Score") +
                 theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
                 theme(plot.title = element_text(hjust = 0.5, size = 16)))
        
        ggsave(filename = paste0("graphs/enrichment-plot-", i, "-", levels(as.factor(MNP.genes$group))[2], "-",
                                 levels(as.factor(MNP.genes$group))[1], "-", v, ".tiff"), 
               dpi = 300, width = 10, height = 5, units = "in",
               plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_m, desc(NES))[1]$pathway]],
                                     ranks_m) + 
                 labs(title = arrange(fgseaRes_m, desc(NES))[1]$pathway) +
                 xlab("Rank") + ylab("Enrichment Score") +
                 theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
                 theme(plot.title = element_text(hjust = 0.5, size = 16)))
        
        #Tidying the data:
        fgseaResTidy_m <- fgseaRes_m %>%
          as_tibble() %>%
          arrange(desc(NES)) %>%
          mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ", ", character(1L)))
        
        write.table(fgseaResTidy_m, file = paste0(i, "-", levels(as.factor(diff.genes$group))[2], "-critical-", v, ".tsv"), sep = "\t", col.names = NA)
        
                
        fgseaResTidy_c <- fgseaRes_c %>%
          as_tibble() %>%
          arrange(desc(NES)) %>%
          mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ", ", character(1L)))
        
        write.table(fgseaResTidy_c, file = paste0(i, "-critical-",  v, levels(as.factor(diff.genes$group))[2],".tsv"), sep = "\t", col.names = NA)
        
                
        # choose top upregulated genes
        topUpcrit <- MNP.genes %>% 
          dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
          filter(rank > 0) %>% 
          top_n(50, wt=-padj)
        
        topDowncrit <- MNP.genes %>% 
          dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
          filter(rank < 0) %>% 
          top_n(50, wt=-padj)
        
        topUpmod <- MNP.genes %>% 
          dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
          filter(rank > 0) %>% 
          top_n(60, wt=-padj)
        
        topDownmod <- MNP.genes %>%
          dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
          filter(rank < 0) %>%
          top_n(50, wt=-padj)
        
        topPathways <- bind_rows(topUpcrit, topDowncrit, topUpmod, topDownmod)
        
        top <- topPathways %>% group_by(group) 
        
        #Enrichment table graphs:
        tiff(filename = paste0("graphs/gsea-table-", i, "-", v, "-", levels(as.factor(MNP.genes$group))[2], "-", levels(as.factor(MNP.genes$group))[1], ".tiff"), res = 300,
             width = 10, height = 5, units = "in")
        plotGseaTable(fgsea_sets[arrange(fgseaRes_m, desc(NES))[1]$pathway], 
                      ranks_m, 
                      fgseaResTidy_m, 
                      gseaParam = 0.5)
        dev.off()
        
        tiff(filename = paste0("graphs/gsea-table-", i, "-", v, "-", levels(as.factor(MNP.genes$group))[1], "-", levels(as.factor(MNP.genes$group))[2], ".tiff"), res = 300,
             width = 10, height = 5, units = "in")
        plotGseaTable(fgsea_sets[arrange(fgseaRes_c, desc(NES))[1]$pathway], 
                      ranks_c, 
                      fgseaResTidy_c, 
                      gseaParam = 0.5)
        dev.off()
        
        
        # create a scale.data slot for the selected genes in subset data
        alldata <- ScaleData(object = obj[,obj$severity == levels(as.factor(MNP.genes$group))[1] | obj$severity == levels(as.factor(MNP.genes$group))[2]], 
                             features = as.character(unique(top$feature), assay = "RNA"))
        
        DefaultAssay(alldata) <- "RNA"
        
        ggsave(filename = paste0("graphs/wilcox-rrho-hm-", i, "-", v, ".tiff"), dpi = 300,
               width = 15, height = 10, units = "in",
               plot = DoHeatmap(alldata,
                                features = unique(top$feature), group.by = "sample",
                                size = 4, raster=TRUE, label=FALSE) + 
                 labs(title = i) +
                 theme(axis.text.y = element_text(size = 9),
                       plot.title = element_text(hjust = 0.5, size = 16)))
        
        #Plotting pathways:
        fgseaResTidy_c$adjPvalue <- ifelse(fgseaResTidy_c$padj <= 0.05, "significant", "non-significant")
        fgseaResTidy_c$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_c$pathway)
        cols <- c("significant" = "red") #"non-significant" = "grey", 
        plot_c <- ggplot(fgseaResTidy_c, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
          geom_col() +
          scale_fill_manual(values = cols) +
          coord_flip() +
          labs(x="Pathway", y="Normalized Enrichment Score",
               title = paste0(i, " ", levels(as.factor(diff.genes$group))[1])) + 
          theme(axis.text.x = element_text(vjust = 0.5, size = 10),
                axis.text.y = element_text(size = 10, vjust = 0.5),
                plot.title = element_text(hjust = 0.5, size = 14),
                axis.title = element_text(size = 10),
                legend.position='right', legend.key.size = unit(0.25, 'cm'),
                legend.key.height = unit(0.25, 'cm'),
                legend.key.width = unit(0.25, 'cm'), 
                legend.title = element_text(size=10),
                legend.text = element_text(size=10)) 
        ggsave(filename = paste0("graphs/gsea-",  i, "-", levels(as.factor(diff.genes$group))[1], "-", v, ".tiff"),
               width = 10, height = 10, dpi = 300, units = "in", plot = plot_c)
        
        fgseaResTidy_m$adjPvalue <- ifelse(fgseaResTidy_m$padj <= 0.05, "significant", "non-significant")
        fgseaResTidy_m$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_m$pathway)
        cols <- c("significant" = "red") #"non-significant" = "grey", 
        plot_m <- ggplot(fgseaResTidy_m, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
          geom_col() +
          scale_fill_manual(values = cols) +
          coord_flip() +
          labs(x="Pathway", y="Normalized Enrichment Score",
               title = paste0(i, " ", levels(as.factor(diff.genes$group))[2])) + 
          theme(axis.text.x = element_text(vjust = 0.5, size = 10),
                axis.text.y = element_text(size = 10, vjust = 0.5),
                plot.title = element_text(hjust = 0.5, size = 14),
                axis.title = element_text(size = 10),
                legend.position='none') 
        ggsave(filename = paste0("graphs/gsea-",  i, "-", levels(as.factor(diff.genes$group))[2], "-", v, ".tiff"),
               width = 10, height = 10, dpi = 300, units = "in", plot = plot_m)
        
      }
      
    }

  )
  
}
}


################################For all B cells#################################
#Running the "quick" Wilcoxon test to extract DEGs:
diff.genes <- wilcoxauc(BCR, "severity", c("critical", "moderate"),
                        seurat_assay='RNA', assay = "data") %>%
  subset(padj < 0.05) %>% arrange(padj)

#Separating moderate genes from critical genes:
MNP.genes <- diff.genes %>%
  mutate(rank = -log10(pval) * sign(logFC)) %>%
  arrange(-rank)  

#Selecting only the feature and rank columns of DGE data for fgsea run:
MNP.genes_c <- MNP.genes %>%
  dplyr::filter(group == levels(as.factor(diff.genes$group))[1]) %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)

MNP.genes_m <- MNP.genes %>%
  dplyr::filter(group == levels(as.factor(diff.genes$group))[2]) %>%
  arrange(-rank) %>% 
  dplyr::select(feature, rank)

#Converting matrices to tables:
ranks_c <- deframe(MNP.genes_c)
ranks_m <- deframe(MNP.genes_m)

#Running fgsea based on gene set ranking (stats = ranks):
fgseaRes_c <- fgseaMultilevel(fgsea_sets, stats = ranks_c, nPermSimple = 1000, 
                              maxSize = 200) %>% arrange(padj)
fgseaRes_m <- fgseaMultilevel(fgsea_sets, stats = ranks_m, nPermSimple = 1000, 
                              maxSize = 200) %>% arrange(padj)

#Enrichment plots:
ggsave(filename = paste0("graphs/enrichment-plot-all-", levels(as.factor(MNP.genes$group))[1], "-",
                         levels(as.factor(MNP.genes$group))[2], "-", v, ".tiff"), 
       dpi = 300, width = 15, height = 10, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_c, desc(NES))[1]$pathway]],
                             ranks_c) + 
         labs(title = arrange(fgseaRes_c, desc(NES))[1]$pathway) + #Suggested title improvement: 
         xlab("Rank") + ylab("Enrichment Score") +                 #arrange(fgseaRes_c, desc(NES))[1]$pathway %>% sub(pattern = "HALLMARK_", replacement = "") %>% sub(pattern = "_", replacement = " ")
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))

ggsave(filename = paste0("graphs/enrichment-plot-all-", levels(as.factor(MNP.genes$group))[2], "-",
                         levels(as.factor(MNP.genes$group))[1], "-", v, ".tiff"), 
       dpi = 300, width = 15, height = 10, units = "in",
       plot = plotEnrichment(fgsea_sets[[arrange(fgseaRes_m, desc(NES))[1]$pathway]],
                             ranks_m) + 
         labs(title = arrange(fgseaRes_m, desc(NES))[1]$pathway) +
         xlab("Rank") + ylab("Enrichment Score") +
         theme(text = element_text(size = 12), axis.text = element_text(size = 12)) +
         theme(plot.title = element_text(hjust = 0.5, size = 16)))

#Tidying the data:
fgseaResTidy_m <- fgseaRes_m %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ", ", character(1L)))

write.table(fgseaResTidy_m, file = "all-b-cells-moderate-vs-critical.tsv", sep = "\t", col.names = NA)


fgseaResTidy_c <- fgseaRes_c %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ", ", character(1L)))

write.table(fgseaResTidy_c, file = "all-b-cells-critical-vs-moderate.tsv", sep = "\t", col.names = NA)


# choose top upregulated genes
topUpcrit <- MNP.genes %>% 
  dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
  filter(rank > 0) %>% 
  top_n(50, wt=-padj)

topDowncrit <- MNP.genes %>% 
  dplyr::filter(group == levels(as.factor(MNP.genes$group))[1]) %>%
  filter(rank < 0) %>% 
  top_n(50, wt=-padj)

topUpmod <- MNP.genes %>% 
  dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
  filter(rank > 0) %>% 
  top_n(60, wt=-padj)

topDownmod <- MNP.genes %>%
  dplyr::filter(group == levels(as.factor(MNP.genes$group))[2]) %>%
  filter(rank < 0) %>%
  top_n(50, wt=-padj)

topPathways <- bind_rows(topUpcrit, topDowncrit, topUpmod, topDownmod)

top <- topPathways %>% group_by(group) 

#Enrichment table graphs:
plotGseaTable(fgsea_sets[arrange(fgseaRes_m, desc(NES))[1]$pathway], 
              ranks_m, 
              fgseaResTidy_m, 
              gseaParam = 0.5)

plotGseaTable(fgsea_sets[arrange(fgseaRes_c, desc(NES))[1]$pathway], 
              ranks_c, 
              fgseaResTidy_c, 
              gseaParam = 0.5)

# create a scale.data slot for the selected genes in subset data
alldata <- ScaleData(object = BCR, 
                     features = as.character(unique(top$feature), assay = "RNA"))

DefaultAssay(alldata) <- "RNA"

ggsave(filename = paste0("graphs/wilcox-rrho-hm-", i, "-", v, ".tiff"), dpi = 300,
       width = 15, height = 10, units = "in",
       plot = DoHeatmap(alldata,
                        features = unique(top$feature), group.by = "sample",
                        size = 4, raster=TRUE, label=FALSE) + 
         labs(title = i) +
         theme(axis.text.y = element_text(size = 9),
               plot.title = element_text(hjust = 0.5, size = 16)))

#Plotting pathways:
fgseaResTidy_c$adjPvalue <- ifelse(fgseaResTidy_c$padj <= 0.05, "significant", "non-significant")
fgseaResTidy_c$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_c$pathway)
cols <- c("significant" = "red") #"non-significant" = "grey", 
plot_c <- ggplot(fgseaResTidy_c, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = paste0(i, " ", levels(as.factor(diff.genes$group))[1])) + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='right', legend.key.size = unit(0.25, 'cm'),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) 
ggsave(filename = paste0("graphs/gsea-",  i, "-", levels(as.factor(diff.genes$group))[1], "-", v, ".tiff"),
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_c)

fgseaResTidy_m$adjPvalue <- ifelse(fgseaResTidy_m$padj <= 0.05, "significant", "non-significant")
fgseaResTidy_m$pathway <- gsub(pattern = "HALLMARK_", replacement = "", x = fgseaResTidy_m$pathway)
cols <- c("significant" = "red") #"non-significant" = "grey", 
plot_m <- ggplot(fgseaResTidy_m, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title = paste0(i, " ", levels(as.factor(diff.genes$group))[2])) + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 10),
        legend.position='none') 
ggsave(filename = paste0("graphs/gsea-",  i, "-", levels(as.factor(diff.genes$group))[2], "-", v, ".tiff"),
       width = 10, height = 10, dpi = 300, units = "in", plot = plot_m)
