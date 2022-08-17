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


#Choosing the human Hallmark annotated gene set from the MsigDB:
m_df <- msigdbr(species = "Homo sapiens", category = "H")


#Preparing a list of the gene sets, for fgsea:
fgsea_sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#Printing out the most abundant V-J combinations for each sample:
hiclo <- NULL
for (i in levels(factor(BCR@meta.data$patient[BCR$severity != "healthy"]))) {
  BCR@meta.data[BCR$patient == i & !is.na(BCR$v_gene),] %>%
    group_by(v_gene, patient)  %>% dplyr::count() %>% na.omit() %>%
    arrange(desc(n)) %>% as.data.frame() -> matrix
  hiclo <- rbind(hiclo, matrix[1,])
  
}

rm(list = setdiff(ls(), c("BCR", "hiclo", "m_df", "fgsea_sets")))

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


################################FGSEA###########################################





################################Normalized Enrichment Score Plots###############























