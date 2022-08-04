####################################Loading#####################################
library(Seurat)
library(tidyverse)
library(Matrix)
library(ggplot2)

covid <- readRDS("04-covid-clustered.rds")
BCR <- readRDS("05-BCR-combined.rds")
BCR <- BCR[,BCR$patient != "Control 4"]
BCR@meta.data$sample <- droplevels(BCR@meta.data$sample)

####################################UMAPs#######################################
#UMAP plots for control samples:
ggsave(filename = "azimuth-healthy-all.jpeg", path = "./graphs/", 
       dpi = "print", width = 10, height = 5,
       plot = DimPlot(covid[,covid$severity == "healthy"], reduction = "umap",
                      repel = T) + ggtitle("Controls") +
         theme(plot.title = element_text(hjust = 0.5)))
ggsave(filename = "azimuth-healthy-grid.tiff", path = "./graphs/", 
       dpi = "print", width = 15, height = 5,
       plot = DimPlot(covid[,covid$severity == "healthy" & covid$patient != "Control 4"], reduction = "umap",
                      repel = T, split.by = "patient") + NoLegend())
DimPlot(covid[,covid$sample == i], reduction = "umap",
        repel = T, label = T) + NoLegend()

#UMAP plots for individual patients:
#Patient1:
pt1 <- covid[,covid$patient == "Patient1"]
pt1$severity <- factor(pt1$severity, levels = c("moderate", "critical"))

ggsave(filename = "azimuth-pt1.tiff", path = "./graphs/", 
       dpi = "print", width = 10, height = 5,
       plot = DimPlot(pt1, split.by = "severity",
                      reduction = "umap", repel = T) + 
         NoLegend() +
         ggtitle("Patient 1") + theme(plot.title = element_text(hjust = 0.5)))

#Patient2:
pt2 <- covid[,covid$patient == "Patient2"]
pt2$severity <- factor(pt2$severity, levels = c("moderate", "critical"))

ggsave(filename = "azimuth-pt2.jpeg", path = "./graphs/", 
       dpi = "print", width = 10, height = 5,
       plot = DimPlot(pt2, split.by = "severity",
                      reduction = "umap", repel = T) + 
         NoLegend() +
         ggtitle("Patient 2") + theme(plot.title = element_text(hjust = 0.5)))

#Patient3:
pt3 <- covid[,covid$patient == "Patient3"]
pt3$severity <- factor(pt3$severity, levels = c("mild", "critical"))

ggsave(filename = "azimuth-pt3.jpeg", path = "./graphs/", 
       dpi = "print", width = 10, height = 5,
       plot = DimPlot(pt3, split.by = "severity",
                      reduction = "umap", repel = T) + 
         NoLegend() +
         ggtitle("Patient 3") + theme(plot.title = element_text(hjust = 0.5)))

#Patient4:
pt4 <- covid[,covid$patient == "Patient4"]
pt4$severity <- factor(pt4$severity, levels = c("mild", "critical"))

ggsave(filename = "azimuth-pt4.jpeg", path = "./graphs/", 
       dpi = "print", width = 10, height = 5,
       plot = DimPlot(pt4, split.by = "severity",
                      reduction = "umap", repel = T) + 
         NoLegend() +
         ggtitle("Patient 4") + theme(plot.title = element_text(hjust = 0.5)))

#Patient5:
pt5 <- covid[,covid$patient == "Patient5"]
pt5$severity <- factor(pt5$severity, levels = c("critical", "severe", "moderate"))

ggsave(filename = "azimuth-pt5.jpeg", path = "./graphs/", 
       dpi = "print", width = 10, height = 5,
       plot = DimPlot(pt5, split.by = "severity",
                      reduction = "umap", repel = T) + 
         NoLegend() +
         ggtitle("Patient 5") + theme(plot.title = element_text(hjust = 0.5)))

#Patient6:
pt6 <- covid[,covid$patient == "Patient6"]
pt6$severity <- factor(pt6$severity, levels = c("critical", "severe", "moderate"))

ggsave(filename = "azimuth-pt6.jpeg", path = "./graphs/", 
       dpi = "print", width = 10, height = 5,
       plot = DimPlot(pt6, split.by = "severity",
                      reduction = "umap", repel = T) + 
         NoLegend() +
         ggtitle("Patient 6") + theme(plot.title = element_text(hjust = 0.5)))



for (i in levels(as.factor(covid[,covid$severity == "healthy"]$patient))) {
  
  ggsave(filename = paste0("azimuth-healthy-", i,".jpeg" ), path = "./graphs/", 
         dpi = "print", width = 10, height = 10,
         plot = DimPlot(covid[,covid$patient == i], reduction = "umap",
                        label = T, repel = T) + NoLegend() + ggtitle(i) + theme(plot.title = element_text(hjust = 0.5)))
  
}



####################################BCR#########################################
#B cell distribution bar graphs:
for (i in levels(factor(BCR$sample))) {
  
  ggsave(filename = paste0(i, "-b-cells-distribution.jpeg"), path = "./graphs/", dpi = "print",
         plot = ggplot(data = BCR@meta.data[BCR$sample == i,], aes(x = severity, fill = azimuthNames)) +
           geom_histogram(stat = "count", position = "fill") +
           scale_y_continuous(name = "",
                              breaks = c(0, 0.25, 0.5, 0.75, 1),
                              labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) +
           xlab("") + guides(x = "none"))
  
  
}

#Histograms for V gene usage in all:
BCR@meta.data %>%
        as.data.frame() -> df
df$v_gene <- vapply(strsplit(df$v_gene, "[V]"), "[", "", 2)
df$j_gene <- vapply(strsplit(df$j_gene, "[J]"), "[", "", 2)

      
plot <- ggplot(data = df, aes(x = v_gene, fill = j_gene)) +
geom_bar() + 
xlab("V Gene") + ylab("Count") + ggtitle("All") +
  theme(axis.text.x = element_text(angle = 90, size = 24),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20))+ 
  scale_fill_discrete(name = "J Gene", type = c("#FFA700", "#E206CB", "#DC143C",
                                                "#4666FF", "#EEE600", "#39FF14"))

ggsave(filename = paste0("v-j-gene-bar-all-stacked.jpeg"), plot = plot,
       dpi = "print", height = 10, width = 20, path = "./graphs/")

#Histograms for V gene usage in healthy controls:
BCR@meta.data[BCR$severity == "healthy" & !is.na(BCR$v_gene),] %>%
        as.data.frame() -> df
df$v_gene <- vapply(strsplit(df$v_gene, "[V]"), "[", "", 2)
df$j_gene <- vapply(strsplit(df$j_gene, "[J]"), "[", "", 2)
df$patient <- factor(df$patient, levels = c("Control 1", "Control 2", "Control 3"))
      
plot <- ggplot(data = df, aes(x = v_gene, fill = j_gene)) +
geom_bar() + facet_wrap(~patient, ncol = 1) +
xlab("V Gene") + ylab("Count") + 
  theme(axis.text.x = element_text(angle = 90, size = 24),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20)) + 
        scale_fill_discrete(name = "J Gene", type = c("#FFA700", "#E206CB", "#DC143C",
                                                      "#4666FF", "#EEE600", "#39FF14"))

ggsave(filename = paste0("v-j-gene-bar-hc-stacked.jpeg"), plot = plot,
       dpi = "print", height = 10, width = 20, path = "./graphs/")

#Histograms for V gene usage in Patient 1:
BCR@meta.data[BCR$patient == "Patient 1",] %>%
      na.omit() %>% 
      as.data.frame() -> df
df$v_gene <- vapply(strsplit(df$v_gene, "[V]"), "[", "", 2)
df$j_gene <- vapply(strsplit(df$j_gene, "[J]"), "[", "", 2)
df$severity <- factor(df$severity, levels = c("moderate", "critical"))
      
plot <- ggplot(data = df, aes(x = v_gene, fill = j_gene)) +
geom_bar() + facet_wrap(~severity, ncol = 1) +
xlab("V Gene") + ylab("Count") + ggtitle(df$patient) +
  theme(axis.text.x = element_text(angle = 90, size = 24),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20)) + 
  scale_fill_discrete(name = "J Gene", type = c("#FFA700", "#E206CB", "#DC143C",
                                                "#4666FF", "#EEE600", "#39FF14"))
    
ggsave(filename = paste0("v-j-gene-bar-pt1-stacked.jpeg"), plot = plot,
             dpi = "print", height = 10, width = 20, path = "./graphs/")
      
#Histograms for V gene usage in Patient 2:
BCR@meta.data[BCR$patient == "Patient 2",] %>%
      na.omit() %>% 
      as.data.frame() -> df
df$v_gene <- vapply(strsplit(df$v_gene, "[V]"), "[", "", 2)
df$j_gene <- vapply(strsplit(df$j_gene, "[J]"), "[", "", 2)
df$severity <- factor(df$severity, levels = c("moderate", "critical"))
      
plot <- ggplot(data = df, aes(x = v_gene, fill = j_gene)) +
geom_bar() + facet_wrap(~severity, ncol = 1) +
xlab("V Gene") + ylab("Count") + ggtitle(df$patient) +
  theme(axis.text.x = element_text(angle = 90, size = 24),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20)) + 
  scale_fill_discrete(name = "J Gene", type = c("#FFA700", "#E206CB", "#DC143C",
                                                "#4666FF", "#EEE600", "#39FF14")) 
    
ggsave(filename = paste0("v-j-gene-bar-pt2-stacked.jpeg"), plot = plot,
             dpi = "print", height = 10, width = 20, path = "./graphs/")
      
#Histograms for V gene usage in Patient 3:
BCR@meta.data[BCR$patient == "Patient 3",] %>%
      na.omit() %>% 
      as.data.frame() -> df
df$v_gene <- vapply(strsplit(df$v_gene, "[V]"), "[", "", 2)
df$j_gene <- vapply(strsplit(df$j_gene, "[J]"), "[", "", 2)
df$severity <- factor(df$severity, levels = c("mild", "critical"))
      
plot <- ggplot(data = df, aes(x = v_gene, fill = j_gene)) +
geom_bar() + facet_wrap(~severity, ncol = 1) +
xlab("V Gene") + ylab("Count") + ggtitle(df$patient) +
  theme(axis.text.x = element_text(angle = 90, size = 24),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20)) + 
  scale_fill_discrete(name = "J Gene", type = c("#FFA700", "#E206CB", "#DC143C",
                                                "#4666FF", "#EEE600", "#39FF14"))
    
ggsave(filename = paste0("v-j-gene-bar-pt3-stacked.jpeg"), plot = plot,
             dpi = "print", height = 10, width = 20, path = "./graphs/")
      
#Histograms for V gene usage in Patient 4:
BCR@meta.data[BCR$patient == "Patient 4",] %>%
      na.omit() %>% 
      as.data.frame() -> df
df$v_gene <- vapply(strsplit(df$v_gene, "[V]"), "[", "", 2)
df$j_gene <- vapply(strsplit(df$j_gene, "[J]"), "[", "", 2)
df$severity <- factor(df$severity, levels = c("mild", "critical"))
      
plot <- ggplot(data = df, aes(x = v_gene, fill = j_gene)) +
geom_bar() + facet_wrap(~severity, ncol = 1) +
xlab("V Gene") + ylab("Count") + ggtitle(df$patient) +
  theme(axis.text.x = element_text(angle = 90, size = 24),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20)) + 
  scale_fill_discrete(name = "J Gene", type = c("#FFA700", "#E206CB", "#DC143C",
                                                "#4666FF", "#EEE600", "#39FF14")) 
    
ggsave(filename = paste0("v-j-gene-bar-pt4-stacked.jpeg"), plot = plot,
             dpi = "print", height = 10, width = 20, path = "./graphs/")
      
#Histograms for V gene usage in Patient 5:
BCR@meta.data[BCR$patient == "Patient 5",] %>%
      na.omit() %>% 
      as.data.frame() -> df
df$v_gene <- vapply(strsplit(df$v_gene, "[V]"), "[", "", 2)
df$j_gene <- vapply(strsplit(df$j_gene, "[J]"), "[", "", 2)
df$severity <- factor(df$severity, levels = c("critical", "severe", "moderate"))
      
plot <- ggplot(data = df, aes(x = v_gene, fill = j_gene)) +
geom_bar() + facet_wrap(~severity, ncol = 1) +
xlab("V Gene") + ylab("Count") + ggtitle(df$patient) +
  theme(axis.text.x = element_text(angle = 90, size = 24),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20)) + 
  scale_fill_discrete(name = "J Gene", type = c("#FFA700", "#E206CB", "#DC143C",
                                                "#4666FF", "#EEE600", "#39FF14"))
    
ggsave(filename = paste0("v-j-gene-bar-pt5-stacked.jpeg"), plot = plot,
             dpi = "print", height = 15, width = 20, path = "./graphs/")
      
#Histograms for V gene usage in Patient 6:
BCR@meta.data[BCR$patient == "Patient 6",] %>%
      na.omit() %>% 
      as.data.frame() -> df
df$v_gene <- vapply(strsplit(df$v_gene, "[V]"), "[", "", 2)
df$j_gene <- vapply(strsplit(df$j_gene, "[J]"), "[", "", 2)
df$severity <- factor(df$severity, levels = c("critical", "severe", "moderate"))
      
plot <- ggplot(data = df, aes(x = v_gene, fill = j_gene)) +
geom_bar() + facet_wrap(~severity, ncol = 1) +
xlab("V Gene") + ylab("Count") + ggtitle(df$patient) +
theme(axis.text.x = element_text(angle = 90, size = 24),
      axis.text.y = element_text(size = 15),
      axis.title = element_text(size = 20),
      strip.text = element_text(size = 20)) + 
  scale_fill_discrete(name = "J Gene", type = c("#FFA700", "#E206CB", "#DC143C",
                                                "#4666FF", "#EEE600", "#39FF14"))
    
ggsave(filename = paste0("v-j-gene-bar-pt6-stacked.jpeg"), plot = plot,
             dpi = "print", height = 15, width = 20, path = "./graphs/")
      
#Heatmaps for V and J gene usage in B cells:

for (i in levels(as.factor(BCR$sample))) {
  for (j in levels(as.factor(BCR@meta.data[BCR$sample == i,]$azimuthNames))) {
    try({
      BCR@meta.data[BCR$sample == i & BCR$azimuthNames == j,] %>%
        group_by(v_gene, j_gene) %>% na.omit() %>% dplyr::count() %>% spread(key = j_gene, value = n) %>%
        as.data.frame() -> matrix
      rownames(matrix) <- matrix$v_gene
      matrix$v_gene <- NULL
      matrix <- as.matrix(matrix)
      matrix[is.na(matrix)] <- 0
      plot <- pheatmap(matrix, cluster_cols = F, cluster_rows = F, height = 26, width = 10, 
                       labels_row = vapply(strsplit(rownames(matrix), "[V]"), "[", "", 2),
                       labels_col = vapply(strsplit(colnames(matrix), "[J]"), "[", "", 2),
                       cellwidth = 45, cellheight = 45, angle_col = 0,
                       fontsize_col = 60, fontsize_row = 60, fontsize = 28)
      ggsave(filename = paste0("v-j-heatmap-", i, "-", j, ".jpeg"), plot = plot,
             dpi = "print", height = 26, width = 10, path = "./graphs/")
      
    })
  }
}

#For "moderate303_Patient2":
#    v_gene IGHJ4
#1 IGHV3-33     1
#2 IGHV4-34     1

#For "mild186_Patient3":
#         IGHJ5
#IGHV1-18     1

#For"mild227_Patient4":
#    v_gene IGHJ4
#1 IGHV3-74     1
#2 IGHV4-59     1


#For each sample:

for (i in levels(factor(BCR@meta.data$sample))) {
  BCR@meta.data[BCR$sample == i,] %>%
    group_by(v_gene, j_gene) %>% na.omit() %>% dplyr::count() %>% spread(key = j_gene, value = n) %>%
    as.data.frame() -> matrix
  rownames(matrix) <-  matrix$v_gene
  matrix$v_gene <- NULL
  matrix <- as.matrix(matrix)
  matrix[is.na(matrix)] <- 0
  
  ggsave(filename = paste0("v-j-heatmap-", i, ".jpeg"), path = "./graphs/",
         height = 28, width = 10, dpi = "print",
         plot = pheatmap(matrix, cluster_cols = F, cluster_rows = F, angle_col = 0,
                         labels_row = vapply(strsplit(rownames(matrix), "[V]"), "[", "", 2),
                         labels_col = vapply(strsplit(colnames(matrix), "[J]"), "[", "", 2),
                         cellwidth = 45, cellheight = 45, name = "Frequency", 
                         fontsize = 26, fontsize_col = 60, fontsize_row = 60))
}

#Printing out the most abundant V-J combinations for each sample:
hiclo <- NULL
for (i in levels(factor(BCR@meta.data$sample))) {
  BCR@meta.data[BCR$sample == i,] %>%
    group_by(v_gene, j_gene, sample) %>% na.omit() %>% dplyr::count() %>% arrange(desc(n)) %>%
    as.data.frame() -> matrix
    hiclo <- rbind(hiclo, matrix[1,])
    
}
#The interesting v-j combinations:
#IGHV4-34.IGHJ5
#IGHV1-18.IGHJ3
#IGHV4-39.IGHJ4
#IGHV4-59.IGHJ6
#IGHV4-59.IGHJ4
#IGHV3-48.IGHJ4


Idents(BCR) <- "v_gene"
for (i in levels(factor(BCR[,!is.na(BCR$v_gene)]$azimuthNames))) {
  for (j in unique(all$v_gene)) {
    try({
      degs <- FindMarkers(BCR[, !is.na(BCR$v_gene) & BCR$azimuthNames == i], logfc.threshold = 0.58, ident.1 = j) %>%
      subset(p_val_adj < 0.05) %>% arrange(p_val_adj)
    write.table(x = degs, file = paste0("DEGs-", i, "-", j, ".tsv"), sep = "\t")
    ggsave(filename = paste0("graphs/DEGs-", i, "-", j, ".jpeg"), dpi = "print",
           height = 10, width = 15, units = "in",
           plot = DoHeatmap(BCR[, !is.na(BCR$v_gene)], group.by = "sample", size = 6,
                            features = rownames(degs)))
      
    })
    
    
    
  }
  
}
Idents(BCR) <- "v.j"

degs <- FindMarkers(BCR[, !is.na(BCR$v_gene) & BCR$azimuthNames == "B intermediate" & BCR$v_gene == "IGHV1-18"], 
                    logfc.threshold = 0.58) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)
write.table(x = degs, file = paste0("DEGs-", i, "-", j, ".tsv"), sep = "\t")
ggsave(filename = paste0("graphs/DEGs-", i, "-", j, ".jpeg"), dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(BCR[, !is.na(BCR$v_gene) & BCR$azimuthNames == i], group.by = "sample", size = 6,
                        features = rownames(degs)))
Idents(BCR) <- "v.j"

deg.all <- FindAllMarkers(BCR[, !is.na(BCR$v_gene) & BCR$azimuthNames == "B intermediate"],
               only.pos = T, logfc.threshold = 0.58) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)
deg.all[!grepl("IGHV", deg.all$gene),]
deg.all[!grepl("IGH", deg.all$gene) & !grepl("IGL", deg.all$gene) & !grepl("IGK", deg.all$gene),] %>% view()

degs <- FindMarkers(BCR[, !is.na(BCR$v_gene) & BCR$azimuthNames == "B intermediate" & BCR$v_gene == "IGHV1-18"], 
                    logfc.threshold = 0.58) %>%
  subset(p_val_adj < 0.05) %>% arrange(p_val_adj)
write.table(x = degs, file = paste0("DEGs-", i, "-", j, ".tsv"), sep = "\t")
ggsave(filename = paste0("graphs/DEGs-", i, "-", j, ".jpeg"), dpi = "print",
       height = 10, width = 15, units = "in",
       plot = DoHeatmap(BCR[, !is.na(BCR$v_gene) & BCR$azimuthNames == i], group.by = "sample", size = 6,
                        features = rownames(degs)))

#"healthy1_control1"
#    v_gene j_gene  n
#1 IGHV3-23  IGHJ4 21

#"healthy2_control2"
#    v_gene j_gene  n
#1 IGHV3-23  IGHJ4 14

#"healthy3_control3"
#    v_gene j_gene  n
#1 IGHV3-33  IGHJ4 25

#"healthy4_control4"
#    v_gene j_gene  n
#1 IGHV3-23  IGHJ4 10

#"moderate272_Patient1"
#    v_gene j_gene n
#1 IGHV3-33  IGHJ4 8

#"critical293_Patient1"
#    v_gene j_gene  n
#1 IGHV1-18  IGHJ3 57

#"moderate303_Patient2"
#    v_gene j_gene n
#1 IGHV4-59  IGHJ4 6

#"critical308_Patient2"
#    v_gene j_gene   n
#1 IGHV4-34  IGHJ5 220

#"mild186_Patient3"
#    v_gene j_gene n
#1 IGHV3-48  IGHJ4 4

#"critical213_Patient3"
#    v_gene j_gene  n
#1 IGHV4-39  IGHJ4 42

#"mild227_Patient4"
#    v_gene j_gene n
#1 IGHV3-33  IGHJ4 8

#"critical238_Patient4"
#    v_gene j_gene  n
#1 IGHV3-23  IGHJ4 43

#"critical119_Patient5"
#    v_gene j_gene  n
#1 IGHV3-30  IGHJ4 20

#"severe123_Patient5"
#    v_gene j_gene  n
#1 IGHV3-23  IGHJ4 16

#"moderate138_Patient5"
#    v_gene j_gene  n
#1 IGHV3-23  IGHJ4 15

#"critical120_Patient6"
#    v_gene j_gene  n
#1 IGHV4-59  IGHJ6 35

#"severe122_Patient6"
#    v_gene j_gene  n
#1 IGHV3-23  IGHJ4 10

#"moderate124_Patient6"
#    v_gene j_gene n
#1 IGHV3-23  IGHJ4 6



#Isotypes for each sample:
for (i in levels(factor(BCR@meta.data$sample))) {
  BCR@meta.data[complete.cases(BCR@meta.data) & BCR$sample == i & !is.na(BCR$c_gene),] %>%
    group_by(v_gene, c_gene) %>%  dplyr::count() %>% na.omit() %>% spread(key = c_gene, value = n) %>%
    as.data.frame() -> matrix
  rownames(matrix) <-  matrix$v_gene
  matrix$v_gene <- NULL
  matrix[9] <- NULL
  matrix <- as.matrix(matrix)
  matrix[is.na(matrix)] <- 0
  
  ggsave(filename = paste0("v-c-heatmap-", i, ".jpeg"), path = "./graphs",
         height = 28, width = 12, dpi = "print",
         plot = pheatmap(matrix, cluster_cols = F, cluster_rows = F,
                         cellwidth = 45, cellheight = 45, name = "Frequency", 
                         labels_row = vapply(strsplit(rownames(matrix), "[H]"), "[", "", 2),
                         labels_col = vapply(strsplit(colnames(matrix), "[H]"), "[", "", 2), column_names_rot = 45,
                         fontsize = 26, fontsize_col = 60, fontsize_row = 60))
}
