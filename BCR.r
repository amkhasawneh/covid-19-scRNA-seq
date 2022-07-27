################################Loading#########################################
library(scRepertoire) #remotes::install_github("ncborcherding/scRepertoire@dev")
library(Seurat)
library(scater)
library(SingleCellExperiment)
library(stringr)
library(circlize)
library(scales)
library(dplyr)
library(pals)
library(tidyverse)
library(pheatmap)
library(grid)
library(gridExtra)
library(cowplot)
library(SCpubr)
library(vegan)

#Loading the data:
critical119 <- read.csv("from_cellranger/critical119/vdj_b/filtered_contig_annotations.csv")
critical120 <- read.csv("from_cellranger/critical120/vdj_b/filtered_contig_annotations.csv")
critical213 <- read.csv("from_cellranger/critical213/vdj_b/filtered_contig_annotations.csv")
critical238 <- read.csv("from_cellranger/critical238/vdj_b/filtered_contig_annotations.csv")
critical293 <- read.csv("from_cellranger/critical293/vdj_b/filtered_contig_annotations.csv")
critical308 <- read.csv("from_cellranger/critical308/vdj_b/filtered_contig_annotations.csv")
mild186 <- read.csv("from_cellranger/mild186/vdj_b/filtered_contig_annotations.csv")
mild227 <- read.csv("from_cellranger/mild227/vdj_b/filtered_contig_annotations.csv")
moderate124 <- read.csv("from_cellranger/moderate124/vdj_b/filtered_contig_annotations.csv")
moderate138 <- read.csv("from_cellranger/moderate138/vdj_b/filtered_contig_annotations.csv")
moderate272 <- read.csv("from_cellranger/moderate272/vdj_b/filtered_contig_annotations.csv")
moderate303 <- read.csv("from_cellranger/moderate303/vdj_b/filtered_contig_annotations.csv")
severe122 <- read.csv("from_cellranger/severe122/vdj_b/filtered_contig_annotations.csv")
severe123 <- read.csv("from_cellranger/severe123/vdj_b/filtered_contig_annotations.csv")
hc1 <- read.csv("from_cellranger/hc1/vdj_b/filtered_contig_annotations.csv")
hc2 <- read.csv("from_cellranger/hc2/vdj_b/filtered_contig_annotations.csv")
hc3 <- read.csv("from_cellranger/hc3/vdj_b/filtered_contig_annotations.csv")
hc4 <- read.csv("from_cellranger/hc4/vdj_b/filtered_contig_annotations.csv")


#Combinging the data frames:
contig.list <- list(critical119, critical120, critical213,
                    critical238, critical293, critical308,
                    mild186, mild227, 
                    moderate124, moderate138, moderate272, moderate303,
                    severe122, severe123, 
                    hc1, hc2, hc3, hc4)

B.combined <- combineBCR(contig.list, 
                         samples = c("critical119", "critical120", "critical213",
                                     "critical238", "critical293", "critical308",
                                     "mild186","mild227",
                                     "moderate124", "moderate138", "moderate272", "moderate303",
                                     "severe122", "severe123",
                                     "healthy1", "healthy2", "healthy3", "healthy4"),
                         ID = c("Patient5", "Patient6", "Patient3",
                                "Patient4", "Patient1", "Patient2",
                                "Patient3", "Patient4", 
                                "Patient6", "Patient5", "Patient1", "Patient2",
                                "Patient6", "Patient5",
                                "control1", "control2", "control3", "control4"))

B.combined <- addVariable(B.combined, name = "severity", 
                          variables = c("critical", "critical", "critical",
                                        "critical", "critical", "critical",
                                        "mild", "mild",
                                        "moderate", "moderate", "moderate", "moderate",
                                        "severe", "severe",
                                        "healthy", "healthy", "healthy", "healthy"))
head(B.combined)

#Saving the contigs:
saveRDS(B.combined, "BCR-combined-contigs.rds")

remove(critical119, critical120, critical213,
       critical238, critical293, critical308,
       mild186, mild227, 
       moderate124, moderate138, moderate272, moderate303,
       severe122, severe123,
       hc1, hc2, hc3, hc4, contig.list)
gc()

################################Preliminary BCR analysis########################

B.combined <- readRDS("BCR-combined-contigs.rds")

#Graphs:
abundanceContig(df = B.combined, cloneCall = "gene", group = "sample", scale = F)
abundanceContig(df = B.combined, cloneCall = "gene", group = "sample", scale = F, exportTable = T) %>%
  arrange(desc(Abundance)) -> abundace
write.table(abundace, "abundance.tsv", row.names = F, col.names = T, sep = "\t")

#Calculating IGHV abundance:
quantContig(df = B.combined, cloneCall="gene", chain = "IGH", group="sample", exportTable = T)

lengthContig(df = B.combined, cloneCall="aa", group="sample", scale=T) 

lengthContig(df = B.combined, cloneCall="aa", group="sample",  scale=F) 

lengthContig(df = B.combined, cloneCall = "aa", group = "IGH", scale = F, exportTable = T) %>%
  arrange(desc(length))  -> longest 
longest %>% 
  ggplot(aes(x = IGH, y = length, fill = values)) + geom_bin2d() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

B.combined <-  B.combined[c("healthy1_control1", "healthy2_control2", "healthy3_control3","healthy4_control4",
                            "moderate272_Patient1", "critical293_Patient1", 
                            "moderate303_Patient2", "critical308_Patient2",
                            "mild186_Patient3", "critical213_Patient3",
                            "mild227_Patient4", "critical238_Patient4",
                            "critical119_Patient5", "severe123_Patient5", "moderate138_Patient5",
                            "critical120_Patient6", "severe122_Patient6", "moderate124_Patient6")]

compareClonotypes(df = B.combined, samples = c("healthy1_control1", "healthy2_control2", "healthy3_control3","healthy4_control4",
                                                            "moderate272_Patient1", "critical293_Patient1", 
                                                            "moderate303_Patient2", "critical308_Patient2",
                                                            "mild186_Patient3", "critical213_Patient3",
                                                            "mild227_Patient4", "critical238_Patient4",
                                                            "critical119_Patient5", "severe123_Patient5", "moderate138_Patient5",
                                                            "critical120_Patient6", "severe122_Patient6", "moderate124_Patient6"),
                  cloneCall="gene", graph = "alluvial") 

compareClonotypes(df = B.combined, samples = c("healthy1_control1", "healthy2_control2", "healthy3_control3","healthy4_control4",
                                               "moderate272_Patient1", "critical293_Patient1", 
                                               "moderate303_Patient2", "critical308_Patient2",
                                               "mild186_Patient3", "critical213_Patient3",
                                               "mild227_Patient4", "critical238_Patient4",
                                               "critical119_Patient5", "severe123_Patient5", "moderate138_Patient5",
                                               "critical120_Patient6", "severe122_Patient6", "moderate124_Patient6"),
                  cloneCall="gene", exportTable = T) %>%
  arrange(desc(Proportion))   -> top.clonotypes
 
write.table(top.clonotypes, "top-clonotypes.tsv", sep = "\t", row.names = F) 

clonalHomeostasis(df = B.combined, cloneCall = "gene",
                  cloneTypes = c(Rare = 1e-04, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded= 1))

clonalProportion(df = B.combined, cloneCall = "gene", split = c(10, 100, 1000, 10000, 30000, 1e+05)) 

clonalOverlap(df = B.combined, cloneCall = "gene", method = "overlap") +
  theme(axis.text = element_text(size = 12), axis.text.x = element_text(angle = 90))


clonalDiversity(df = B.combined[c(11, 5)], cloneCall = "gene", group = "sample")
clonalDiversity(df = B.combined, cloneCall = "gene", group = "sample", exportTable = T) %>%
  write.table("clonal-diversity.tsv", sep = "\t", row.names = F)

saveRDS(B.combined, "BCR-combined-contigs.rds")
################################Integration with Seurat object##################

B.combined <- readRDS("BCR-combined-contigs.rds")

#Loading clustered Seurat object:
covid <- readRDS("04-covid-clustered.rds")

#Changing the current Idents:
Idents(covid) <- "azimuthNames"

#Subsetting B-cells:
BCR <- covid[,grep("(B)|(Plasmablast)", covid$azimuthNames)]

#Combining Seurat object with BCR table:
BCR <- combineExpression(B.combined, BCR, cloneCall="gene", filterNA = T)
head(BCR@meta.data)
gc()

#Adjusting the sample variable:
BCR$sample <- factor(BCR$sample, levels = c("healthy1_control1", "healthy2_control2", "healthy3_control3","healthy4_control4",
                                            "moderate272_Patient1", "critical293_Patient1", 
                                            "moderate303_Patient2", "critical308_Patient2",
                                            "mild186_Patient3", "critical213_Patient3",
                                            "mild227_Patient4", "critical238_Patient4",
                                            "critical119_Patient5", "severe123_Patient5", "moderate138_Patient5",
                                            "critical120_Patient6", "severe122_Patient6", "moderate124_Patient6"))

#Adding V and J gene usage:
BCR$v_gene <- vapply(strsplit(BCR$CTgene, "[.]"), "[", "", 1)
BCR$v_gene <- sub("(NA)", NA, BCR$v_gene)
BCR$j_gene <- vapply(strsplit(BCR$CTgene, "[.]"), "[", "", 2)
BCR$j_gene <- gsub("(.*[KL].*)", NA, BCR$j_gene)
BCR$c_gene <- vapply(strsplit(vapply(strsplit(BCR$CTgene, "[_]"), "[", "", 1), "[.]"), "[", "", 4)

saveRDS(BCR, "05-BCR-combined.rds")

BCR <- readRDS("05-BCR-combined.rds")

remove(covid, B.combined, abundace, top.clonotypes, longest)
gc()

#Saving some tables:
write.table(BCR@meta.data[order(BCR@meta.data[["Frequency"]],decreasing=TRUE),],
            file = "BCR-frequency.tsv", sep="\t", append = FALSE, quote=FALSE, 
            row.names = FALSE, col.names = TRUE)

#Tables with cell numbers for each sample:
cells <- BCR@meta.data %>%
  group_by(sample, azimuthNames) %>% count() %>% 
  spread(key = sample, value = n) %>% as.data.frame()
rownames(cells) <- cells$azimuthNames
cells$azimuthNames <- NULL
cells <- as.matrix(cells)
cells <- proportions(as.matrix(cells), margin = 2) * 100
write.table(cells, file = "cell-proportions-samples.tsv", sep = "\t", col.names = NA)


################################Simple gene usage###############################
BCR <- readRDS("05-BCR-combined.rds")

heavy.chains <- BCR@meta.data %>%
  group_by(v_gene, j_gene, c_gene, azimuthNames, patient, severity) %>% dplyr::count() %>% arrange(desc(n)) %>% 
  as.data.frame() %>% na.omit()


ggplot(data = heavy.chains[heavy.chains$n > 1,], aes(x= v_gene, y= j_gene, color = azimuthNames, size = n)) +
  geom_count(na.rm = T) +
  labs(x="V Gene", y="J Gene", color="Cells", size="Count") + ggtitle("Clonotypes & Cells") +
  theme(axis.text.x = element_text(angle = 90, size = 20), axis.title = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.y = element_text(size = 20), legend.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 20))
ggplot(data = heavy.chains[heavy.chains$n > 1,], aes(x=v_gene, fill = azimuthNames)) +
  geom_histogram(stat = "count", position = "fill") + 
  facet_wrap(~patient) +
  scale_y_continuous(labels = percent_format()) +
  labs(x="V Gene", y="Proportion", color="Cells") + ggtitle("V Genes & Cells") +
  theme(axis.text.x = element_text(angle = 90, size = 12), axis.title = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.y = element_text(size = 20), legend.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 20))

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
    plot <- pheatmap(matrix, cluster_cols = F, cluster_rows = F, height = 20, width = 8,
             cellwidth = 32, cellheight = 32, name = paste0(BCR$patient, "-", BCR$severity),
             fontsize_col = 42, fontsize_row = 42, fontsize = 25)
    ggsave(filename = paste0("v-j-heatmap-", i, "-", j, "2.jpeg"), plot = plot, path = "./graphs/",
           dpi = "print", height = 20, width = 8)
    
    })
  }
}


for (i in levels(factor(BCR$patient))) {
  if (i == "Patient3" | i == "Patient4"){
  
    ggsave(filename = paste0(i, "-b-cells.jpeg"), path = "./graphs/", dpi = "print",
           plot = ggplot(data = BCR@meta.data[BCR$patient == i,], aes(x = severity, fill = azimuthNames)) +
  geom_histogram(stat = "count", position = "fill") +
    xlim(c("critical", "severe", "moderate")) +
    scale_y_continuous(name = "Percentage",
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) +
      xlab("Severity"))  
  
  
  } else {
    ggsave(filename = paste0(i, "-b-cells.jpeg"), path = "./graphs/", dpi = "print",
           plot = ggplot(data = BCR@meta.data[BCR$patient == i,], aes(x = severity, fill = azimuthNames)) +
      geom_histogram(stat = "count", position = "fill") +
      scale_y_continuous(name = "Percentage",
                         breaks = c(0, 0.25, 0.5, 0.75, 1),
                         labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) +
      xlab("Severity"))
        }
}

for (i in levels(factor(BCR$sample))) {
  
    ggsave(filename = paste0(i, "-b-cells-distribution.jpeg"), path = "./graphs/", dpi = "print",
           plot = ggplot(data = BCR@meta.data[BCR$sample == i,], aes(x = severity, fill = azimuthNames)) +
  geom_histogram(stat = "count", position = "fill") +
    scale_y_continuous(name = "",
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) +
    xlab("") + guides(x = "none"))
  
  
  }

ggplot(data = heavy.chains[heavy.chains$patient == i,], aes(x = severity, fill = azimuthNames)) +
  geom_bar() +
  geom_text(aes(label = n, y = n), color="green", position = position_stack(vjust = 0.5)) 

#Re-ordering abundance from highest to lowest:
abundance <- abundanceContig(BCR, cloneCall = "gene", split.by = "sample", 
                             exportTable = T, scale = F) %>% arrange(desc(Abundance)) %>%
  head(40)

abundanceContig(BCR, cloneCall = "gene", split.by = "sample",
                exportTable = F, scale = F)


write.table(BCR@meta.data, "seurat-clonotypes.tsv", sep = "\t", quote = F, row.names = F)

saveRDS(BCR, "05-BCR-combined")

BCR@meta.data[complete.cases(BCR@meta.data),] %>%
  group_by(v_gene, j_gene) %>% dplyr::count(sort = T)

#V-J gene usage heatmaps:
#For all:
BCR@meta.data[complete.cases(BCR@meta.data),] %>%
  group_by(v_gene, j_gene) %>% dplyr::count() %>% spread(key = v_gene, value = n) %>% as.data.frame() -> matrix
rownames(matrix) <-  matrix$j_gene
matrix$j_gene <- NULL
matrix <- as.matrix(matrix)
matrix[is.na(matrix)] <- 0
ggsave(filename = paste0("v-j-heatmap.jpeg"), path = "./graphs",
       height = 6 , width = 26,
       plot = pheatmap(matrix, cluster_cols = F, cluster_rows = F,main = "All",
                       cellwidth = 30, cellheight = 30, name = "Frequency", fontsize = 20))

#For each severity group:
for (i in levels(factor(BCR$severity))) {
  BCR@meta.data[complete.cases(BCR@meta.data) & BCR$severity == i,] %>%
    group_by(v_gene, j_gene) %>% dplyr::count() %>% spread(key = j_gene, value = n) %>% as.data.frame() -> matrix
  rownames(matrix) <-  matrix$v_gene
  matrix$v_gene <- NULL
  matrix <- as.matrix(matrix)
  matrix[is.na(matrix)] <- 0
  
  ggsave(filename = paste0("v-j-heatmap-", i, ".jpeg"), path = "./graphs",
         height = 24, width = 10,
         plot = pheatmap(matrix, cluster_cols = F, cluster_rows = F,
                         cellwidth = 30, cellheight = 30, name = "Frequency", fontsize = 20))
}

#For each sample:
for (i in levels(factor(BCR@meta.data[complete.cases(BCR$sample),]$sample))) {
  BCR@meta.data[complete.cases(BCR@meta.data) & BCR$sample == i,] %>%
    group_by(v_gene, j_gene) %>% dplyr::count() %>% spread(key = j_gene, value = n) %>% as.data.frame() -> matrix
  rownames(matrix) <-  matrix$v_gene
  matrix$v_gene <- NULL
  matrix <- as.matrix(matrix)
  matrix[is.na(matrix)] <- 0
  
  ggsave(filename = paste0("v-j-heatmap-", i, ".jpeg"), path = "./graphs",
         height = 22, width = 10,
         plot = pheatmap(matrix, cluster_cols = F, cluster_rows = F,
                         main = paste0(BCR$patient[BCR$sample == i], " ", BCR$severity[BCR$sample == i]), 
                         cellwidth = 30, cellheight = 30, name = "Frequency", fontsize = 20))
}

#V gene usage by isotype:
#For all:
BCR@meta.data[complete.cases(BCR@meta.data),] %>%
  group_by(v_gene, c_gene) %>% dplyr::count() %>% spread(key = c_gene, value = n) %>% as.data.frame() -> matrix
rownames(matrix) <-  matrix$v_gene
matrix$v_gene <- NULL
matrix <- as.matrix(matrix)
matrix[is.na(matrix)] <- 0
ggsave(filename = paste0("v-c-heatmap.jpeg"), path = "./graphs",
       height = 24 , width = 10,
       plot = pheatmap(matrix, cluster_cols = F, cluster_rows = F, main = "All",
                       cellwidth = 30, cellheight = 30, name = "Frequency", fontsize = 20))

#For each sample:
for (i in levels(factor(BCR@meta.data$sample))) {
  BCR@meta.data[complete.cases(BCR@meta.data) & BCR$sample == i,] %>%
    group_by(v_gene, c_gene) %>% dplyr::count() %>% spread(key = c_gene, value = n) %>% as.data.frame() -> matrix
  rownames(matrix) <-  matrix$v_gene
  matrix$v_gene <- NULL
  matrix <- as.matrix(matrix)
  matrix[is.na(matrix)] <- 0
  
  ggsave(filename = paste0("v-c-heatmap-", i, ".jpeg"), path = "./graphs",
                  height = 22, width = 10,
         plot = pheatmap(matrix, cluster_cols = F, cluster_rows = F,
                         main = paste0(BCR$patient[BCR$sample == i], " ", BCR$severity[BCR$sample == i]),
                         cellwidth = 30, cellheight = 30, name = "Frequency", fontsize = 20))
}

#Isotype by patient:
BCR@meta.data[complete.cases(BCR@meta.data),] %>%
  ggplot(aes(x=patient, fill=as.factor(c_gene))) + geom_bar(stat = "count", position = "dodge") + facet_wrap(~severity) +
  labs(x="V Gene", y="Count", fill="Isotype") + ggtitle("Isotype by Patient") +
  theme(axis.text.x = element_text(angle = 90, size = 20), axis.title = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.y = element_text(size = 20), legend.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, color = "blue", size = 20))

#CDR3 length by isotype:
ggplot(data = BCR@meta.data[complete.cases(BCR@meta.data) ,], aes(x=nchar(CTaa), fill = c_gene)) +
  geom_bar(stat = "count", position = "dodge") +
  labs(x="CDR3 Length", y="Frequency", fill="Isotype") + ggtitle("CDR3 Length by Isotype") +
  theme(axis.text.x = element_text(size = 20), axis.title = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.y = element_text(size = 20), legend.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, color = "blue", size = 20))

#Isotype by cell type:
ggplot(data = BCR@meta.data[complete.cases(BCR@meta.data),], aes(x=c_gene, fill=azimuthNames)) +
  geom_bar(stat = "count", position = "dodge") + 
  labs(x="Isotype", y="Frequency", fill="Cells") + ggtitle("Isotypes by Cell type") +
  theme(axis.text.x = element_text(size = 20, angle = 90), axis.title = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text.y = element_text(size = 20), legend.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, color = "blue", size = 20))

#Isotype per severity group:
BCR@meta.data[complete.cases(BCR@meta.data),] %>%
  group_by(severity, c_gene) %>% dplyr::count() %>% spread(key = c_gene, value = n) %>% as.data.frame() -> matrix
rownames(matrix) <-  matrix$severity
matrix$severity <- NULL
matrix <- as.matrix(matrix)
matrix[is.na(matrix)] <- 0
ggsave(filename = "isotype-severity.jpeg", path = "./graphs",
       height = 4, width = 8,
       plot = pheatmap(matrix, cellwidth = 30, cellheight = 30, name = "Frequency", fontsize = 20))

#Isotype per patient:
BCR@meta.data[complete.cases(BCR@meta.data),] %>%
  group_by(sample, c_gene) %>% dplyr::count() %>% spread(key = c_gene, value = n) %>% as.data.frame() -> matrix
rownames(matrix) <-  matrix$sample
matrix$sample <- NULL
matrix <- as.matrix(matrix)
matrix[is.na(matrix)] <- 0
ggsave(filename = "isotype-sample.jpeg", path = "./graphs",
       height = 8, width = 8,
       plot = pheatmap(matrix, labels_row = c("HC1", "HC2", "HC3","HC4",
                                              "Patient1 moderate", "Patient1 critical",
                                              "Patient2 mild", "Patient2 critical",
                                              "Patient3 moderate", "Patient3 critical",
                                              "Patient4 mild", "Patient4 critical", 
                                              "Patient5 critical", "Patient5 severe", "Patient5 moderate",
                                              "Patient6 critical", "Patient6 severe", "Patient6 moderate"),
                       main = "Isotype per Sample", cluster_rows = F, cluster_cols = F,
                       cellwidth = 30, cellheight = 30, name = "Frequency", fontsize = 20))

################################Complicated stuff###############################

BCR <- readRDS("05-BCR-combined.rds")


#Creating UMAPs of B-cell clonotypes based on Azimuth cell group predictions:
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
DimPlot(BCR, group.by = "cloneType",reduction = "umap", pt.size=1.1) +
  scale_color_manual(values = c(colorblind_vector(5)), na.value="grey") 
DimPlot(BCR, group.by = "orig.ident",reduction = "umap", pt.size=1.1) +
  scale_color_manual(values = c(colorblind_vector(10)), na.value="grey")
DimPlot(BCR, group.by = "azimuthNames",reduction = "umap", pt.size=1.1)
DimPlot(BCR, group.by = "severity",reduction = "umap", pt.size=1.1)

#Re-ordering abundance from highest to lowest:
abundance <- abundanceContig(BCR, cloneCall = "gene", split.by = "sample",
                exportTable = T, scale = F) %>% arrange(desc(Abundance)) 

abundanceContig(BCR, cloneCall = "gene", split.by = "sample",
                exportTable = F, scale = F)

#Adding V and J gene usage:
abundance$v_gene <- vapply(strsplit(abundance$CTgene, "[.]"), "[", "", 1)
abundance$v_gene <- sub("(NA)", NA, abundance$v_gene)
abundance$j_gene <- vapply(strsplit(abundance$CTgene, "[.]"), "[", "", 2)
abundance$j_gene <- gsub("(.*[KL].*)", NA, abundance$j_gene)
abundance$c_gene <- vapply(strsplit(vapply(strsplit(abundance$CTgene, "[_]"), "[", "", 1), "[.]"), "[", "", 4)
ggplot(data = head(abundance, 10), aes(color = values, x = CTgene, y = Abundance)) + geom_jitter()

head(abundance, 57)

abundance[abundance$CTgene == "IGHV1-18.IGHJ3..IGHA1_IGLV1-47.IGLJ2.IGLC2",]
abundance[abundance$values == "moderate272_Patient1",]
top20 <- abundance[order(abundance$Abundance,decreasing = TRUE),][1:20,] 
write.table(abundance, "abundant-clonotypeses.tsv", sep = "\t", quote = F, row.names = F)

#Most abundant CT gene:
#IGHV1-18.IGHJ3..IGHA1_IGLV1-47.IGLJ2.IGLC2

Idents(BCR) <- "azimuthNames"
DimPlot(BCR, pt.size = 2, order = TRUE)

vizGenes(BCR, gene = "V", chain = "IGH", scale = T,
         y.axis = "J", order = "gene") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 22),
          axis.text.y = element_text(vjust = 0.5, hjust=1, size = 10))

#To see the top 10 genes:
hold <- B.combined[[as.character(abundance[order(abundance$Abundance,decreasing = TRUE),][1,3])]]
hold[hold$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][1,1], "CTgene"][1]

hold2 <- B.combined[[as.character(abundance[order(abundance$Abundance,decreasing = TRUE),][2,3])]]
hold2[hold2$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][2,1], "CTgene"][1]

hold3 <- B.combined[[as.character(abundance[order(abundance$Abundance,decreasing = TRUE),][3,3])]]
hold3[hold3$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][3,1], "CTgene"][1]

hold4 <- B.combined[[as.character(abundance[order(abundance$Abundance,decreasing = TRUE),][4,3])]]
hold4[hold4$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][4,1], "CTgene"][2]

hold5 <- B.combined[[as.character(abundance[order(abundance$Abundance,decreasing = TRUE),][5,3])]]
hold5[hold5$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][5,1], "CTgene"][1]

hold6 <- B.combined[[as.character(abundance[order(abundance$Abundance,decreasing = TRUE),][6,3])]]
hold6[hold6$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][6,1], "CTgene"][1]

hold7 <- B.combined[[as.character(abundance[order(abundance$Abundance,decreasing = TRUE),][7,3])]]
hold7[hold7$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][7,1], "CTgene"][1]

hold8 <- B.combined[[as.character(abundance[order(abundance$Abundance,decreasing = TRUE),][8,3])]]
hold8[hold8$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][8,1], "CTgene"][1]

hold9 <- B.combined[[as.character(abundance[order(abundance$Abundance,decreasing = TRUE),][9,3])]]
hold9[hold9$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][9,1], "CTgene"][1]

hold10 <- B.combined[[as.character(abundance[order(abundance$Abundance,decreasing = TRUE),][10,3])]]
hold10[hold9$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][10,1], "CTgene"][1]

#Top 10 genes: 
#"IGHV3-23.IGHJ4..IGHA2_IGKV1-5.IGKJ1.IGKC"
#"IGHV1-46.IGHJ6..IGHA1_IGLV1-44.IGLJ2.IGLC3"
#"IGHV1-18.IGHJ4..IGHG4_IGKV3-20.IGKJ2.IGKC"
#"IGHV3-72.IGHJ3..IGHA1_IGKV3-20.IGKJ1.IGKC"
#"IGHV4-34.IGHJ5..IGHA1_IGKV1D-39.IGKJ2.IGKC"
#"IGHV4-30-4.IGHJ6.IGHD6-6.IGHG1_IGKV1-16.IGKJ4.IGKC"
#"IGHV3-33.IGHJ5..IGHG1_IGLV2-14.IGLJ3.IGLC2"
#"IGHV2-5.IGHJ4..IGHG1_IGKV1-12.IGKJ1.IGKC"
#"IGHV3-33.IGHJ4.IGHD6-19.IGHG1_IGLV3-10.IGLJ3.IGLC2"
#"IGHV4-39..IGHJ4.IGHG1_IGKV1-5.IGKJ1.IGKC"

#By clonotype gene:
DimPlot(BCR, pt.size =3,order = TRUE, shape.by = "azimuthNames") +
  scale_color_manual(labels = c(
    hold[hold$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][1,1],"CTgene"][1],
    hold2[hold2$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][2,1],"CTgene"][1],
    hold3[hold3$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][3,1],"CTgene"][1],
    hold4[hold4$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][4,1],"CTgene"][2],
    hold5[hold5$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][5,1],"CTgene"][1],
    hold6[hold6$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][6,1],"CTgene"][1],
    hold7[hold7$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][7,1],"CTgene"][1],
    hold8[hold8$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][8,1],"CTgene"][1],
    hold9[hold9$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][9,1],"CTgene"][1],
    hold9[hold9$CTgene==abundance[order(abundance$Abundance,decreasing = TRUE),][10,1],"CTgene"][1]),
    values = c("red","pink","orange","blue","lightblue","green","lightgreen","yellow","purple")) +
  labs(color = "CTgene") + theme(legend.text = element_text(size = 22))

#Clonal proportions:
occ.m <- occupiedscRepertoire(BCR[,BCR$severity == "mild"], x.axis = "azimuthNames") +
  theme(legend.text = element_text(size = 22)) + ggtitle("Mild") +
  geom_text(aes(label = value), color="green", position = position_stack(vjust = 0.5)) 
occ.s <- occupiedscRepertoire(BCR[,BCR$severity == "severe"], x.axis = "azimuthNames") +
  theme(legend.text = element_text(size = 22)) + ggtitle("Severe") +
  geom_text(aes(label = value), color="green", position = position_stack(vjust = 0.5)) 
occ.mo <- occupiedscRepertoire(BCR[,BCR$severity == "moderate"], x.axis = "azimuthNames") +
  theme(legend.text = element_text(size = 22)) + ggtitle("Moderate") +
  geom_text(aes(label = value), color="green", position = position_stack(vjust = 0.5)) 
occ.c <- occupiedscRepertoire(BCR[,BCR$severity == "critical"], x.axis = "azimuthNames") +
  theme(legend.text = element_text(size = 22)) + ggtitle("Critical") +
  geom_text(aes(label = value), color="green", position = position_stack(vjust = 0.5)) 
ggsave(filename = "clonal-proportions.jpeg", path = "graphs/",
       width = 1920, height = 1080, units = "px", dpi = 120,
       plot = grid.arrange(occ.m, occ.mo, occ.s, occ.c, top= textGrob("Clonal Proportion", gp=gpar(fontsize=30,font=3)))
          )

occupiedscRepertoire(BCR, x.axis = "sample")

#Tracking clonotypes:
ggsave("patient-severity-cell.jpeg", width = 1920, height = 1080, dpi = 300, units = "px",
       plot = alluvialClonotypes(BCR, cloneCall = "gene", 
                   y.axes = c("patient", "severity", "azimuthNames"), color = "azimuthNames", facet = NULL))

#Chord diagram:
circles <- getCirclize(BCR, group.by = "cluster")
grid.cols <- hue_pal()(length(unique(BCR$orig.ident)))
names(grid.cols) <- as.character(unique(BCR$orig.ident))
chordDiagram(circles, self.link = 1, grid.col = grid.cols)

#Chord diagram:
circles <- getCirclize(BCR, group.by = c("v_gene", "j_gene"))
grid.cols <- hue_pal()(length(unique(BCR$orig.ident)))
names(grid.cols) <- as.character(unique(BCR$orig.ident))
chordDiagram(circles, self.link = 1, grid.col = grid.cols)

for (i in levels(BCR$patient)) {
  circles <- getCirclize(BCR[,BCR$severity == i], group.by = "cluster")
  grid.cols <- hue_pal()(length(unique(BCR$sample)))
  names(grid.cols) <- as.character(unique(BCR$sample))
  jpeg(filename = paste0("graphs/shared-clonotypes-", i, ".jpeg"),
       width = 1920, height = 1080,
       units = "px", res = 300)
  chordDiagram(circles, self.link = 1, grid.col = grid.cols)
  dev.off()
}

#Chord diagram for shared clonotypes:
circles <- getCirclize(BCR, group.by = "sample")
grid.cols <- hue_pal()(length(unique(BCR$sample)))
names(grid.cols) <- as.character(unique(BCR$sample))
jpeg(filename = "graphs/shared-clonotypes.jpeg", width = 4000, height = 4000,
     units = "px", res = 300)
chordDiagram(circles, self.link = 1, grid.col = grid.cols)
dev.off()

#The number of cells in each severity category:
BCR@meta.data %>%
  group_by(sample, orig.ident) %>%
  count(sort = T) %>%
  write.table(file = "b-cells-patients.tsv", sep = "\t", row.names = F)

#The number of cells in each category:
BCR@meta.data %>%
  group_by(azimuthNames) %>%
  count(sort = T) %>%
  write.table(file = "b-cell-counts.tsv", sep = "\t", row.names = F)

ggplot(data = BCR@meta.data, aes(x = sample, fill = azimuthNames)) + geom_bar()

clonalNetwork(BCR, 
              reduction = "umap", 
              identity = "cluster",
              filter.clones = NULL,
              filter.identity = NULL,
              cloneCall = "aa", exportTable = T) %>% 
  spread(key = "to", value = "weight") 

BCR@meta.data[complete.cases(BCR@meta.data$CTgene),] %>%
  group_by(CTgene) %>%
  count(sort = T) %>% head(20)

contig[contig$chain == "IGH",] %>%
  group_by(v_gene, j_gene) %>%
  count(sort = T) %>% head(20)

Idents(BCR) <- "CTgene"
FindMarkers(object = BCR, ident.1 = "IGHV1-18.IGHJ3..IGHA1_IGLV1-47.IGLJ2.IGLC2",
            group.by = "CTgene") %>%
  arrange(p_val_adj) %>% filter(p_val_adj < 0.05) -> clonotypeMarkerGenes



################################Diversity Testing###############################


#V gene diversity, using the Inverse Simpson index:
div <- data.frame()
for (i in levels(BCR$sample)) {
  v.div <- c(i, diversity(BCR[,BCR$sample == i]$v_gene %>% table, index = "invsimpson"), levels(factor(BCR$outcome[BCR$sample == i])))
  div <- rbind(div, v.div)
  rownames(div) <- div[,1]
  colnames(div) <- c("sample", "Inv.Simpson.score", "outcome")
}

#Comparing diversity between outcome groups:
#Recovered vs. Healthy (p-value = 0.02424):
wilcox.test(x = as.numeric(div$Inv.Simpson.score[div$outcome == "Recovered"]), 
            y = as.numeric(div$Inv.Simpson.score[div$outcome == "Healthy"]))
#Recovered vs. Deceased (p-value = 0.662):
wilcox.test(x = as.numeric(div$Inv.Simpson.score[div$outcome == "Recovered"]), 
            y = as.numeric(div$Inv.Simpson.score[div$outcome == "Deceased"]))
#Deceased vs. Healthy (p-value = 0.02381):
wilcox.test(x = as.numeric(div$Inv.Simpson.score[div$outcome == "Deceased"]), 
            y = as.numeric(div$Inv.Simpson.score[div$outcome == "Healthy"]))
#Healhty vs. not Healthy(p-value = 0.005882):
wilcox.test(x = as.numeric(div$Inv.Simpson.score[div$outcome == "Healthy"]), 
            y = as.numeric(div$Inv.Simpson.score[div$outcome != "Healthy"]))
#I think we can conclude that there is no difference between recovered and dead patients,
#but that healthy controls have significantly higher diversity.
