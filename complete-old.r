#################################Loading########################################
#Loading libraries:
library(Seurat)
library(SeuratDisk)
library(Azimuth)
library(tidyverse)
library(patchwork)
library(Matrix)
library(pals)
library(future)

#Loading data:
mild227 <- Read10X("from_cellranger\\mild227\\sample_feature_bc_matrix")
moderate124 <- Read10X("from_cellranger\\moderate124\\sample_feature_bc_matrix")
moderate138 <- Read10X("from_cellranger\\moderate138\\sample_feature_bc_matrix")
moderate272 <- Read10X("from_cellranger\\moderate272\\sample_feature_bc_matrix")
severe122 <- Read10X("from_cellranger\\severe122\\sample_feature_bc_matrix")
severe123 <- Read10X("from_cellranger\\severe123\\sample_feature_bc_matrix")
critical119 <- Read10X("from_cellranger\\critical119\\sample_feature_bc_matrix")
critical120 <- Read10X("from_cellranger\\critical120\\sample_feature_bc_matrix")
critical238 <- Read10X("from_cellranger\\critical238\\sample_feature_bc_matrix")
critical293 <- Read10X("from_cellranger\\critical293\\sample_feature_bc_matrix")

#Creating Seurat objects:
mild227 <- CreateSeuratObject(counts = mild227, min.cells = 3, min.features = 100)
moderate124 <- CreateSeuratObject(counts = moderate124, min.cells = 3, min.features = 100)
moderate138 <- CreateSeuratObject(counts = moderate138, min.cells = 3, min.features = 100)
moderate272 <- CreateSeuratObject(counts = moderate272, min.cells = 3, min.features = 100)
severe122 <- CreateSeuratObject(counts = severe122, min.cells = 3, min.features = 100)
severe123 <- CreateSeuratObject(counts = severe123, min.cells = 3, min.features = 100)
critical119 <- CreateSeuratObject(counts = critical119, min.cells = 3, min.features = 100)
critical120 <- CreateSeuratObject(counts = critical120, min.cells = 3, min.features = 100)
critical238 <- CreateSeuratObject(counts = critical238, min.cells = 3, min.features = 100)
critical293 <- CreateSeuratObject(counts = critical293, min.cells = 3, min.features = 100)

covid <- merge(x = critical119, y = c(critical120, critical238, critical293, mild227,
                                      moderate124, moderate138, moderate272,
                                      severe122, severe123), 
               add.cell.ids = c("critical119_Patient3", "critical120_Patient4", "critical238_Patient2", "critical293_Patient1", 
                                "mild227_Patient2", "moderate124_Patient4", "moderate138_Patient3", "moderate272_Patient1",
                                "severe122_Patient4", "severe123_Patient3"),
               project = "covid-19")

head(covid@meta.data)

#Saving file for upload to Azimuth:
saveRDS(covid, "00-covid-raw.rds")

remove(critical119, critical120, critical238, critical293, mild227,
                          moderate124, moderate138, moderate272,
                          severe122, severe123)
gc()

#################################Environment####################################
#Preparing the environment:
#Adjusting the limit for allowable R object sizes: 
options(future.globals.maxSize = 8000 * 1024^2)

#Expanding memory:
memory.limit(size=5600000)

#Enabling parallelization:
plan("multicore", workers=8)

gc()

#################################QC#############################################

covid <- readRDS("00-covid-raw.rds")

#Adding the number of genese per UMI:
covid$log10GenesPerUMI <- log10(covid$nFeature_RNA)/log10(covid$nCount_RNA)


#Cell-level filtering#
#Adding mitochondrial gene percentages:
covid$mitoRatio <- PercentageFeatureSet(covid, pattern = "^MT-")/100
VlnPlot(covid, features = c("nCount_RNA", "nFeature_RNA", "mitoRatio"), ncol = 3)
plot1 <- FeatureScatter(covid, feature1 = "nCount_RNA", feature2 = "mitoRatio")
plot2 <- FeatureScatter(covid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Separating the metadat for a while:
metadata <- covid@meta.data

#Adding cell IDs to the metadata:
metadata$cell <- rownames(metadata)
head(metadata)

#Renaming columns:
metadata <- metadata %>%
  rename(seq.folder = orig.ident,
         nUMI = nCount_RNA,
         nGene = nFeature_RNA)
head(metadata)

#Creating a severity column:
metadata$severity <- NA
metadata$severity[which(str_detect(metadata$cell, ".*mild"))] <- "mild"
metadata$severity[which(str_detect(metadata$cell, ".*moderate"))] <- "moderate"
metadata$severity[which(str_detect(metadata$cell, ".*severe"))] <- "severe"
metadata$severity[which(str_detect(metadata$cell, ".*critical"))] <- "critical"

#Creating a sample column:
metadata$sample <- NA
metadata$sample <- sub("(.*?)_{1}(.*?)($|-.*)", "\\1", rownames(metadata))

#Creating a patient column:
metadata$patient <- sub("(.*?)_", "\\2", metadata$sample)

#Creating an outcome column:
metadata$outcome <- "Recovered"
metadata$outcome[metadata$patient == "Patient1",] <- "Deceased"

#Putting metadata back into Seurat object:
covid@meta.data <- metadata

head(covid@meta.data)

#Visualizing the number of cell counts per severity group:
covid@meta.data %>% 
  ggplot(aes(x = severity, fill = severity)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Cell counts")

#Visualizing UMIs/transcript per cell:
covid@meta.data %>% 
  ggplot(aes(color = severity, x = nUMI, fill = severity)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_vline(xintercept = 500) +
  ggtitle("UMIs/transcript per cell")

#Visualizing the distribution of genes detected per cell via histogram:
covid@meta.data %>% 
  ggplot(aes(color = severity, x = nGene, fill = severity)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Genes per cell")

#Visualizing the distribution of genes detected per cell via boxplot:
covid@meta.data %>% 
  ggplot(aes(x = severity, y = log10(nGene), fill = severity)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

#Visualizing the correlation between genes and UMIs and determining whether 
#there's a strong presence of cells with low numbers of genes/UMIs:
covid@meta.data %>% 
  ggplot(aes(x = nUMI, y = nGene, color = mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method = "lm") +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~severity) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Correlation between genes detected and number of UMIs")

#Visualizing the distribution of mitochondrial gene expression per cell:
covid@meta.data %>% 
  ggplot(aes(color = severity, x = mitoRatio, fill = severity)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Mitochondrial gene expression detected per cell")

#Visualizing the overall gene expression complexity by visualizing the genes detected per UMI:
covid@meta.data %>%
  ggplot(aes(x = log10GenesPerUMI, color = severity, fill = severity)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Overall gene expression complexity")

#Sub-setting based on gene and UMI counts and mitochondrial genes:
covid <- subset(covid, subset = (nUMI >= 500) & (nGene > 250) &  
                  (log10GenesPerUMI > 0.80) & (mitoRatio < 0.2))

#Removing genes with zero counts:
counts <- GetAssayData(covid, slot = "counts")
nonzero <- counts > 0
head(nonzero)

#Summing up non-zeros and returning genes with 10 or more values:
keep <- rowSums(counts) >= 10

#Creating final(?) Seurat object:
covid <- CreateSeuratObject(counts = counts[keep,],
                               meta.data = covid@meta.data)
head(covid@meta.data)
covid$nCount_RNA <- NULL
covid$nFeature_RNA <- NULL

remove(counts, metadata, nonzero, keep, plot1, plot2)
gc()

#Saving filtered file:
saveRDS(covid, "01-covid-filtered.rds")

covid <- readRDS("01-covid-filtered.rds")

#Cell cycle scoring#
#Loading cell cycle markers:
load("cycle.rda")

#Normalizing the counts:
covid <- NormalizeData(covid)

#Cell cycle scoring:
covid <- CellCycleScoring(covid, s.features = s_genes, g2m.features = g2m_genes)
head(covid@meta.data)

#Identifying the most variable genes:
covid <- FindVariableFeatures(covid, selection.method = "vst",
                              nfeatures = 5000)
VariableFeatures(covid) %>% head(20) -> top20
LabelPoints(plot = VariableFeaturePlot(covid), points = top20, repel = T)

remove(top20, s_genes, g2m_genes)
gc()

saveRDS(covid, "01-covid-filtered.rds")


#################################Split-Scaling & Azimuth########################

covid <- readRDS("01-covid-filtered.rds")

#Splitting Seurat object by sample to perform scaling on all samples:
#First splitting into samples:
split.covid <- SplitObject(covid, split.by = "sample")

#Scaling:
for (i in 1:length(split.covid)) {
  split.covid[[i]] <- NormalizeData(split.covid[[i]], verbose = T)
  split.covid[[i]] <- FindVariableFeatures(split.covid[[i]], verbose = T)
  split.covid[[i]] <- ScaleData(split.covid[[i]], verbose = T)
}

#Saving separate RDS files for each sample:
for (i in names(split.covid)) {
  saveRDS(split.covid[[i]], file = paste0(names(split.covid[i]), ".rds"))
  names(split.covid[i])
}

#After uploading the separate RDS files to Azimuth, incorporating the results
#into the object:
Patient1_critical <- readRDS("Patient1_critical.rds")
Patient1_moderate <- readRDS("Patient1_moderate.rds")
Patient2_critical <- readRDS("Patient2_critical.rds")
Patient2_mild <- readRDS("Patient2_mild.rds")
Patient3_critical <- readRDS("Patient3_critical.rds")
Patient3_severe <- readRDS("Patient3_severe.rds")
Patient3_moderate <- readRDS("Patient3_moderate.rds")
Patient4_critical <- readRDS("Patient4_critical.rds")
Patient4_moderate <- readRDS("Patient4_moderate.rds")
Patient4_severe <- readRDS("Patient4_severe.rds")

#Importing Azimuth's results for each sample:
Patient1_critical$azimuthNames <- read.table("Patient1_critical_azimuth_pred.tsv", sep = "\t", header = T)$predicted.celltype.l2
Patient1_moderate$azimuthNames<- read.table("Patient1_moderate_azimuth_pred.tsv", sep = "\t", header = T)$predicted.celltype.l2
Patient2_critical$azimuthNames<- read.table("Patient2_critical_azimuth_pred.tsv", sep = "\t", header = T)$predicted.celltype.l2
Patient2_mild$azimuthNames<- read.table("Patient2_mild_azimuth_pred.tsv", sep = "\t", header = T)$predicted.celltype.l2
Patient3_critical$azimuthNames<- read.table("Patient3_critical_azimuth_pred.tsv", sep = "\t", header = T)$predicted.celltype.l2
Patient3_severe$azimuthNames<- read.table("Patient3_severe_azimuth_pred.tsv", sep = "\t", header = T)$predicted.celltype.l2
Patient3_moderate$azimuthNames<- read.table("Patient3_moderate_azimuth_pred.tsv", sep = "\t", header = T)$predicted.celltype.l2
Patient4_critical$azimuthNames<- read.table("Patient4_critical_azimuth_pred.tsv", sep = "\t", header = T)$predicted.celltype.l2
Patient4_moderate$azimuthNames<- read.table("Patient4_moderate_azimuth_pred.tsv", sep = "\t", header = T)$predicted.celltype.l2
Patient4_severe$azimuthNames<- read.table("Patient4_severe_azimuth_pred.tsv", sep = "\t", header = T)$predicted.celltype.l2

Patient1_critical@reductions$umap_azimuth <- readRDS("Patient1_critical_azimuth_umap.rds")
Patient1_moderate@reductions$umap_azimuth <- readRDS("Patient1_moderate_azimuth_umap.rds")
Patient2_critical@reductions$umap_azimuth <- readRDS("Patient2_critical_azimuth_umap.rds")
Patient2_mild@reductions$umap_azimuth <- readRDS("Patient2_mild_azimuth_umap.rds")
Patient3_critical@reductions$umap_azimuth <- readRDS("Patient3_critical_azimuth_umap.rds")
Patient3_severe@reductions$umap_azimuth <- readRDS("Patient3_severe_azimuth_umap.rds")
Patient3_moderate@reductions$umap_azimuth <- readRDS("Patient3_moderate_azimuth_umap.rds")
Patient4_critical@reductions$umap_azimuth <- readRDS("Patient4_critical_azimuth_umap.rds")
Patient4_moderate@reductions$umap_azimuth <- readRDS("Patient4_moderate_azimuth_umap.rds")
Patient4_severe@reductions$umap_azimuth <- readRDS("Patient4_severe_azimuth_umap.rds")

gc()

split.covid <- list(Patient1_critical, Patient1_moderate, Patient2_critical, Patient2_mild,
                    Patient3_critical, Patient3_moderate, Patient3_severe,
                    Patient4_critical, Patient4_moderate, Patient4_severe)


#Saving the split object:
saveRDS(split.covid, "02-split-scaled.rds")
remove(Patient1_critical, Patient1_moderate, Patient2_critical, Patient2_mild, 
       Patient3_critical, Patient3_moderate,
       Patient3_severe, Patient4_critical, Patient4_moderate, Patient4_severe, i)

#################################Integration####################################

#Importing the split file:
split.covid <- readRDS("02-split-scaled.rds")

#Selecting the most variable features for integration:
integ.feat <- SelectIntegrationFeatures(object.list = split.covid, nfeatures = 5000)

#Performing canonical correlation analysis (CCA), and identifying and filtering anchors:
anchors <- FindIntegrationAnchors(object.list = split.covid, anchor.features = integ.feat,
                                  verbose = T)

saveRDS(anchors, "anchors.rds")

remove(covid, split.covid, integ.feat)
gc()

#Integrating across conditions:
covid <- IntegrateData(anchorset = anchors) #This step had to be done on the Linux machine,
                                            #as an RScript, due to memory restrictions.

remove(anchors)
gc()

saveRDS(covid, "03-integrated.rds")
covid <- readRDS("03-integrated.rds")

DefaultAssay(covid) <- "integrated"

#Scaling the data:
covid <- ScaleData(covid, features = rownames(covid))

#Linear dimensional reduction:
covid <- RunPCA(covid, features = VariableFeatures(object = covid), verbose = T)
covid <- RunUMAP(covid)
?RunUMAP
#Examining PCA results:
print(covid[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(covid, dims = 1:2, reduction = "pca")
DimPlot(covid, reduction = "pca", split.by = "severity")
DimPlot(covid, reduction = "pca")
DimHeatmap(covid, dims = 1:15, cells = 500, balanced = TRUE)

#Determining the dimensionality of the data:
ElbowPlot(covid) #This implies going for around 18 might be OK.

saveRDS(covid, "03-integrated.rds")

#################################UMAP Clustering################################

covid <- readRDS("03-integrated.rds")
gc()

#Graph-based clustering:
covid <- FindNeighbors(covid, dims = 1:18, verbose = T)
covid <- FindClusters(covid, verbose = T,
                      resolution = 1.0)

Idents(covid) <- "integrated_snn_res.1" #This one seems to make the most sense.

#Non-linear dimensional reduction:
covid <- RunUMAP(covid, dims = 1:18, return.model = T)
DimPlot(covid, reduction = "umap", label = T, repel = T) + NoLegend()
DimPlot(covid, reduction = "umap", split.by = "severity") + NoLegend()
Idents(covid) <- "azimuthNames"
gc()

saveRDS(covid, "04-clustered.rds")

#################################Finding Cluster Bio-markers#####################

covid <- readRDS("04-clustered.rds")

#Finding marker genes:
covid.markers <- FindAllMarkers(covid, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(covid.markers, "covid-markers.rds")
covid.markers <- readRDS("covid-markers.rds")

covid.markers %>%
  group_by(cluster) %>%
 slice_max(n = 10, order_by = avg_log2FC) -> top10

#################################Azimuth Code Clustering########################
covid <- readRDS("04-clustered.rds")

#Loading the reference:
reference <- LoadH5Seurat("pbmc_multimodal.h5seurat")

#Preprocessing:
covid <- SCTransform(
  object = covid,
  verbose = T
)
gc()

saveRDS(covid, "04-clustered.rds")

#Finding anchors between covid object and reference:
transfer.anchors <- FindTransferAnchors(
  reference = reference,
  query = covid,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:40
)
gc()

saveRDS(transfer.anchors, "transfer-anchors.rds")

covid <- MapQuery(
  anchorset = transfer.anchors,
  query = covid,
  reference = reference,
  refdata = list(
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)
gc()

saveRDS(covid, "05-mapped.rds")

covid <- readRDS("05-mapped.rds")

gc()

Idents(covid) <- "azimuthNames"

tiff("graphs/azimuth2.tiff", width = 3237, height = 2176, units = "px", res = 300)
dev.off()
#Changing the order of severity levels:
covid$severity <-  factor(covid$severity, levels = c("mild", "moderate", "severe", "critical"))
covid$sample <- factor(covid$sample, levels = c("moderate272_Patient1", "critical293_Patient1",
                                            "mild227_Patient2", "critical238_Patient2",
                                            "critical119_Patient3", "severe123_Patient3", "moderate138_Patient3",
                                            "critical120_Patient4", "severe122_Patient4", "moderate124_Patient4"))
ggsave(path = "graphs/", filename = "azimuth.jpeg", dpi = 300,
         plot = DimPlot(covid, label = T, repel = T, label.size = 5, reduction = "umap") +
           NoLegend() + ggtitle("Azimuth Clusters") + theme(plot.title = element_text(hjust = 0.5, size = 20)))
ggsave(path = "graphs/", filename = "azimuth-by-severity.jpeg", 
       plot = DimPlot(covid, label = T, split.by = "severity", repel = T, reduction = "umap", label.size = 3) +
         NoLegend() + ggtitle("Clusters by Severity") + theme(plot.title = element_text(hjust = 0.5, size = 20)))
ggsave(path = "graphs/", filename = "azimuth-by-sample.jpeg", 
       plot = DimPlot(covid, label = T, split.by = "sample", repel = T, reduction = "umap", label.size = 2) +
         NoLegend() + ggtitle("Clusters by Sample") + theme(plot.title = element_text(hjust = 0.5, size = 20)))


for (i in levels(factor(covid$sample))) {
    DimPlot(covid[,covid$sample == as.character(i)], reduction = "umap", label = T, repel = T, label.size = 2) +
    NoLegend() + ggtitle(i) + theme(plot.title = element_text(hjust = 0.5, size = 20))
    ggsave(filename = paste0("graphs/", i, ".jpeg"))
} 
for (i in levels(factor(covid$severity))) {
  DimPlot(covid[,covid$severity == as.character(i)], reduction = "umap", label = T, repel = T, label.size = 2) +
    NoLegend() + ggtitle(i) + theme(plot.title = element_text(hjust = 0.5, size = 20))
  ggsave(filename = paste0("graphs/", i, ".jpeg"))
} 

DimPlot(covid, reduction = "umap", split.by = "severity") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20), legend.text=element_text(size=15))

#Cell numbers per patient per severity:
counts <- table(covid$sample, covid$azimuthNames) %>% addmargins()
write.table(counts, "cell-counts-per-sample.tsv", sep = "\t", col.names = NA)

counts.severity <- table(covid$severity, covid$azimuthNames) %>% addmargins()
write.table(counts.severity, "cell-counts-per-severity.tsv", sep = "\t", col.names = NA)


ggplot(data = covid@meta.data[covid$patient == "Patient1",], aes(x = sample, fill = azimuthNames)) + geom_bar(position = "fill") +
    scale_y_continuous(name = "Percentage",
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab("Sample") +
    scale_fill_manual(values =  glasbey(30)) + labs(fill = "Cells") + ggtitle("Patient 1") +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_discrete(labels = c("moderate", "critical")) 
ggsave(path = "graphs/", filename =  "patient1-cells.jpeg")

ggplot(data = covid@meta.data[covid$patient == "Patient2",], aes(x = sample, fill = azimuthNames)) + geom_bar(position = "fill") +
        scale_y_continuous(name = "Percentage",
                           breaks = c(0, 0.25, 0.5, 0.75, 1),
                           labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab("Sample") +
        scale_fill_manual(values =  glasbey(30)) + labs(fill = "Cells") + ggtitle("Patient 2") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_x_discrete(labels = c("mild", "critical")) 
ggsave(path = "graphs/", filename =  "patient2-cells.jpeg")

ggplot(data = covid@meta.data[covid$patient == "Patient3",], aes(x = sample, fill = azimuthNames)) + geom_bar(position = "fill") +
  scale_y_continuous(name = "Percentage",
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab("Sample") +
  scale_fill_manual(values =  glasbey(30)) + labs(fill = "Cells") + ggtitle("Patient 3") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("critical", "severe", "moderate")) 
ggsave(path = "graphs/", filename =  "patient3-cells.jpeg")

ggplot(data = covid@meta.data[covid$patient == "Patient4",], aes(x = sample, fill = azimuthNames)) + geom_bar(position = "fill") +
  scale_y_continuous(name = "Percentage",
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab("Sample") +
  scale_fill_manual(values =  glasbey(30)) + labs(fill = "Cells") + ggtitle("Patient 4") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("critical", "severe", "moderate")) 
ggsave(path = "graphs/", filename =  "patient4-cells.jpeg")

ggplot(data = covid@meta.data, aes(x = sample, fill = azimuthNames)) + geom_bar(position = "fill") +
  scale_y_continuous(name = "Percentage",
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) + xlab("Patient") +
  scale_fill_manual(values =  glasbey(30)) + labs(fill = "Azimuth Names")
  
ggplot(data = covid@meta.data, aes(x = severity, fill = covid$azimuthNames)) + geom_bar(position = "fill") +
  scale_y_continuous(name = "Percentage",
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) +
  xlab("Severity") +
  scale_fill_manual(values =  glasbey(30)) + NoLegend()

#Adding in the Azimuth prediction scores:
files <- list.files(pattern = "Patient[[:digit:]]{1}.*.tsv")
scores <- data.frame(NULL)
for (i in files) {
  
  scores <- rbind(scores, read.table(file = i, header = T, sep = "\t"))
  
}
covid <- AddMetaData(covid, metadata = scores$predicted.celltype.l2.score, col.name = "az.pred.score")

saveRDS(covid, "05-mapped.rds")

#################################Manual Clustering##############################

covid <- readRDS("05-mapped.rds")
#Looking at some known PBMC markers:
FeaturePlot(covid, features = c("MS4A1", "GNLY", "NKG7", "CD3E", 
                                "CD14", "FCER1A", "CST3", "FCGR3A",
                                "MS4A7", "FCGR3A", "LYZ", "PPBP",
                                "CD8A", "IL7R", "S100A4", "CCR7"), label = T)

##Comparing well known cell type markers with our cells:
#Using CD14+ monocyte markers to determine monocytes:
FeaturePlot(covid, reduction = "umap",
            features = c("CD14", "FCGR3A", "FCGR3B", "CD68", "S100A12"),
            order = T, min.cutoff = "q10", label = T)

#Cluster(s): 

#Doing the same for FCGR3A+ monocyte markers:
FeaturePlot(covid, reduction = "umap",
            features = c("FCGR3A", "MS4A7"),
            order = T, min.cutoff = "q10", label = T)
#Cluster(s):

#Memory B cells:
FeaturePlot(covid, reduction = "umap",
            features = c("CD27", "CD38", "IGHG1", "IGHA1"),
            order = T, min.cutoff = "q10", label = T)
#Cluster(s): 29, 33

#Plasmablasts:
FeaturePlot(covid, reduction = "umap",
            features = c("XBP1", "MZB1", "IGHA", "IGHG", "CD38", "MS4A1", "CD20"),
            order = T, min.cutoff = "q10", label = T)
#Cluster(s): 29 

#Macrophages:
FeaturePlot(covid, reduction = "umap",
            features = c("MARCO", "ITGAM", "ADGRE1", "CCR5", "CD68"),
            order = T, min.cutoff = "q10", label = T)
#Cluster(s):

#Conventional dendritic cells:
FeaturePlot(covid, reduction = "umap", 
            features = c("FCER1A", "CST3"),
            order = T, min.cutoff = "q10", label = T)
#Cluster(s):

#Plasmacytoid dendritic cells:
FeaturePlot(covid, reduction = "umap",
            features = c("IL3RA", "GZMB", "SERPINF1", "ITM2C"),
            order = T, min.cutoff = "q10", label = T)
#Cluster(s):

#B cells:
FeaturePlot(covid, reduction = "umap",
            features = c("CD79A", "MS4A1"),
            order = T, min.cutoff = "q10", label = T)
#Cluster(s):

#T cells:
FeaturePlot(covid, reduction = "umap",
            features = "CD3D",
            order = T, min.cutoff = "q10", label = T)
#Cluster(s):

#Helper 1 CD4+ T cells:
FeaturePlot(covid, reduction = "umap",
            features = c("CD3D", "CD4", "CXCR3", "IL7R", "CCR7", "IFNG"),
            order = T, min.cutoff = "q10", label = T)
#Cluster(s): 2, 3, 6, 8, 12

#Helper 2 CD4+ T cells:
FeaturePlot(covid, reduction = "umap",
            features = c("IL3", "IL4", "IL5", "IL6", "IL10", "IL13"),
            order = T, min.cutoff = "q10", label = T)
#Cluster(s): 2, 8, 16

#CD8+ T cells:
FeaturePlot(covid, reduction = "umap",
            features = c("CD3D", "CD8A"),
            order = T, min.cutoff = "q10", label = T)
#Cluster(s):

#CD8+ naive T cells:
FeaturePlot(covid, reduction = "umap",
            features = c("CD3E", "CD4", "CD8A", "CD8B", "CCR7", "LEF1", "TCF7",
                         "CD27", "CD28", "SELL", "S1PR1"),
            order = T, min.cutoff = "q10", label = T)
#Cluster(s):

#CD8+ T cells:
FeaturePlot(covid, reduction = "umap",
            features = c("CD3D", "CD8A"),
            order = T, min.cutoff = "q10", label = T)
#Cluster(s):


#NK cells:
FeaturePlot(covid, reduction = "umap",
            features = c("GNLY", "NKG7"),
            order = T, min.cutoff = "q10", label = T)
#Cluster(s):

#Megakaryocytes:
FeaturePlot(covid, reduction = "umap",
            features = "PPBP",
            order = T, min.cutoff = "q10", label = T)
#Cluster(s):

#Erythrocytes:
FeaturePlot(covid, reduction = "umap",
            features = c("HBB", "HBA2"),
            order = T, min.cutoff = "q10", label = T)
#Cluster(s):
