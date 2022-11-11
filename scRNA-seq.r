#################################Loading########################################
#Loading libraries:
library(Seurat)
library(SeuratDisk)
library(Azimuth)
library(tidyverse)
library(patchwork)
library(Matrix)
library(pals)
library(cowplot)

#Loading data:
mild227 <- Read10X("from_cellranger/mild227/sample_feature_bc_matrix")
moderate124 <- Read10X("from_cellranger/moderate124/sample_feature_bc_matrix")
mild186 <- Read10X("from_cellranger/mild186/sample_feature_bc_matrix")
critical213<- Read10X("from_cellranger/critical213/sample_feature_bc_matrix")
moderate303 <- Read10X("from_cellranger/moderate303/sample_feature_bc_matrix")
critical308 <- Read10X("from_cellranger/critical308/sample_feature_bc_matrix")
moderate138 <- Read10X("from_cellranger/moderate138/sample_feature_bc_matrix")
moderate272 <- Read10X("from_cellranger/moderate272/sample_feature_bc_matrix")
critical119 <- Read10X("from_cellranger/critical119/sample_feature_bc_matrix")
critical120 <- Read10X("from_cellranger/critical120/sample_feature_bc_matrix")
critical238 <- Read10X("from_cellranger/critical238/sample_feature_bc_matrix")
critical293 <- Read10X("from_cellranger/critical293/sample_feature_bc_matrix")
healthy1 <- Read10X("from_cellranger/healthy1/sample_feature_bc_matrix")
healthy2 <- Read10X("from_cellranger/healthy2/sample_feature_bc_matrix")
healthy3 <- Read10X("from_cellranger/healthy3/sample_feature_bc_matrix")

#Creating Seurat objects:
mild227 <- CreateSeuratObject(counts = mild227, min.cells = 3, min.features = 100)
mild186 <- CreateSeuratObject(counts = mild186, min.cells = 3, min.features = 100)
moderate124 <- CreateSeuratObject(counts = moderate124, min.cells = 3, min.features = 100)
moderate138 <- CreateSeuratObject(counts = moderate138, min.cells = 3, min.features = 100)
moderate272 <- CreateSeuratObject(counts = moderate272, min.cells = 3, min.features = 100)
moderate303 <- CreateSeuratObject(counts = moderate303, min.cells = 3, min.features = 100)
critical119 <- CreateSeuratObject(counts = critical119, min.cells = 3, min.features = 100)
critical120 <- CreateSeuratObject(counts = critical120, min.cells = 3, min.features = 100)
critical213 <- CreateSeuratObject(counts = critical213, min.cells = 3, min.features = 100)
critical238 <- CreateSeuratObject(counts = critical238, min.cells = 3, min.features = 100)
critical293 <- CreateSeuratObject(counts = critical293, min.cells = 3, min.features = 100)
critical308 <- CreateSeuratObject(counts = critical308, min.cells = 3, min.features = 100)
healthy1 <- CreateSeuratObject(counts = healthy1, min.cells = 3, min.features = 100)
healthy2 <- CreateSeuratObject(counts = healthy2, min.cells = 3, min.features = 100)
healthy3 <- CreateSeuratObject(counts = healthy3, min.cells = 3, min.features = 100)

covid <- merge(x = critical119, y = c(critical120, critical213, critical238,
                                      critical293, critical308,
                                      mild186, mild227, 
                                      moderate124, moderate138, moderate272, moderate303,
                                      healthy1, healthy2, healthy3), 
               add.cell.ids = c("critical119_Patient5", "critical120_Patient6", "critical213_Patient3",
                                "critical238_Patient4", "critical293_Patient1", "critical308_Patient2", 
                                "mild186_Patient3", "mild227_Patient4", 
                                "moderate124_Patient6", "moderate138_Patient5", "moderate272_Patient1", "moderate303_Patient2",
                                "healthy1_control1", "healthy2_control2", "healthy3_control3"),
               project = "covid-19")

head(covid@meta.data)

#Saving file for upload to Azimuth:
saveRDS(covid, "00-covid-raw.rds")

remove(critical119, critical120, critical238, critical213, critical293, critical308,
      mild186, mild227, moderate124, moderate138, moderate272, moderate303,
      healthy1, healthy2, healthy3)
gc()

#################################Environment####################################
#Preparing the environment:
#Adjusting the limit for allowable R object sizes: 
options(future.globals.maxSize = 9000 * 1024^2)

#Clearing memory:
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
head(metadata[which(str_detect(metadata$cell, ".*healthy")),])

#Creating a severity column:
metadata$severity <- NA
metadata$severity[which(str_detect(metadata$cell, ".*mild"))] <- "mild"
metadata$severity[which(str_detect(metadata$cell, ".*moderate"))] <- "moderate"
metadata$severity[which(str_detect(metadata$cell, ".*critical"))] <- "critical"
metadata$severity[which(str_detect(metadata$cell, ".*healthy"))] <- "healthy"
metadata$severity <- as.factor(metadata$severity)
#Creating a sample column:
metadata$sample <- NA
metadata$sample <- sub("(.*?)_{1}(.*?)($|-.*)", "\\1", rownames(metadata))
#Adjusting the sample variable:
metadata$sample <- factor(metadata$sample, levels = c("healthy1_control1", "healthy2_control2", "healthy3_control3",
                                            "moderate272_Patient1", "critical293_Patient1",
                                            "moderate303_Patient2", "critical308_Patient2",
                                            "mild186_Patient3", "critical213_Patient3",
                                            "mild227_Patient4", "critical238_Patient4",
                                            "critical119_Patient5", "moderate138_Patient5",
                                            "critical120_Patient6", "moderate124_Patient6"))

#Creating a patient column:
metadata$patient <- sub("(.*?)_", "\\2", metadata$sample)
metadata$patient[metadata$patient == "control1"] <- "Control 1"
metadata$patient[metadata$patient == "control2"] <- "Control 2"
metadata$patient[metadata$patient == "control3"] <- "Control 3"
metadata$patient[metadata$patient == "Patient1"] <- "Patient 1"
metadata$patient[metadata$patient == "Patient2"] <- "Patient 2"
metadata$patient[metadata$patient == "Patient3"] <- "Patient 3"
metadata$patient[metadata$patient == "Patient4"] <- "Patient 4"
metadata$patient[metadata$patient == "Patient5"] <- "Patient 5"
metadata$patient[metadata$patient == "Patient6"] <- "Patient 6"

#Creating an outcome column:
metadata$outcome <- "Recovered"
metadata$outcome[metadata$severity == "healthy"] <- "Healthy"
metadata$outcome[metadata$patient == "Patient1" | metadata$patient == "Patient2" | metadata$patient == "Patient3"] <- "Deceased"

#Creating a sample collection date column:
metadata$date[metadata$sample == "moderate272_Patient1"] <- "2021-05-24" 
metadata$date[metadata$sample == "critical293_Patient1"] <- "2021-06-03"
metadata$date[metadata$sample == "moderate303_Patient2"] <- "2021-06-09"
metadata$date[metadata$sample == "critical308_Patient2"] <- "2021-06-18"
metadata$date[metadata$sample == "mild186_Patient3"] <- "2021-03-29"
metadata$date[metadata$sample == "critical213_Patient3"] <- "2021-04-15"
metadata$date[metadata$sample == "mild227_Patient4"] <- "2021-04-22"
metadata$date[metadata$sample == "critical238_Patient4"] <- "2021-04-30"
metadata$date[metadata$sample == "critical213_Patient3"] <- "2021-04-15"
metadata$date[metadata$sample == "critical119_Patient5"] <- "2020-11-27"
metadata$date[metadata$sample == "moderate138_Patient5"] <- "2020-12-25"
metadata$date[metadata$sample == "critical213_Patient3"] <- "2021-04-15"
metadata$date[metadata$sample == "critical120_Patient6"] <- "2020-11-27"
metadata$date[metadata$sample == "moderate124_Patient6"] <- "2020-12-11"
metadata$date[metadata$sample == "healthy1_control1"] <- "2021-07-19"
metadata$date[metadata$sample == "healthy2_control2"] <- "2021-09-01"
metadata$date[metadata$sample == "healthy3_control3"] <- "2021-09-02"

metadata$date <- metadata$date %>% as.Date(origin='1970-01-01')

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
LabelPoints(plot = VariableFeaturePlot(covid[,covid$severity == "healthy"]), points = top20, repel = T)

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

#Saving separate RDS files for each sample to upload to Azimuth:
for (i in names(split.covid)) {
  saveRDS(split.covid[[i]], file = paste0(names(split.covid[i]), ".rds"))
  names(split.covid[i])
}

#After uploading the separate RDS files to Azimuth, incorporating the results
#by importing Azimuth's results for each sample:
for (i in names(split.covid)) {
	split.covid[[i]]$azimuthNames <- read.table(paste0(i, "_azimuth_pred.tsv"), sep = "\t", header = T)$predicted.celltype.l2
}

gc()

#Saving the split object:
saveRDS(split.covid, "02-covid-split-scaled.rds")
remove(critical119_Patient5, critical120_Patient6, critical238_Patient4, critical293_Patient1,
       critical213_Patient3, critical308_Patient2,
       mild227_Patient4, mild186_Patient3,
       moderate124_Patient6, moderate138_Patient5, moderate272_Patient1, moderate303_Patient2,
       healthy1_control1, healthy2_control2, healthy3_control3)

gc()

#################################Integration####################################

#Importing the split file:
split.covid <- readRDS("02-split-scaled.rds")

#Selecting the most variable features for integration:
integ.feat <- SelectIntegrationFeatures(object.list = split.covid, nfeatures = 5000)

#Performing canonical correlation analysis (CCA), and identifying and filtering anchors:
anchors <- FindIntegrationAnchors(object.list = split.covid, anchor.features = integ.feat,
                                  verbose = T)                  #
                                                                #These two steps had to be done
saveRDS(anchors, "anchors.rds")                                 #on the Linux machine,
                                                                #as an RScript,
                                                                #due to memory restrictions.
remove(covid, split.covid, integ.feat)                          ###(They need a ton of RAM!)
gc()                                                            
#
                                                                #
anchors <- readRDS("anchors.rds")                               #
                                                              #  
#Integrating across conditions:                             #    
covid <- IntegrateData(anchorset = anchors, verbose = T)###             

remove(anchors)
gc()

saveRDS(covid, "03-covid-integrated.rds")
covid <- readRDS("03-covid-integrated.rds")

DefaultAssay(covid) <- "integrated"

#Scaling the data:
covid <- ScaleData(covid, features = rownames(covid))

#Linear dimensional reduction:
covid <- RunPCA(covid, features = VariableFeatures(object = covid), verbose = T)

#Examining PCA results:
print(covid[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(covid, dims = 1:2, reduction = "pca")
DimPlot(covid, reduction = "pca", split.by = "severity")
DimPlot(covid, reduction = "pca")
DimHeatmap(covid, dims = 1:15, cells = 500, balanced = TRUE)

#Determining the dimensionality of the data:
ElbowPlot(covid, ndims = 40) #This implies going for around 30 might be OK.

saveRDS(covid, "03-covid-integrated.rds")

#################################UMAP Clustering################################

covid <- readRDS("03-covid-integrated.rds")
gc()

#Graph-based clustering:
covid <- FindNeighbors(covid, dims = 1:40, verbose = T)
covid <- FindClusters(covid, verbose = T,
                      resolution = 0.8)

Idents(covid) <- "integrated_snn_res.0.8" #This one seems to make the most sense.

#Non-linear dimensional reduction:
covid <- RunUMAP(covid, dims = 1:40, return.model = T)
Idents(covid) <- "integrated_snn_res.0.8"
DimPlot(covid, reduction = "umap", label = T, repel = T, raster = T, group.by = "azimuthNames") + NoLegend()
DimPlot(covid, reduction = "umap", split.by = "severity") + NoLegend()
DimPlot(covid, reduction = "umap", split.by = "sample") + NoLegend()
gc()

saveRDS(covid, "04-covid-clustered.rds")
