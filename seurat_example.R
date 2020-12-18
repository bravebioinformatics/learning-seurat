R.Version()
getwd()
setwd("/Users/labmac1/haramari/seurat")
### Seurat tutorial using PBMC from 10X Genomics

##installing seurat packages
install.packages("Seurat")
install.packages("jlmelville/uwot")
install.packages("dplyr")
install.packages("patchwork")

library(dplyr)
library(Seurat)
library(patchwork)

#Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
#Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

#Examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

dense.size <- object.size(as.matrix(pbmc.data))
dense.size

sparse.size <- object.size(pbmc.data)
sparse.size

dense.size/sparse.size

#The [[ operator can add columns to object metadata. This is a great place to stash QC stats]]
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#FeatureScatter is typically used to visualize feature-feature relationships, but can be used
#for anything calculated by the object, i.e. columns in object etadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

##Normalization
#Normalizing the data using "LogNormalize"
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- NormalizeData(pbmc)

##Feature selection
#Identification of highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly vaiable genes 
top10 <- head(VariableFeatures(pbmc), 10)

#plot variable features with and without labels
var_plot1 <- VariableFeaturePlot(pbmc)
var_plot2 <- LabelPoints(plot = var_plot1, points = top10, repel = TRUE)
var_plot1 + var_plot2
