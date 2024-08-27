#!/usr/bin/env Rscript
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: scRNAseq_step2.R
##
## Description: 
#~
## This script will and harmonize multiple scRNAseq datasets using Harmony
##
## Authors: 
#~
## Pierluigi Di Chiaro
##
## License: 
#~
## GNU GPL v3
## Copyright 2022-2024 
## Copyright Pierluigi Di Chiaro
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Notes:
#~
## h5 files  downloaded from EGA archive (EGAC00001000710)
## This script requires umap-learn python package to run RunUMAP() with umap.method="umap-learn"
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



library(Seurat)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ComplexHeatmap) 
library(RColorBrewer)
library(viridis)
library(tidyverse)
library(data.table)
library(sctransform)
library(reshape2)
library(grid)
library(gridExtra)
library(harmony)

#  
dataset_folder <- "/.../EGAD00010001811/" 
Output_folder <- "4.scRNAseq/"

SAVE_IN <- paste0(Output_folder,"/Step2/")
dir.create(SAVE_IN)


##--- dataset ----###
#Cluster_Data <- read_rds(paste0(Output_folder,"/Step1/seurat1.rds"))
#Data.list <- SplitObject(Cluster_Data, split.by = "Sample.ID")
#13 PDAC samples
path <- list.files(dataset_folder,pattern="EGA*",full.names=TRUE) 
h5_files <- list.files(path,pattern=paste0("*_filtered_gene_bc_matrices_h5.h5","$"),full.names=TRUE) 

samples <- sub("_filtered_gene_bc_matrices_h5.h5","",basename(h5_files))
samples <- gsub("_",".",samples)

# create seurat object:
h5_read <- lapply(h5_files, Read10X_h5)
h5_seurat <- lapply(h5_read, function(x) CreateSeuratObject(x,min.cells = 3, min.features = 200))
Data <- merge(h5_seurat[[1]], y = h5_seurat[2:length(h5_seurat)], add.cell.ids = samples, project = "scRNAseq_Chan")

sample_metadata <- as.data.frame(Data[["orig.ident"]],drop=F)
sample_metadata$Sample.ID <- sub("\\_.*","",rownames(sample_metadata))
metadata <- sample_metadata[,-1] 
Data <- AddMetaData(object = Data, metadata = metadata, col.name = 'Sample.ID')

Data[["orig.ident"]] <- Data[["Sample.ID"]] 
Data[["percent.mt"]] <- PercentageFeatureSet(Data, pattern = "^MT-")
Idents(Data) <- "orig.ident"

#Removing cells with unique feature counts less than 1000 and <25% mitochondrial counts
Data <- subset(Data, subset = nFeature_RNA > 1000 & percent.mt < 25)   #downsample = 1000 
Data <- subset(Data, idents = c("COMP.0158.P","90209.CMP"), invert=TRUE)

# Remove MT-genes from count.data
MT.index <- grep(pattern = "^MT-", x = rownames(Data), value = FALSE) 
Data <- Data[-MT.index, ]

cell_type <- read.delim(paste0(Output_folder,"/Step1","/sample_metadata_celltype.txt"),sep="\t",stringsAsFactors=FALSE,check.names=FALSE)  
Data <- AddMetaData(object = Data, metadata = cell_type, col.name = 'celltype')
Idents(Data) = "celltype"

#Subset tumor cells 
Data <- subset(x = Data, idents = c("malignant"),invert = FALSE)  

# split the dataset into a list of seurat objects
Data.list <- SplitObject(Data, split.by = "Sample.ID")


##--- Normalizing the data ----###
# SCTtransform for each dataset independently
#set.seed(1448145)
#Data.list <- lapply(X = Data.list, FUN = SCTransform, verbose = FALSE, variable.features.n = 3000, seed.use = 1448145)  #vars.to.regress = c("percent.mt") --> remove confounding sources of variation

uncorrected <- merge(Data.list[[1]], y=c(Data.list[-1]))
uncorrected <- SCTransform(uncorrected, verbose = FALSE, variable.features.n = 3000, seed.use = 1448145)

##--- Integrating the data ----###
dir.create(paste0(SAVE_IN,"/Integration/"))

##PCA
set.seed(12345)
uncorrected <- RunPCA(uncorrected, features = rownames(uncorrected), seed.use=12345, verbose = FALSE)

pct <- uncorrected[["pca"]]@stdev / sum(uncorrected[["pca"]]@stdev) * 100   # Determine percent of variation associated with each PC
co <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1  # Determine the difference between variation of PC and subsequent PC

##Check the ‘dimensionality’ of the dataset
pdf(paste0(SAVE_IN,"/Integration/","/ElbowPlot.pdf"))
plot <- ElbowPlot(uncorrected, ndims=30) + ggtitle(paste0("uncorrected","    ","PCs: ",co))  
print(plot)
dev.off()

## HARMONY: integration of single cell data
# we embed cells into PCA space, use fuzzy clustering to assign each cell to multiple cluster, calculate a global centroid for each cluster, calculate a correction factor for each dataset based on the centroids and correct each cell with a cell-specific factor (iterative process until covergence)
integrated <- RunHarmony(uncorrected, group.by.vars = c("Sample.ID"), assay.use = "SCT")       #add covariates

##UMAP
set.seed(12345)
uncorrected <- RunUMAP(uncorrected, umap.method = "umap-learn", dims = 1:co, seed.use=12345)     
integrated <- RunUMAP(integrated, reduction = "harmony", umap.method = "umap-learn",dims = 1:co, seed.use=12345)     

pdf(paste0(SAVE_IN,"/Integration/","/UMAP_patientID.pdf"),width=10,height=10)
plot1 <- DimPlot(uncorrected, reduction = "umap", group.by = "Sample.ID", label = TRUE, repel = TRUE) + ggtitle("uncorrected")
plot2 <- DimPlot(integrated, reduction = "umap", group.by = "Sample.ID", label = TRUE, repel = TRUE) + ggtitle("harmony")
print(plot1)
print(plot2)
dev.off()

pdf(paste0(SAVE_IN,"/Integration/","/UMAP_QC.pdf"),width=12,height=12)
plot1 <- FeaturePlot(uncorrected, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
plot2 <- FeaturePlot(integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
print(plot1)
print(plot2)
dev.off()

write_rds(integrated, paste0(SAVE_IN,"Seurat_Harmony.rds"))

