#!/usr/bin/env Rscript
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: scRNAseq_step1.R
##
## Description: 
#~
## This script will process scRNAseq data for each patient 
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
library(RColorBrewer)
library(viridis)
library(tidyverse)
library(data.table)
library(sctransform)
library(reshape2)
library(grid)
library(gridExtra)
library(infercnv)

#  
dataset_folder <- "/.../EGAD00010001811/" 
Output_folder <- ".../scRNAseq_output/"

SAVE_IN <- paste0(Output_folder,"/Step1/")
dir.create(SAVE_IN)

##--- QC dataset ----###
#13 PDAC samples
path <- list.files(dataset_folder,pattern="EGA*",full.names=TRUE) 
h5_files <- list.files(path,pattern=paste0("*_filtered_gene_bc_matrices_h5.h5","$"),full.names=TRUE) 

samples <- sub("_filtered_gene_bc_matrices_h5.h5","",basename(h5_files))
samples <- gsub("_",".",samples)

# create seurat object:
h5_read <- lapply(h5_files, Read10X_h5)
h5_seurat <- lapply(h5_read, function(x) CreateSeuratObject(x,min.cells = 3, min.features = 200))
Data  <- merge(h5_seurat[[1]], y = h5_seurat[2:length(h5_seurat)], add.cell.ids = samples, project = "scRNAseq_patients")

sample_metadata <- as.data.frame(Data[["orig.ident"]],drop=F)
sample_metadata$Sample.ID <- sub("\\_.*","",rownames(sample_metadata))
sample_metadata <- sample_metadata[,-1] 
Data <- AddMetaData(object = Data, metadata = sample_metadata, col.name = 'Sample.ID')

Data[["orig.ident"]] <- Data[["Sample.ID"]] 
Data[["percent.mt"]] <- PercentageFeatureSet(Data, pattern = "^MT-")
Idents(Data) <- "orig.ident"

pdf(paste0(SAVE_IN,"/QC_violinPlot.pdf"))
plot <- VlnPlot(Data, group.by="Sample.ID", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)  
print(plot)
dev.off()

#Removing cells with unique feature counts less than 1000 and <25% mitochondrial counts
Data <- subset(Data, subset = nFeature_RNA > 1000 & percent.mt < 25)   #downsample = 1000 
Data <- subset(Data, idents = c("COMP.0158.P","90209.CMP"), invert=TRUE)

# Remove MT-genes from count.data
MT.index <- grep(pattern = "^MT-", x = rownames(Data), value = FALSE) 
Data <- Data[-MT.index, ]

##--- Normalizing the data ----###
dir.create(paste0(SAVE_IN,"/Normalization/"))

# split the dataset into a list of seurat objects
Data.list <- SplitObject(Data, split.by = "Sample.ID")

# SCTtransform for each dataset independently
Data.list <- lapply(X = Data.list, FUN = SCTransform, verbose = FALSE)  #vars.to.regress = c("percent.mt") --> remove confounding sources of variation

##Check cell cycle markers and MT percent as cofounding factors
cc.genes <- read.delim(".../Tirosh_2015_CellCycleMarkers.txt", header=TRUE) 
s.genes <- cc.genes[which(cc.genes$CellCycle=="S"),"Gene"]
g2m.genes <- cc.genes[which(cc.genes$CellCycle=="G2M"),"Gene"]
cc <- lapply(X = Data.list, FUN = CellCycleScoring, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

pdf(paste0(SAVE_IN,"/Normalization/","/RidgePlot_cellCycle.pdf"),width=10,height=10)
lapply(X = cc, FUN = function(x) {
  plot <- RidgePlot(x, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2) + ggtitle(x$Sample.ID)
  print(plot)
})
dev.off()

set.seed(12345)
cc <- lapply(X = cc, FUN = RunPCA, features = c(s.genes, g2m.genes), seed.use=12345, verbose = FALSE)
pdf(paste0(SAVE_IN,"/Normalization/","/PCA_cell_cycle.pdf"))
lapply(X = cc, FUN = function(x) {
  plot <- DimPlot(x, reduction = "pca", label = TRUE, repel = TRUE) + ggtitle(x$Sample.ID)
  print(plot)
})
dev.off()

##PCA
set.seed(12345)
Norm <- lapply(X = cc, FUN = RunPCA, features = rownames(cc), seed.use=12345, verbose = FALSE)     
pdf(paste0(SAVE_IN,"/Normalization/","/PCA.pdf"))
lapply(X = Norm, FUN = function(x) {
  plot <- DimPlot(x,reduction = "pca", group.by = "Sample.ID", cols = c("grey50")) + NoLegend() + ggtitle(x$Sample.ID)  #label = TRUE, repel = TRUE
  print(plot)
})
dev.off()

#Determine the point where the percent change in variation between the consecutive PCs is less than 0.1%.
pct <- lapply(X = Norm, FUN = function(x) {x[["pca"]]@stdev / sum(x[["pca"]]@stdev) * 100 })  # Determine percent of variation associated with each PC
co <- lapply(X = pct, FUN = function(x) {sort(which((x[1:length(x) - 1] - x[2:length(x)]) > 0.1), decreasing = T)[1] + 1})  # Determine the difference between variation of PC and subsequent PC

#Check the ‘dimensionality’ of the dataset
pdf(paste0(SAVE_IN,"/Normalization/","/ElbowPlot.pdf"))
lapply(X = Norm, FUN = function(x) {
  w <- unique(x$Sample.ID)
  plot <- ElbowPlot(x, ndims=20) + ggtitle(paste0(x$Sample.ID,"    ","PCs: ",co[[w]]))  
  print(plot)
})
dev.off()

pdf(paste0(SAVE_IN,"/Normalization/","/DimReductionGenes.pdf"),width=15,height=15)
lapply(X = Norm, FUN = function(x) {
  plot <- VizDimLoadings(x, dims = 1:10, nfeatures = 30, reduction = "pca") + ggtitle(x$Sample.ID)  
  print(plot)
})
dev.off()

pdf(paste0(SAVE_IN,"/Normalization/","/Heatmap_PCs.pdf"),width=10,height=10)
lapply(X = Norm, FUN = function(x) {
  plot <- DimHeatmap(x, dims = 1:10, cells = 500, nfeatures = 30, balanced = TRUE)
  print(plot)
})
dev.off()

##--- Clustering the data ----###
dir.create(paste0(SAVE_IN,"/Clustering/"))

#graph-based clustering approach: it embed cells in a graph structure -KNN graph- with edges drawn between cells with similar feature expression patterns and then attempt to partition this graph into highly interconnected ‘communities’.
set.seed(12345)
Cluster_Data <- lapply(X = Norm, FUN = function(x) {
  w <- unique(x$Sample.ID)
  x <- RunUMAP(x, umap.method = "umap-learn", dims = 1:co[[w]], verbose = FALSE) #10 PCs as default; similar result with higher settings for this parameter!
  x <- FindNeighbors(x, dims = 1:co[[w]], k.param = 20, verbose = FALSE) 
  x <- FindClusters(x, resolution = 0.8, algorithm = 4, verbose = FALSE)  #it groups cells together; resolution parameter sets the ‘granularity’ of the clustering with increased values leading to a greater number of clusters.
})  

pdf(paste0(SAVE_IN,"/Clustering/","/UMAP.pdf"))
lapply(X = Cluster_Data, FUN = function(x) {
  plot <- DimPlot(x, reduction = "umap", label = TRUE, repel=TRUE) + ggtitle(x$Sample.ID)  
  print(plot)
})
dev.off()

pdf(paste0(SAVE_IN,"/Clustering/","/topVarFeatures.pdf"),width=10,height=10)
lapply(X = Cluster_Data, FUN = function(x) {
  plot <- DimHeatmap(x, dims = 1:15, cells = 500, balanced = TRUE)  
  print(plot)
})
dev.off()


##--- Finding cluster biomarkers and Assigning cell type identity ----###
# find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- lapply(X = Cluster_Data, FUN = FindAllMarkers,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc")   #test.use = "wilcox"     #min.pct=min fraction of cells
for(t in names(Cluster_Data)){
  marker.table <- markers[[t]] %>% group_by(cluster) %>% arrange(cluster, desc(avg_log2FC))     
  write.table(marker.table , file=paste0(SAVE_IN, "/Clustering/",t,"_Marker.matrix.txt"), quote=F, sep="\t")
}

top <- lapply(X = markers, FUN = function(x) {x %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)})
pdf(paste0(SAVE_IN,"/Clustering/","/topMarkers_heatmap.pdf"),width=15,height=15)
lapply(X = Cluster_Data, FUN = function(x) {
  w <- unique(x$Sample.ID)
  plot <- DoHeatmap(x, features = top[[w]]$gene) + NoLegend() + ggtitle(x$Sample.ID)  
  print(plot)
})
dev.off()

features1 <- lapply(X = markers, FUN = function(x) {x %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC) %>% pull(gene)})
#KRT19,KRT10,EPCAM --> epithelial;  GCG,INS,SST --> endocrine;  SPP1,CFTR --> exocrine;  CD33,CD68,PTPRC --> immune;  LUM,ACTA2,THY1 --> fibroblast ; FLT1 --> endothelial
features2 = c("KRT19","KRT10","EPCAM", "SPP1","CFTR", "LUM","ACTA2","THY1", "FLT1", "CD33","CD68","PTPRC", "GCG","INS","SST")    

pdf(paste0(SAVE_IN,"/Clustering/","/Markers_violinPlot.pdf"),width=12,height=12)
lapply(X = Cluster_Data, FUN = function(x) {
  w <- unique(x$Sample.ID)
  plot <- VlnPlot(x, features = features1[[w]]) + labs(title = x$Sample.ID)
  print(plot)
})
dev.off()

pdf(paste0(SAVE_IN,"/Clustering/","/Markers_UMAP.pdf"),width=12,height=12)
lapply(X = Cluster_Data, FUN = function(x) {
  w <- unique(x$Sample.ID)
  plot <- FeaturePlot(x, features = features1[[w]]) 
  print(plot)
})
dev.off()

pdf(paste0(SAVE_IN,"/Clustering/","/SelectedMarkers_violinPlot.pdf"),width=12,height=12)
lapply(X = Cluster_Data, FUN = function(x) {
  plot <- VlnPlot(x, features = features2) 
  print(plot)
})
dev.off()

pdf(paste0(SAVE_IN,"/Clustering/","/SelectedMarkers_UMAP.pdf"),width=12,height=12)
lapply(X = Cluster_Data, FUN = function(x) {
  plot <- FeaturePlot(x, features = features2, slot = "data") 
  print(plot)
})
dev.off()


###----- Identify tumor cells by infering CNV -----###
# Infer CNV
dir.create(paste0(SAVE_IN,"/inferCNV/"))

# Assign normal cell types as reference
ref <- list()
for(t in 1:13){
  cluster <- c("11","12","11","8","12","14","15","12","8","13","12","11","9")
  ref[[t]] <- cluster[t]
}
names(ref) <- names(Data.list)

gene_file <- read.delim(".../gencode_v19_gen_pos.txt",sep="\t",stringsAsFactors=FALSE,check.names=FALSE, header = FALSE)  
rownames(gene_file) <- gene_file[,1]
gene_file <- gene_file[,-1]

counts_matrix <- lapply(X = Cluster_Data, FUN = function(x) {
  x = GetAssayData(x, slot="counts")
})

sample_annotation <- lapply(X = Cluster_Data, FUN = function(x) {
  x = as.data.frame(Idents(object = x),drop=F)
})

set.seed(12345)
infercnv_obj <- lapply(X = Cluster_Data, FUN = function(x) {
  w <- unique(x$Sample.ID)
  #dir.create(paste0(SAVE_IN,"/inferCNV/",w))
  counts_matrix <- GetAssayData(x, slot="counts")
  sample_annotation <- as.data.frame(Idents(object = x),drop=F)
  x <- CreateInfercnvObject(raw_counts_matrix=counts_matrix, annotations_file=sample_annotation, delim="\t", gene_order_file=gene_file, ref_group_names=ref[[w]])
  x <- infercnv::run(x,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics  ###Genes with a mean number of counts across cells will be excluded
    out_dir=paste0(SAVE_IN,"/inferCNV/",w),
    cluster_by_groups=TRUE, 
    plot_steps=FALSE,
    denoise=TRUE,
    HMM=TRUE,   #CNV predictions
    no_prelim_plot=TRUE,
    save_rds = FALSE,
    save_final_rds = FALSE,
    num_threads = 6,
    output_format = "pdf")
})  


###----- Assign cell types  -----###
dir.create(paste0(SAVE_IN,"/CellTypeAssignment/"))

cluster.ids <- list(c(rep("malignant",6),"immune",rep("malignant",3),"islets","endothelium","fibroblast"),                                                             #94930
                    c("fibroblast",rep("malignant",2),"fibroblast","duct","fibroblast","fibroblast",rep("malignant",3),"duct","islets","malignant","endothelium"),     #87235
                    c(rep("malignant",8),"duct","immune","islets"),                                                                                                    #91610
                    c(rep("malignant",3),"duct",rep("malignant",2),"fibroblast","islets","endothelium"),                                                               #96460
                    c(rep("malignant",11),"islets","immune","fibroblast"),                                                                                             #91412
                    c(rep("malignant",4),"fibroblast","malignant","fibroblast","fibroblast","duct","malignant","duct","duct","endothelium","islets","fibroblast"),     #G9903
                    c(rep("malignant",11),"fibroblast","endothelium","immune","islets"),                                                                               #100070
                    c(rep("malignant",5),"duct","duct",rep("malignant",2),"fibroblast","duct","immune","endothelium","malignant","fibroblast"),                        #85948
                    c(rep("malignant",7),"fibroblast","malignant","duct"),                                                                                             #87784
                    c(rep("malignant",10),"duct","malignant","immune","endothelium","fibroblast"),                                                                     #95092
                    c(rep("malignant",7),"immune",rep("malignant",2),"duct","islets","duct","fibroblast"),                                                             #95373
                    c(rep("malignant",6),"duct","malignant","immune","fibroblast","islets","endothelium","malignant"),                                                 #91706
                    c(rep("malignant",7),"islets","islets","endothelium","duct",rep("malignant",3)))                                                                   #97727
names(cluster.ids) <- names(Data.list)

Cluster_Data2 <- lapply(X = Cluster_Data, FUN = function(x) {
  w <- unique(x$Sample.ID)
  names(cluster.ids[[w]]) <- levels(x)
  x <- RenameIdents(x,cluster.ids[[w]])  
  #x[["celltype"]] <- Idents(object = x)
})   

pdf(paste0(SAVE_IN,"/CellTypeAssignment/","/UMAP_SinglePatient.pdf"),width=12,height=12)
lapply(X = Cluster_Data2, FUN = function(x) {
  plot <- DimPlot(x, reduction = "umap", group.by = "celltype", label = TRUE, repel=TRUE) + ggtitle(x$Sample.ID) 
print(plot)
})
dev.off()


# Save final R object
Cluster_Data2 <- merge(Cluster_Data2[[1]], y=c(Cluster_Data2[-1]))
Cluster_Data2[["celltype"]] <- Idents(object = Cluster_Data2)

sample_metadata <- as.data.frame(Cluster_Data2[["celltype"]],drop=F)
write.table(sample_metadata, file=paste0(SAVE_IN,"/sample_metadata_celltype.txt"), row.names = TRUE, col.names=TRUE, quote=F,sep="\t") 

#write_rds(Cluster_Data2, paste0(SAVE_IN,"seurat1.rds"))



