#!/usr/bin/env Rscript
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: scRNAseq_step3.R
##
## Description: 
#~
## This script will cluster and score tumor cells after harmonization
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
##
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
Output_folder <- ".../scRNAseq_output/"
signature <- ".../Gene_Signature/"

SAVE_IN <- paste0(Output_folder,"/Step3/")
dir.create(SAVE_IN)


##--- Tumor cells with integration ----###
integrated <- read_rds(paste0(Output_folder,"/Step2/Seurat_Harmony.rds"))

Data_tumor <- integrated
#norm.counts <- GetAssayData(object = Data_tumor, slot = "data")
#write.table(norm.counts, file=paste0(SAVE_IN,"/norm.count.matrix_Harmony.txt"), quote=F, sep="\t") 


###----- Gene Signature Score -----###
## AddModuleScore: Calculate module scores for feature expression programs 
#Calculate the average expression levels of each program subtracted by the aggregated expression of control feature sets. 
#All analyzed features are binned based on averaged expression, and the control features are randomly selected from each bin.
dir.create(paste0(SAVE_IN,"/Score/"))

files <- list.files(signature,full.names = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = TRUE)
files  <- gsub("//","/",files)
files <- files[c(1)]  #neuronal_signature

for(sig in files){
	name <- gsub(".txt","",basename(sig))
	#dir.create(paste0(SAVE_IN,"/",name))
	cat("processing Gene Signature: ",name,"\n")
	DF <- read.delim(sig,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
	for(ll in unique(DF[,"cell_type"])){
		cat(" Processing :",ll,"\n")
    	features <- list(c(DF[DF[,"cell_type"] == ll,"Gene"]))
		common <- intersect(unlist(features), rownames(Data_tumor))
		if(length(common) > 10) {ctrl=10} else{ctrl=length(common)}  
    	Data_tumor <- AddModuleScore(object = Data_tumor,features = features, ctrl = ctrl,name = paste0(name,"_",ll))   #nbin=24
		Data_tumor[[paste0(name,"_",ll,".Scaled")]] <- scale(Data_tumor[[paste0(name,"_",ll,"1")]],center=TRUE,scale=TRUE)
  	}
}

pdf(paste0(SAVE_IN,"/Score/","/violinPlot_Score_Neuronal.pdf"),width=15,height=15)
plot1 <- VlnPlot(Data_tumor, group.by = "Sample.ID", 
			pt.size = 0,  #adjust=1.5  
			features = c("neuronal_signature_SynapticTransmission1","neuronal_signature_NeuronalDifferentiation1")) &
			geom_violin(colour="white",draw_quantiles = 0.5) 
print(plot1)
dev.off()


##--- Clustering the tumor data ----###
dir.create(paste0(SAVE_IN,"/Clustering/"))

set.seed(12345)
co=10
Data_tumor <- FindNeighbors(Data_tumor, reduction = "harmony", dims = 1:co, k.param = 50, verbose = FALSE)      
Data_tumor <- FindClusters(Data_tumor, resolution=seq(0.1,0.8,0.1), algorithm = 4, verbose = FALSE)   # resolution=seq(0.2,1.2,0.2)

resolution <- c(paste0("SCT_snn_res.",seq(0.1,0.8,0.1)))

pdf(paste0(SAVE_IN,"/Clustering/","/UMAP_tumor_seq.pdf")) 
lapply(X = resolution, FUN = function(x) {
DimPlot(Data_tumor,reduction = "umap", group.by = x, label = TRUE,repel = TRUE) + NoLegend() 
})
dev.off()

pdf(paste0(SAVE_IN,"/Clustering/","/UMAP_tumor_split_seq.pdf")) 
lapply(X = resolution, FUN = function(x) {
DimPlot(Data_tumor,reduction = "umap", group.by = x, split.by = x, label = TRUE,repel = TRUE) + NoLegend()   
})
dev.off()

pdf(paste0(SAVE_IN,"/Clustering/","/violinPlot_Score_Neuronal_seq.pdf"),width=15,height=15)
lapply(X = resolution, FUN = function(x) {
VlnPlot(Data_tumor, same.y.lims = TRUE, group.by = x, 
			pt.size = 0,  #adjust=1.5
		 	features = c("neuronal_signature_SynapticTransmission1","neuronal_signature_NeuronalDifferentiation1")) &
			scale_y_continuous(limits =c(-0.10,0.25)) &
		 	geom_violin(colour="white",draw_quantiles = 0.5) 
})
dev.off()


