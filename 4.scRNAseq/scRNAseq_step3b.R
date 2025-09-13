#!/usr/bin/env Rscript
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: scRNAseq_step3b.R
##
## Description: 
#~
## This script will score tumor cells using alternative scoring methods
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
Output_folder <- "4.scRNAseq/"
signature <- ".../Gene_Signature/"

SAVE_IN <- paste0(Output_folder,"/Step3/")
dir.create(SAVE_IN)


##--- Tumor cells with integration ----###
integrated <- read_rds(paste0(Output_folder,"/Step2/Seurat_Harmony.rds"))

Data_tumor <- integrated
#norm.counts <- GetAssayData(object = Data_tumor, slot = "data")
#write.table(norm.counts, file=paste0(SAVE_IN,"/norm.count.matrix_Harmony.txt"), quote=F, sep="\t") 


###----- Gene Signature Score -----###
#GSVA package: GSVA or ssGSEA as enrichment score
library(GSVA)

files <- list.files(signature,full.names = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = TRUE)
files  <- gsub("//","/",files)
files <- files[c(1)]  #neuronal_signature

neuronal_DF = read.delim(files,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)

Synaptic_list <- neuronal_DF[neuronal_DF$cell_type == "SynapticTransmission","Gene"]
NeuronalDifferentiation_list <- neuronal_DF[neuronal_DF$cell_type == "NeuronalDifferentiation","Gene"]
neuronal_list <- list(SynapticTransmission=Synaptic_list,NeuronalDifferentiation=NeuronalDifferentiation_list)

GSVA.table <- as.data.frame(Data_tumor@assays$SCT@data[,])
GSVA.table <- as.matrix(GSVA.table)

gsva <- gsva(GSVA.table, neuronal_list, method="ssgsea", min.sz=3, verbose=FALSE)
gsva.score <- as.data.frame(gsva)
gsva.score <- as.data.frame(t(gsva.score))

Data_tumor <- AddMetaData(object = Data_tumor, metadata = gsva.score, col.name = c('SynapticTransmission','NeuronalDifferentiation'))

pdf(paste0(SAVE_IN,"/Score/","/violinPlot_Score_Neuronal_GSVA.pdf"),width=15,height=15)
plot2 <- VlnPlot(Data_tumor, same.y.lims = TRUE, group.by = "SCT_snn_res.0.3",
			pt.size = 0,  #adjust=1.5
		 	features = c("SynapticTransmission","NeuronalDifferentiation")) &
		 	geom_violin(colour="white",draw_quantiles = 0.5) 
print(plot1)		
dev.off()

