#!/usr/bin/env Rscript
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: promoters_targets_Goggle.R
##
## Description: 
#~
## This tool will run a R script for the extraction of promoter sequences of the targets of master regulators previously identified with VIPER
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
## custom script using VIPER results previously generated
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library("knitr")
library("tidyverse")
library("rmarkdown")
library("GenomeInfoDb")
library("Rsamtools")
library("GenomicAlignments")
library("BiocParallel")
library("Rsubread")
library("GenomicFeatures")
library("GenomicRanges")
library("BSgenome.Hsapiens.UCSC.hg38")
library('rtracklayer')


#
reference_folder <- ".../reference_folder/"
ARACNE <- paste0(Output_folder,"/ARACNE_output/")
Output_folder <- "3.GRN/"

SAVE_IN <- paste0(Output_folder,"/Targets/")
dir.create(SAVE_IN)
#

# Import network from ARCACNE
network <- read.delim(paste0(ARACNE,"/","network_viper.txt"),sep="\t",header=FALSE,stringsAsFactors=FALSE,check.names=FALSE) 
colnames(network) <- c("TF","Target","Weight","Weight2")


# Import the reference file - try to perform the same using Ensambl data GTFs
Ref_files <- list.files(reference_folder,full.names = TRUE, recursive = FALSE,ignore.case = FALSE, include.dirs = TRUE)
Ref_Txdb <- loadDb(Ref_files[grep(".sqlite",Ref_files)])
seqlevels(Ref_Txdb) <- seqlevels(Ref_Txdb)[grep("random|alt|chrUn|fix",seqlevels(Ref_Txdb),invert=T)]

# Get promoters
PR <- promoters(Ref_Txdb, upstream=500, downstream=50) #TSS-proximal - 500bp + 50bp

# Loop thorugh individual Master regulatros from VIPER
name <- c("FOXA2","HNF1B","MYRF","ZEB1")

Master_reg <- unique(network$TF)
Master_reg <- Master_reg[which(Master_reg %in% name)]

for(x in seq_along(Master_reg)){
    nm <-Master_reg[x]
    cat("processing ",nm,"\n")

    # get the targets for the master regulator
    target_file <- network[which(network$TF == nm),]
    
    targets <- unique(target_file$Target)
    cat ("Number of targets: ",length(targets),"\n")
    
    PR_ALL <- PR[which(names(PR) %in% targets)]
    mcols(PR_ALL)$name <- names(PR_ALL)  # Add gene name as metadata
    
    #export to bed 
    rtracklayer::export.bed(PR_ALL, paste0(SAVE_IN, "/", nm, "_promoters.bed"))
}






