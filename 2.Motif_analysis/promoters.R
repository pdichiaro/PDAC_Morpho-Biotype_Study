#!/usr/bin/env Rscript
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: promoters.R
##
## Description: 
#~
## This script will extract promoters in fasta format using gene expression data (LMD-seq) 
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
Output_folder <- ".../Output/"

SAVE_IN <- paste0(Output_folder,"/PROMOTERS/")
dir.create(SAVE_IN)

#
DEG_genes <- read.delim("/.../DEG_folder/DEG_list.txt",sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
rownames(DEG_genes) <- DEG_genes[,1]
colnames(DEG_genes) <- c("Gene","cluster")
#

# Import the reference file - try to perform the same using Ensambl data GTFs
Ref_files <- list.files(reference_folder,full.names = TRUE, recursive = FALSE,ignore.case = FALSE, include.dirs = TRUE)
Ref_Txdb <- loadDb(Ref_files[grep(".sqlite",Ref_files)])
seqlevels(Ref_Txdb) <- seqlevels(Ref_Txdb)[grep("random|alt|chrUn|fix",seqlevels(Ref_Txdb),invert=T)]

Anno_Txdb <- read.delim(Ref_files[grep("_merged_annotation.txt",Ref_files)])
Keep <- c("protein_coding","protein-coding", "protein-coding gene")
ww <- which(Anno_Txdb$type_of_gene %in% Keep)

PR <- promoters(Ref_Txdb, upstream=500, downstream=50)

#PR_BK <- PR[which(names(PR) %in% Anno_Txdb[ww,"Symbol"])]
#promoter_seq <- getSeq(Hsapiens, PR_BK)  # get sequences for all the promoters 
#writeXStringSet(promoter_seq, paste0(SAVE_IN,"/BK_promoters_500_50.fasta"), append=FALSE,compress=FALSE,compression_level=NA, format="fasta")

PR_BK <- PR[which(names(PR) %in% DEG_genes$Gene)]
promoter_seq <- getSeq(Hsapiens, PR_BK)  # get sequences for the promoters of DEGs
writeXStringSet(promoter_seq, paste0(SAVE_IN,"/BK_DEG_promoters_500_50.fasta"), append=FALSE,compress=FALSE,compression_level=NA, format="fasta")


# Loop thorugh individual stat file separate UP/DOWN p<=0.01 ABS(LOG2FC)>=2
DEG_folder <- ".../DEG_folder/"
DEG <- list.files(DEG_folder,pattern="DEG_RESULT.txt",full.names=TRUE)
LOOG_DEG <- DEG[grep("GL_vs_Rest|TR_vs_Rest|UN_vs_Rest",DEG)] # MODIFY
LOOG_DEG <- LOOG_DEG[order(LOOG_DEG,decreasing=T)]

for(x in seq_along(LOOG_DEG)){
    nm <- gsub("/","",gsub("_DEG_RESULT.txt","",gsub(DEG_folder,"",LOOG_DEG[x])))
    cat("Import DEGS ",nm,"\n")
    name <- gsub("_vs_Rest","", nm)
    file=read.delim(LOOG_DEG[x],sep="\t",row.names=1)
    file <- file[order(as.numeric(as.character(file$padj)),decreasing=F),]
    rownames(file) <- file$Symbol

    sel_up <- file[which(as.numeric(as.character(file$log2FoldChange)) >= log2(2) & as.numeric(as.character(file$padj)) <= 0.01),]
    #sel_down <- file[which(as.numeric(as.character(file$log2FoldChange)) <= log2(1/2) & as.numeric(as.character(file$padj)) <= 0.01),]
    
    genes <- rownames(sel_up)
    PR_DEG <- PR[which(names(PR) %in% genes)]

    # get sequences for the promoters 
    promoter_seq <- getSeq(Hsapiens, PR_DEG)

    writeXStringSet(promoter_seq, paste0(SAVE_IN,"/Group_",name,"_promoters_500_50.fa"), append=FALSE, compress=FALSE,compression_level=NA, format="fasta")
}






