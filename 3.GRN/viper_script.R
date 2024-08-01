#!/usr/bin/env Rscript
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: viper_script.R
##
## Description: 
#~
## This tool will run VIPER tool for the identification of master regulators
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
## custom script using DEGs previously generated
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



library("tidyverse")
library("knitr")
library("BiocParallel")
library("ComplexHeatmap")  
library("RColorBrewer")
library("viridis")
library("viper")
library("Rgraphviz")


#
Output_folder <- ".../Output/"
ARACNE <- paste0(Output_folder,"/ARACNE_output/")

SAVE_IN <- paste0(Output_folder,"/viper/")
dir.create(SAVE_IN)


###----- Dataset -----###
MASTER_COUNTS <- read.delim(".../Matrix_rLog_Subgroups.txt",sep="\t",stringsAsFactors=FALSE,check.names=FALSE) 
rownames(MASTER_COUNTS) <- MASTER_COUNTS$Symbol

sample_annotation <- read_tsv(".../sample_annotation.txt")
sample_annotation <- sample_annotation %>% mutate(SUBTYPES=as.character(SUBTYPES))

samples <- colnames(MASTER_COUNTS[,17:length(MASTER_COUNTS)])
samples <- sample_annotation[order(sample_annotation[,"SUBTYPES"],decreasing=TRUE),] %>% pull(Samples)


#########################################################################################################################################
## MASTER REGULATOR ANALYSIS performed by msVIPER 
## msVIPER infers the relative activity of a regulatory gene based on the enrichment of its most closely-regulated targets on a given GES 
#########################################################################################################################################

#Generating the regulon object
network <- paste0(ARACNE,"/network_viper.txt")
regul <- aracne2regulon(network, as.matrix(MASTER_COUNTS[,samples]), format = c("3col"), verbose = FALSE)

#LOOP for each subtype
DEG_folder <- ".../DEG_folder/"
DEG <- list.files(DEG_folder,pattern="DEG_RESULT.txt",full.names=TRUE)

Targets <- unique(as.character(sample_annotation %>% arrange(desc(SUBTYPES)) %>% pull(SUBTYPES)))

for(pp in seq_along(Targets)){
    name <- Targets[pp]
    dir.create(paste0(SAVE_IN,name))
    cat("Processing for: ",name,"\n")
    
    sample_anno <- sample_annotation[order(sample_annotation[,"SUBTYPES"],decreasing=TRUE),] %>% mutate(SUBGROUPS=replace(SUBTYPES, SUBTYPES!=name,"Rest")) %>% 
        mutate(SUBGROUPS=as.character(SUBGROUPS))
    s1 <- sample_anno %>% filter(SUBGROUPS %in% name) %>% pull(Samples)
    s2 <- sample_anno %>% filter(SUBGROUPS %in% "Rest") %>% pull(Samples)

    #Generating the gene expression signatures(GES)
    LOOG_DEG <- DEG[grep(paste0(name,"_vs_Rest"),DEG)] # MODIFY
    cat("Import DEGS ",name,"\n")
    file <- read.delim(LOOG_DEG,sep="\t",row.names=1)
    file <- file[order(as.numeric(as.character(file$padj)),decreasing=F),]
    rownames(file) <- file$Symbol
    sel_up <- file[which(as.numeric(as.character(file$log2FoldChange)) >= log2(2) & as.numeric(as.character(file$padj)) <= 0.01),]
    sel_down <- file[which(as.numeric(as.character(file$log2FoldChange)) <= log2(1/2) & as.numeric(as.character(file$padj)) <= 0.01),]
    sel <- rbind(sel_up,sel_down)

    signature=sel[rownames(sel),"log2FoldChange"]
    names(signature) <- rownames(sel)

    #NULL model by sample permutations
    nullmodel <- ttestNull(as.matrix(MASTER_COUNTS[rownames(sel),s1]), as.matrix(MASTER_COUNTS[rownames(sel),s2]), per = 1000, repos = TRUE, verbose = FALSE)  #rownames(sel)
    
    pdf(paste0(SAVE_IN,name,"/nullmodel.pdf"))
    plot(density(nullmodel))
    dev.off()

    #msVIPER
    mrs <- msviper(signature, regul, nullmodel, verbose = FALSE)
    mrs <- bootstrapmsviper(mrs, "mode")

    pdf(paste0(SAVE_IN,name,"/msVIPER_plot.pdf"))
    plot(mrs, 25, cex = .7)
    dev.off()

    msVIPER <- summary(mrs,100)
    msVIPER$genes <- NA

    #Leading-edge analysis
    mrs <- ledge(mrs)

    for(t in seq_along(msVIPER$Regulon)){
        regulon <- msVIPER$Regulon[t]
        target <- mrs$regulon[[regulon]][["tfmode"]]

        genes <- unique(unlist(mrs$ledge[regulon]))
        target <- target[which(names(target) %in% genes)]
        target <- as.data.frame(target)
        target$gene_target <- rownames(target)
        target$new_target <- paste0(target$gene_target,"_",target$target)

        msVIPER$genes[t] <- paste0(target$new_target,collapse=",")
        #msVIPER$genes[t] <- paste0(genes,collapse=" ")
    }
    write.table(msVIPER,file=paste0(SAVE_IN,name,"/msVIPER.txt"),sep="\t",col.names=NA)

    #Shadow analysis
    mrs <- shadow(mrs, regulators = 100, verbose = FALSE)
    
    shadow.pair <- mrs$shadow
    colnames(shadow.pair) <- c("MR1","MR2")
    write.table(shadow.pair,file=paste0(SAVE_IN,name,"/msVIPER_ShadowPair.txt"),sep="\t",col.names=NA)

    pdf(paste0(SAVE_IN,name,"/msVIPER_shadow_plot.pdf"))
    plot(mrs, 25, cex = .7)
    dev.off()
    
    #Synergy analysis
    mrs <- msviperCombinatorial(mrs, regulators = 25, verbose = FALSE)
    mrs <- msviperSynergy(mrs, verbose = FALSE)

    pdf(paste0(SAVE_IN,name,"/msVIPER_sinergy_plot.pdf"))
    plot(mrs, 25, cex = .7)
    dev.off()

}






