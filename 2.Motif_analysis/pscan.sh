#!/bin/sh
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: pscan.sh
##
## Description: 
#~
## This tool will run pscan for motif analysis using custom PWMs
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


source activate pscan

export PATH=.../pscan:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/.../.conda/envs/pscan/lib

CORES=8
MOTIF_path=.../PWMs/
INPUT=.../PROMOTERS/
OUT_FOLDER=.../Output/

mkdir $OUT_FOLDER/MOTIFS

for input in $INPUT/*.fa
do  
    FOLDER=$(basename $input "_promoters_500_50.fa")
    mkdir $OUT_FOLDER/MOTIFS/$FOLDER
    
    pscan -q $input -p $INPUT/BK_DEG_promoters_500_50.fasta -l $MOTIF_path/20180406_pwms_selected.wil
    mv $INPUT/*.res $OUT_FOLDER/MOTIFS/$FOLDER/"$FOLDER"_BK_DEG.res.txt
    
    #pscan -q $input -p $INPUT/BK_promoters_500_50.fasta -l $MOTIF_path/20180406_pwms_selected.wil
    #mv $INPUT/*.res $OUT_FOLDER/MOTIFS/$FOLDER/"$FOLDER"_BK_all.res.txt
done

source deactivate




