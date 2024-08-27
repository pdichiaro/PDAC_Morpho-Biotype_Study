#!/bin/sh
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: RNAseq_LMD.sh
##
## Description: 
#~
## Master script to run multiple alignements via qsub()
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



## Insert info for the following commands ##
##########
Cores=10
mem=15G
UsefulData=/1.LMDseq/UsefulData/
OUT_FOLDER=/.../Output/
strandness=reverse
##########


#PATH 
export PATH=/opt/pbs/bin/:$PATH

#Run psudoalignment (KALLISTO)
args=()
for f in $OUT_FOLDER/1_fastq/*_R1.fastq
do
x=$(basename "$f" _R1.fastq)
echo "SEND JOB-02_RNAseq."$x""

VAR=Cores=$Cores,OUT_FOLDER=$OUT_FOLDER,strandness=$strandness,UsefulData=$UsefulData,x=$x
ID=$(qsub -N "$x".02_RNAseq -v $VAR -o $OUT_FOLDER/0_LOG/02_RNAseq."$x".log -e $OUT_FOLDER/0_LOG/02_RNAseq."$x".Error.log -l select=1:ncpus=$Cores:mem=$mem /.../RNAseq_LMD.sh)

args+=($ID) 
done


