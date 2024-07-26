#!/bin/sh
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: Aracne.sh
##
## Description: 
#~
## This tool will run ARACNE-AP for for the reconstruction of GRNs using the gene expression matrix (LMD-seq)
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
## ARACNe-AP requires JDK > 1.8
## custom list of TFs
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#PATH
export _JAVA_OPTIONS="-Xmx15G"

Aracne_PATH=.../ARACNe-AP-master/dist/
OUT_FOLDER=.../Output/
MATRIX=.../Matrix_rLog_Subgroups.txt 
TF=.../TF.txt



mkdir $OUT_FOLDER/ARACNE_output/

#Calculate a threshold for Mutual Information
java -jar $Aracne_PATH/aracne.jar -e $MATRIX -o $OUT_FOLDER/ARACNE_output/ --tfs $TF --pvalue 1E-8 --seed 1 --calculateThreshold

#Run ARACNe on bootstraps of the input matrix
for i in {1..200}
do
    java -jar $Aracne_PATH/aracne.jar -e $MATRIX -o $OUT_FOLDER/ARACNE_output/ --tfs $TF --pvalue 1E-8 --seed $i
done

#Consolidate, i.e. combine the bootstraps into a final network file
java -jar $Aracne_PATH/aracne.jar -o $OUT_FOLDER/ARACNE_output/ --consolidate   #--nobonferroni

awk 'NR>1 {print $0}' $OUT_FOLDER/ARACNE_output/network.txt > $OUT_FOLDER/ARACNE_output/network_viper.txt
rm $OUT_FOLDER/ARACNE_output/bootstrapNetwork_*


