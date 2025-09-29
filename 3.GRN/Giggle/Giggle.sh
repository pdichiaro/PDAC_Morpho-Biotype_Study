

#!/bin/sh
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: Giggle.sh
##
## Description: 
#~
## This tool will run GIGGLE to build an index from ChIP-seq data and assess the enrichment of input promoter regions.
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



#
GIGGLE_ROOT=".../3.GRN/Giggle/giggle/giggle-singularity/"
conf_path=".../3.GRN/Giggle/"
OUT_FOLDER="/3.GRN/Giggle/"
#

mkdir $OUT_FOLDER/Results


###########################################################################
# Pull image and Check parameters (first time)

#cd $GIGGLE_ROOT
#$GIGGLE_ROOT/giggle.sh pull -C $GIGGLE_ROOT/sample_data/config.ini
#$GIGGLE_ROOT/giggle.sh check -C $GIGGLE_ROOT/sample_data/config.ini

###########################################################################
#---  Check parameters (every time) ---###
cd $OUT_FOLDER
$GIGGLE_ROOT/giggle.sh check -C $conf_path/config.ini


###---  Sort and index your reference files ---###
if find "$OUT_FOLDER/reference_data/" -type f -name "*.bed" | grep -q .; then
    echo "Found .bed files in ./reference_data/ to sort and index."
    for file in "$OUT_FOLDER/reference_data/"*.bed; do
        bgzip "$file"
    done
else
    echo "No .bed files found in ./reference_data/ to compress."
fi

$GIGGLE_ROOT/giggle.sh index -C $conf_path/config.ini


###--- Enrichment analysis ---###
#effective genome size ~ 2.9Gb (2900000000)
#http://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html

if find "$OUT_FOLDER/Targets/" -type f -name "*.bed" | grep -q .; then
    echo "Found .bed files in ./reference_data/ to sort and index."
    for file in "$OUT_FOLDER/Targets/"*.bed; do
        bgzip "$file"
    done
else
    echo "No .bed files found in ./Targets/ to compress."
fi

for f in $(ls $OUT_FOLDER/Targets | sed 's/.promoters.bed.gz//g');
do
echo "Processing $f"

$GIGGLE_ROOT/giggle.sh search -C $conf_path/config.ini \
    -q $OUT_FOLDER/Targets/"$f"*promoters.bed.gz \
    -g 3095677412 \
    -s > $OUT_FOLDER/Results/"$f"_giggle_output.txt
done


