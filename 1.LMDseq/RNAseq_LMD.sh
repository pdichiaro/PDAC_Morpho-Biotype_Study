#!/bin/sh
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: RNAseq_LMD.sh
##
## Description: 
#~
## This tool will map fastq files generated from LMD samples with Kallisto.
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
## Env trim containing: fastqc, cutadapt, trim_galore
## Env kallisto containing: kallisto, bustools, samtools, bedtools, deeptools, ucsc-wigtobigwig, bc
##
## the reference transcriptome in GTF format is downloaded from Gencode(gencode.v33.basic)
## Kallisto index is created by the cDNA fasta (gencode.v33.transcripts)
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


while getopts "h?C:f:o:s:u:" option; do
    case "$option" in
     h|\?)
          echo Usage: $(basename $0) [-C] [-o] [-s] [-u] [-d]
          echo "options:"
          echo "-h, --help                            show brief help"
          echo "-C, --Cores {1-20}                    specify the number of cores"
          echo "-o, --output-dir                      specify a directory to store output"
          echo "-s, --Strandness=library-strandness   specify if the library has to be treated in a strand specific manner: forward, reverse or NO"
          echo "-u, --usefuldata-DB                   specify the DB directory in UsefulData direcory where files(genome, GTF...) are stored"

          exit 0
          ;;
      C)
          Cores=$OPTARG
          if [ $Cores -le 0 -o $Cores -gt 20 ]; then
              echo "-C: Unknown value"
              exit 1
          fi
          ;;
      o)
          OUT_FOLDER=$OPTARG
          if [ ! -d $OUT_FOLDER ]; then
              echo "-o: Output directory does not exist"
              exit 1
          fi
          ;;
      s)
          strandness=$OPTARG
          if [ $strandness != "forward" -a $strandness != "reverse" -a $strandness != "NO" ]; then
              echo "-s: Unknown value"
              exit 1
          fi
          ;;
      u)
          UsefulData=$OPTARG
          if [ ! -d $UsefulData ]; then
              echo "-u: UsefulData directory does not exist"
              exit 1
          fi
          ;;
     esac
done

if [ $Cores != $Cores -o $OUT_FOLDER != $OUT_FOLDER -o $strandness != $strandness -o $UsefulData != $UsefulData ]; then
    echo 'ERROR: [-C] [-o] [-s] [-u] are required'
fi



#PATH
KALLISTO_path=$UsefulData/.../KALLISTOindex/
gencode_path=$UsefulData/.../hg38_gencode/
Chr_size=$UsefulData/.../Chr_size_info/

#Create folders
mkdir $OUT_FOLDER/2_QC/
mkdir $OUT_FOLDER/3_Trim/
mkdir $OUT_FOLDER/4_Bam/
mkdir $OUT_FOLDER/5_Bw/


#Activate conda env
source activate trim

#QC
echo "QC $x"

mkdir $OUT_FOLDER/2_QC/$x
fastqc $OUT_FOLDER/1_fastq/"$x"_R1.fastq -o $OUT_FOLDER/2_QC/$x

#Adapter Trimming
echo "Trimming $x"

mkdir $OUT_FOLDER/3_Trim/$x

trim_galore -j $Cores --illumina --clip_R2 3 --paired $OUT_FOLDER/1_fastq/"$x"_R1.fastq $OUT_FOLDER/1_fastq/"$x"_R2.fastq -q 20 --stringency 3 --length 20 -o $OUT_FOLDER/3_Trim/$x/ 2> $OUT_FOLDER/3_Trim/$x/"$x".cutAdapter.rep.txt
fastqc $OUT_FOLDER/3_Trim/$x/"$x"_R1_val_1.fq -o $OUT_FOLDER/2_QC/$x
fastqc $OUT_FOLDER/3_Trim/$x/"$x"_R2_val_2.fq -o $OUT_FOLDER/2_QC/$x

#Activate conda env
source activate kallisto

mkdir $OUT_FOLDER/4_Bam/$x

#PseudoAlignment and Quantification in TPM and estimated counts 
##In Kallisto, you can pseudoalign transcripts and then quantify the reads 
if [ "$strandness" = "NO" ]; then
    echo "No strandness"
    echo "Aligning and Quantifying $x"
    kallisto quant -t $Cores -i $KALLISTO_path/gencode.v33.transcripts_hg38_k31.idx --gtf $gencode_path/gencode.v33.basic.annotation.gtf -o $OUT_FOLDER/4_Bam/$x -b 50 --single-overhang --genomebam --chromosomes $Chr_size/hg38.chrom.size.txt $OUT_FOLDER/3_Trim/"$x"/"$x"_R1_val_1.fq $OUT_FOLDER/3_Trim/"$x"/"$x"_R2_val_2.fq

    echo "Creating bigWig $x"
    Scale_to=10000000
    factor=$(samtools view -c $OUT_FOLDER/4_Bam/$x/pseudoalignments.bam)
    factor_pe=$(echo "$factor / 2" | bc)
    scale_f=$(echo "scale=7 ; $Scale_to / $factor_pe" | bc)
    bamCoverage -p $Cores -bs 1 --scaleFactor $scale_f -b $OUT_FOLDER/4_Bam/$x/pseudoalignments.bam -o $OUT_FOLDER/5_Bw/"$x".bw

elif [ "$strandness" = "forward" ]; then
    echo "RNA strand Forward"
    echo "Aligning and Quantifying $x"
    kallisto quant -t $Cores -i $KALLISTO_path/gencode.v33.transcripts_hg38_k31.idx --gtf $gencode_path/gencode.v33.basic.annotation.gtf -o $OUT_FOLDER/4_Bam/$x -b 50 --single-overhang --fr-stranded --genomebam --chromosomes $Chr_size/hg38.chrom.size.txt $OUT_FOLDER/3_Trim/"$x"/"$x"_R1_val_1.fq $OUT_FOLDER/3_Trim/"$x"/"$x"_R2_val_2.fq

    echo "Creating bigWig $x"
    Scale_to=10000000
    factor=$(samtools view -c $OUT_FOLDER/4_Bam/$x/pseudoalignments.bam)
    factor_pe=$(echo "$factor / 2" | bc)
    scale_f=$(echo "scale=7 ; $Scale_to / $factor_pe" | bc)
    bamCoverage -p $Cores -bs 1 --scaleFactor $scale_f -b $OUT_FOLDER/4_Bam/$x/pseudoalignments.bam -o $OUT_FOLDER/5_Bw/"$x".bw
    bamCoverage -p $Cores -bs 1 --scaleFactor $scale_f --filterRNAstrand forward -b $OUT_FOLDER/4_Bam/$x/pseudoalignments.bam -o $OUT_FOLDER/5_Bw/"$x".reverse.bw
    bamCoverage -p $Cores -bs 1 --scaleFactor $scale_f --filterRNAstrand reverse -b $OUT_FOLDER/4_Bam/$x/pseudoalignments.bam -o $OUT_FOLDER/5_Bw/"$x".forward.bw

elif [ "$strandness" = "reverse" ]; then
    echo "RNA strand Reverse"
    echo "Aligning and Quantifying $x"
    kallisto quant -t $Cores -i $KALLISTO_path/gencode.v33.transcripts_hg38_k31.idx --gtf $gencode_path/gencode.v33.basic.annotation.gtf -o $OUT_FOLDER/4_Bam/$x -b 50 --single-overhang --rf-stranded --genomebam --chromosomes $Chr_size/hg38.chrom.size.txt $OUT_FOLDER/3_Trim/"$x"/"$x"_R1_val_1.fq $OUT_FOLDER/3_Trim/"$x"/"$x"_R2_val_2.fq

    echo "Creating bigWig $x"
    Scale_to=10000000
    factor=$(samtools view -c $OUT_FOLDER/4_Bam/$x/pseudoalignments.bam)
    factor_pe=$(echo "$factor / 2" | bc)
    scale_f=$(echo "scale=7 ; $Scale_to / $factor_pe" | bc)
    bamCoverage -p $Cores -bs 1 --scaleFactor $scale_f -b $OUT_FOLDER/4_Bam/$x/pseudoalignments.bam -o $OUT_FOLDER/5_Bw/"$x".bw
    bamCoverage -p $Cores -bs 1 --scaleFactor $scale_f --filterRNAstrand forward -b $OUT_FOLDER/4_Bam/$x/pseudoalignments.bam -o $OUT_FOLDER/5_Bw/"$x".reverse.bw
    bamCoverage -p $Cores -bs 1 --scaleFactor $scale_f --filterRNAstrand reverse -b $OUT_FOLDER/4_Bam/$x/pseudoalignments.bam -o $OUT_FOLDER/5_Bw/"$x".forward.bw
fi

rm $OUT_FOLDER/3_Trim/"$x"/"$x"_R1_val_1.fq $OUT_FOLDER/3_Trim/"$x"/"$x"_R2_val_2.fq

source deactivate 




