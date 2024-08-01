# PDAC_Morpho-Biotype_Study
A repository of code for data analysis/processing generated for the PDAC Morpho-Biotype Study.

The study was performed at the Transcriptional Control in Cancer Lab at European Institute of Oncology (IEO). These data contains transcriptional profiles of multiple tumor areas (200-500 cells) dissected by Laser Micro-dissection (LMD) from primary patients affected by Pancreatic Ductal Adenocarcinoma (PDAC).

Primary findings related to the identification of coexisting Morpho-Biotype in PDAC were published here:
[Di Chiaro et al. Cancer Cell 2024](https://www.cell.com/cancer-cell/fulltext/S1535-6108(24)00079-5)

The GEO SuperSeries can be found here: [GSE209952](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE209952) 

Count matrices, meta data tables and R objects useful for the data analysis are deposited on Zenodo. 
DOI: [10.5281/zenodo.12680172](https://doi.org/10.5281/zenodo.12680172)
![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12680172.svg) 

![moprho_biotypes](https://github.com/user-attachments/assets/f45c1fa8-02a8-4bc7-bab7-07ce1bfb2ac9)

# Directory structure
```
└─ 1. LMDseq/:
|   └─ RNAseq_LMD_master.sh: Master script to run RNAseq_LMD.sh for multiple samples
|   └─ RNAseq_LMD.sh: Script to pseudo-align fastq files generated by laser microdissection (LMD) PDAC samples
|   └─ RNAseq.LMD.txi.import.R: Script for the generationn of a reference database
|   └─ session_R_info.txt: txt file collecting info of R session related to RNAseq.LMD.txi.import.R
|   └─ env_yml/:
|   |   └─ trim.yml/: yml file to set conda envirnoment for RNAseq_LMD.sh (trimming and QC)
|   |   └─ kallisto.yml/: yml file to set conda envirnoment for RNAseq_LMD.sh (pseudo-mapping)
└─ 2. Motif analysis/:
|   └─ promoters.R: Script to extract promoters in fasta format
|   └─ pscan.sh: Script to run pscan for motif analysis using custom PWMs
|   └─ session_R_info.txt: txt file collecting info of R session related to promoters.R
|   └─ env_yml/:
|   |   └─ pscan.yml/: yml file to set conda envirnoment for pscan.sh
└─ 3. GRN/:
|   └─ Aracne.sh: Script to run ARACNE-AP for the reconstruction of GRNs using the gene expression matrix and a list of TFs
|   └─ viper_script.R: Analysis code to infer master regulators usnig VIPER tool
|   └─ session_R_info.txt: txt file collecting info of R session related to viper_script.R
└─ 4. scRNAseq/:
|   └─ scRNAseq_step1.R: Analysis code to create Seurat object and analyse data from multiple PDAC patients 
|   └─ scRNAseq_step2.R: Analysis code to harmonize multiple scRNAseq datasets using Harmony
|   └─ scRNAseq_step3.R: Analysis code to cell cluster and score tumor cells in integrated Seurat object
|   └─ session_R_info.txt: txt file collecting info of R session related to the scripts above
```
