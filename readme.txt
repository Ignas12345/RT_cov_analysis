This repository contains some code and notebooks for analyzing the coverage of different RT_enzymes. Code in this repository can be used to create gene body coverage plots for genes of choice. Also the notebooks can be used for exploring different aspects of the data. Hopefully the use of each notebook is explained in the actual notebook.

Main data files are expected to be stored in such manner:
├── 25_SSCV_KG_01_S1
    ├── 25_SSCV_KG_01_S1_Aligned.sortedByCoord.out.bam #these folders contain the bam files as well as some other stuff 
├── 25_SSCV_KG_02_S2
├── 25_SSCV_KG_03_S3
├── 25_SSCV_KG_04_S4 #data for the maxima enzyme
├── genome_and_annotations
    ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa.modified
    ├── gencode.v41.primary_assembly.annotation.gtf.filtered