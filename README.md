# dna-seq-pipeline
Python-implemented GATK Best Practices SNPs/Indels calling pipeline

Update 20190403. add tumor mutation burden function. 

Reference: 

[1] Pathogenic variant burden in the ExAC database: an empirical approach to evaluating population data for clinical variant interpretation (2017 Y. Kobayashi et al)

[2] Analysis of 100,000 human cancer genomes reveals the landscape of tumor mutational burden (2017 Z. R. Chalmers et al)

[3] Whole exome sequencing for determination of tumor mutation load in liquid biopsy from advanced cancer patients (2017 F. Koeppel et al)

[4] Methods of measurement for tumor mutational burden in tumor tissue (2018 B. Meléndez et al)

==============

1. Introduction

The SNPs/Indels best practice at GATK website is written in WDL. It is mainly suitable for large-scale parallel process on clusters. But one could feel it difficult to deploy and use on desktops and small workstations.

This pipeline written in python is much easier to run, configured and modified. Currently the fastq preprocess and somatic calling process in this pipeline is the same as the "GATK Best Practices" suggested by GATK team (GATK BaseRecalibrator and GATK Mutect2 pipeline). Germline caller is varscan2 but could be changed in the future conveniently.

This script also uses multiple processes to make the pipeline run faster.

2. Dependency

a, common utility, can be download at https://github.com/yaotianran/common/raw/master/common.py

b, python packages: pandas, numpy

c, samtools, varscan, GATK, picard, annovar etc. (please check the [TOOLS] section in config file at https://github.com/yaotianran/dna-seq-pipeline/blob/master/dna-seq-pipeline.ini)

3. Usage

Simply run “python3 dna-seq-pipeline.py dna-seq-pipeline.ini” 
