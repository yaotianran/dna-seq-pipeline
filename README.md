# dna-seq-pipeline
Python-implemented GATK Best Practices SNPs/Indels calling pipeline

1. Introduction

The SNPs/Indels best practice at GATK website is written in WDL. It is mainly suitable for large-scale parallel process on clusters. But one could feel it difficult to deploy and use on desktops and small workstations.

This pipeline written in python is much easier to run, configured and modified. Currently the fastq preprocess and somatic calling process in this pipeline is the same as the pipeline suggested by GATK team (GATK BaseRecalibrator and GATK Mutect2 pipeline). Germline caller is varscan2 but could be changed in the future conveniently.

2. Dependency

a, common utility, can be download at https://github.com/yaotianran/common/blob/master/common.py

b, python packages: pandas, numpy

(to be continued)

