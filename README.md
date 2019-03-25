# dna-seq-pipeline
Python-implemented GATK Best Practices SNPs/Indels calling pipeline

1. Introduction

The SNPs/Indels best practice at GATK website is written in WDL. It is mainly suitable for large-scale parallel process on clusters. But one could feel it difficult to deploy and use on desktop and small workstation.

This pipeline written in python is much easier to run, configured and modified. Currently the preprocess in this pipeline is the same as GATK preprocess. Germline caller is varscan but could be changed conveniently. Somatic caller is GATK mutect2.

2. Dependency

a, common utility, can be download at https://github.com/yaotianran/common/blob/master/common.py
b, python packages: pandas, numpy

(to be continued)

