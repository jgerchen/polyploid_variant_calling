# Variant calling pipeline for mixed-ploidy datasets
## Introduction
This variant calling workflow is similar to the [Fastq-to-vcf pipeline](https://github.com/vlkofly/Fastq-to-vcf) by Jakub Vlcek, however it implements a number of changes including:
* [GATK4](https://www.broadinstitute.org/news/broad-institute-releases-open-source-gatk4-software-genome-analysis-optimized-speed-and)
* use of [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow system

In the following part I will describe how to set up, configure and run the workflow, both in general and specifically for the czech national computing infrastructure [MetaCentrum](https://metavo.metacentrum.cz/).  
