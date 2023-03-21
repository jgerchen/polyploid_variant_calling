# Variant calling pipeline for mixed-ploidy datasets
## Introduction
This pipeline takes raw, paired-end Illumina reads of genomic datasets from plant populations with variable ploidy, aligns them against a (diploid) reference genome using BWA and does variant calling and filtering of results using GATK4. This variant calling workflow is similar to the [Fastq-to-vcf pipeline](https://github.com/vlkofly/Fastq-to-vcf) by Jakub Vlcek, however it implements a number of changes including:
* [GATK4](https://www.broadinstitute.org/news/broad-institute-releases-open-source-gatk4-software-genome-analysis-optimized-speed-and)
* use of [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow system

In the following part I will describe how to set up, configure and run the workflow on HPC computing clusters, both in general and specifically for users of the czech national computing infrastructure [MetaCentrum](https://metavo.metacentrum.cz/).

## General idea
This pipeline is intended to do variant calling on population genomic datasets with variable ploidy. To this there are four general steps:
1. Indexing of reference genome
2. Trimming of Illumina short reads, alignment and marking of PCR duplicates
3. Variant calling on the individual and population level
4. Filtering of variant calls

### Snakemake
The workflow is implemented using [Snakemake](https://snakemake.readthedocs.io/en/stable/). You can find an introduction about the underlying concepts on the [Snakemake GitHub page](https://snakemake.github.io/). In general, individual tasks are defined as rules, based on which Snakemake determines how to generate desired output files from available input files and generates and runs required software.

## Setting up the workflow

### Installing Snakemake

### Installing GATK4

### Installing other software

### Installing the workflow



