# Variant calling pipeline for mixed-ploidy datasets
## Introduction
This pipeline takes raw, paired-end Illumina reads of genomic datasets from plant populations with variable ploidy, aligns them against a (diploid) reference genome using BWA and does variant calling and filtering of results using GATK4. This variant calling workflow is similar to the [Fastq-to-vcf pipeline](https://github.com/vlkofly/Fastq-to-vcf) by Jakub Vlcek, however it implements a number of changes including:
* [GATK4](https://www.broadinstitute.org/news/broad-institute-releases-open-source-gatk4-software-genome-analysis-optimized-speed-and) instead of GATK3
* use of [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow system instead of a collection of Python/Bash scripts

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

The workflow needs a number of software packages. Depending on your cluster configuration, these may be installed in different ways (either by hand, via conda or by loading preconfigured modules for your cluster) and you can adjust the paths and commands for each of the packages individually.

### Installing Snakemake

Snakemake may be installed via [conda or Pip](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). Make sure to install a recent version of Snakemake, since Snakemake is under conti¨nuous development and some of the features used in this workflow may not be supported by older Snakemake versions.
#### On MetaCentrum:

I installed Snakemake using a local Conda installation. First I dowloaded a minimal Conda installation (with the mamba package solver already preinstalled) [here](https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh). I followed the installation instructions and installed it in my home directory. To the last question, if conda should be activated by default I said no. 

### Installing GATK4

### Installing other software

### Installing the workflow



