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
4. Filtering of variant calls (either using GATK/PicardTools or Bcftools)

### Snakemake
The workflow is implemented using [Snakemake](https://snakemake.readthedocs.io/en/stable/). You can find an introduction about the underlying concepts on the [Snakemake GitHub page](https://snakemake.github.io/). In general, individual tasks are defined as rules, based on which Snakemake determines how to generate desired output files from available input files and generates and runs required software.

## Setting up the workflow

The workflow needs a number of software packages. Depending on your cluster configuration, these may be installed in different ways (either by hand, via conda or by loading preconfigured modules for your cluster) and you can adjust the paths and commands for each of the packages individually.

### Installing Snakemake

Snakemake may be installed via [conda or Pip](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). Make sure to install a recent version of Snakemake, since Snakemake is under continuous development and some of the features used in this workflow may not be supported by older Snakemake versions.
#### On MetaCentrum:

I installed Snakemake using a local Conda installation. First I dowloaded a minimal Conda installation (with the mamba package solver already preinstalled) [here](https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh). I followed the installation instructions and installed it in my home directory. To the last question, if conda should be activated by default I said no. You can activate your newly installed [conda environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) using
```
source yourcondadir/bin/activate
```
Create a new conda environment called "Snakemake" using
```
conda create --name Snakemake
```
This environment is activated using
```
conda activate Snakemake
```
Install Snakemake using [mamba](https://github.com/mamba-org/mamba)
```
mamba install -c bioconda snakemake
```
After this you can also directly activate the Snakemake environment using
```
source yourcondadir/bin/activate Snakemake
```
In addition, Snakemake requires a [cluster profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles), which tells it how to submit and queue jobs. I developed a specific cluster profile for MetaCentrum. The profile and instructions for installation and use can be found [here](https://github.com/jgerchen/snakemake_metacentrum).


### Installing GATK4

There are different ways to install GATK4 via [Docker](https://gatk.broadinstitute.org/hc/en-us/articles/360035889991) or with the help of [Conda](https://gatk.broadinstitute.org/hc/en-us/articles/360035889851--How-to-Install-and-use-Conda-for-GATK4). Importantly, there are extensive Python dependencies and it appears that just installing the standard Conda package is insufficient to get a fully functional version of GATK4, but Conda can be used to help [installing the dependencies in combination with a local installation of GATK4](https://gatk.broadinstitute.org/hc/en-us/articles/360035889991).

#### On MetaCentrum:

On metacentrum, I follow the [instructions for installing GATK4 with the help of Conda](https://gatk.broadinstitute.org/hc/en-us/articles/360035889851--How-to-Install-and-use-Conda-for-GATK4).
First I download the [latest GATK4 package](https://github.com/broadinstitute/gatk/releases), extract it somewhere in my homedirectory and enter the unzipped directory
```
unzip gatk-4.4.0.0.zip
cd gatk-4.4.0.0
```
I activate Conda (see above) and create a new Conda environment, in which Python dependencies are installed
```
source yourcondadir/bin/activate
conda env create -n gatk4 -f gatkcondaenv.yml
```
I can activate the Conda environment using
```
source yourcondadir/bin/activate gatk4
```
It can also be convenient to install the Java runtime environment into this conda environment by using
```
mamba install -c conda-forge openjdk
```
With the activated conda environment, GATK4 can be run using the "gatk" wrapper script in the extracted GATK4 folder

### Installing other software

This workflow uses a bunch of other software packages, which are available as modules on MetaCentrum. Specifically, for each part of the workflow they are:

0. Indexing of reference genome
* BWA
* Samtools
* PicardTools
* htslib
1. Trimming of Illumina short reads, alignment and marking of PCR duplicates
* Trimmomatic
* Samtools
* BWA
* PicardTools
* fastQC (optional)
2. Variant calling on the individual level
* GATK4
3. joint genotyping
* Bedtools
* GATK4
4. Filtering of variant calls
* GATK4 and PicardTools or Bcftools
* Bedtools
### Installing the workflow

You can clone the workflow to your local directory using

```
git clone https://github.com/jgerchen/polyploid_variant_calling
```

#### Configuration

General [configuration of Snakemake](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html) is implemented using yaml files. There is one yaml file in the config directory for each of the four parts of the workflow, in which default parameters for some options are put. Importantly, you have to provide an additional custom yaml file for your specific Snakemake run in which you set a number of paths to directories and other config files, for which there are no default settings. You can also add any of the settings from the default config file to this file to override them.
The individual parameters that have to be set are:



Your custom yaml file has to contain the following entries:

#### Pre-run scripts

In addition, there is one script for each of the four parts of the workflow which contain commands, which are run every time before any of the rules in this part of the workflow is run. The idea is that you can load any software modules or conda packages in the way that is appropriate for your computing environment. You can also set a different set of pre-run scripts for your individual Snakemake runs by changing the folder with the cluster_code_dir option.
If your computing environment does not require to run these scripts you can deactivate this functionality by setting the load_cluster_code option to 0.

#### On MetaCentrum:

Here I provide pre-run scripts, which load the conda environments and metacentrum modules for running the scripts. If you want to use these, you'll have to adjust several paths: 

## Running the pipeline





