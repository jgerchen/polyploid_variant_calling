# Variant calling pipeline for mixed-ploidy datasets
## Introduction
This pipeline takes raw, paired-end Illumina reads of genomic datasets from plant populations with variable ploidy, aligns them against a (diploid) reference genome using BWA and does variant calling and filtering of results using GATK4. This variant calling workflow is similar to the [Fastq-to-vcf pipeline](https://github.com/vlkofly/Fastq-to-vcf) by Jakub Vlcek, however it implements a number of changes including:
* [GATK4](https://www.broadinstitute.org/news/broad-institute-releases-open-source-gatk4-software-genome-analysis-optimized-speed-and) instead of GATK3
* use of [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow system instead of a collection of Python/Bash scripts
* optional alternative filtering pipeline using bcftools instead of GATK

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
##### Input files
* **input_fasta** Your reference genome in fasta format
* **sample_list** tab-separated file containing information about individual samples format is as follows:
 ```
sample_name<TAB>file_name<TAB>adaptor_file
```
Here sample_name will determine the naming of other downstream files for each sample (wildcard {sample}) and file name is used to identify individual fastq files containing reads for this sample (file_name must be contained in the name of each fastq file). Finally, adaptor_file determines the adaptor file used for trimmomatic for this sample.
* **sample_ploidies** tab-separated file containing sample name and ploidy of individual
* **sub_intervals** defining genomic intervals over which GATK haplotypecaller and population level genotyping will be parallelized in a [scatter-gather](https://gatk.broadinstitute.org/hc/en-us/articles/360035532012-Parallelism-Multithreading-Scatter-Gather) fashion. Each line defines an interval, over which variant calling will run as a single job. The file is tab-seperated in two columns, the first line defines the name for output file for each interval, the second line is either the name of the contig in the fasta file or, if multiple contigs should be run in a single job, the location of a second file that contains a list of all the contigs, that should be run in a single job.
* **fourfold** BED file containing the location of fourfold-degenerate sites.

##### Pre-existing directories

* **fastq_dir** Directory in which raw illumina data in fastq.gz format are located. The data can be spread across multiple sub-directories and the script will find files with matching file names defined in **sample_list** recursively
* **adapter_dir** directory in which sequencing adapters (used by trimmomatic) defined in **sample_list** are located
##### Output directories (will be automatically created by snakemake if non-existent)
* **fasta_dir** Reference genome will be copied here and associated indexes will be located here
* **fastqc_dir** results of fastqc analysis of reads (before and after trimming) will be put here. FastQC can be deactivated by setting the **run_fastqc** option to 0
* **fastq_trimmed_dir** Fastq files after trimming will be copied here
* **bam_dir** bam files will be copied here
* **gvcf_dir** gvcf files of individual samples (created by GATK HaplotypeCaller) will be put here
* **vcf_dir** unfiltered, population-level vcf files will be put here
* **vcf_filtered** filtered vcf files will be put here
* **log_dir** log files (direct output of stdout and/or stderr or rules and shell scripts) will be put here
* **cluster_code_dir** scripts that will be run at the beginning of each rule (separate scripts for each step of the pipeline). This is used to load software modules or snakemake environments on metacentrum, functionality can be deactivated by setting load_cluster_code to 0
* **depthmask_dir** location where the [depthmask script](https://github.com/jgerchen/polyploid_popgen/tree/main/depth_mask) will put the deptmask bed and additional output files
* **report_dir** location where additional plots and stats are saved, which can later be used to generate a [Snakemake report](https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html), although I did not implement extensive plotting functions for all of the rules yet

#### Pre-run scripts

In addition, there is one script for each of the four parts of the workflow which contain commands, which are run every time before any of the rules in this part of the workflow is run. The idea is that you can load any software modules or conda packages in the way that is appropriate for your computing environment. You can also set a different set of pre-run scripts for your individual Snakemake runs by changing the folder with the cluster_code_dir option.
If your computing environment does not require to run these scripts you can deactivate this functionality by setting the load_cluster_code option to 0.

#### On MetaCentrum:

Here I provide pre-run scripts, which load the conda environments and metacentrum modules for running the scripts. If you want to use these, you'll have to adjust several paths the paths to conda and gatk4 in 2_callvars.sh, 3_genotypeGVCF.sh and 4_filter_GATK.sh (if you want to use GATK for filtering).

## Running the pipeline
In priciple, you have to run Snakemake in the workflow directory (where the main Snakefile is located), giving the desired output file(s) as a parameter. In addition, you'll have to provide your main config file using the --configfile parameter and you'll have to define the number of parallel jobs using the -j parameter. Snakemake will then test if the output can be genrated given the rules and input files. If true, it will run the rules and generate output files. For most downstreaam output files that you're most likely interested in generating, you'll have to set the {species} [wildcard](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards), which will then be automatically used for naming all upstream files.
In general, it is useful to do a dry run, to test if the workflow (with your config file) can be resolved before you actually run it using

```
snakemake output_file -n --configfile config.yaml
```
You can also plot the [DAG](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html#step-4-indexing-read-alignments-and-visualizing-the-dag-of-jobs) of your workflow as a pdf without running the actual workflow
```
snakemake output_file --configfile config.yaml --dag | dot -Tpdf -Gsize="11,8│!" > dag.pdf           
```
For this you need [graphviz](https://graphviz.org/) installed, on metacentrum you can load it as a module using
```
module load graphviz          
```
While it should be possible to run this workflow on any HPC computing platform, I designed some aspects of it with the architecture of MetaCentrum in mind. That means that currently for each job all input files are copied to a temporary directory (defined by the **temp_dir** option, standard is the $SCRATCHDIR variable used to define the local filesystem of computing notes in MetaCentrum) and the results are copied back to the appropriate directory defined in the config file afterwards while all temporary files, which are not specifically designated as output files of rules, are deleted.

### Running it on MetaCentrum
After you installed the  [Metacentrum cluster profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles), you can run snakemake using 
```
snakemake  output_file --configfile config.yaml --profile snakemake_metacentrum -j 100
```
This assumes that your cluster profile is located in $HOME/.config/snakemake/snakemake_metacentrum and will submit up to 100 jobs in parallel.
For longer workflows, you may want to run Snakemake itself on the [oven node](https://wiki.metacentrum.cz/wiki/Oven_node). You can do this by putting your Snakemake command in a [jobscript](https://wiki.metacentrum.cz/wiki/Beginners_guide#Run_batch_jobs) and submit it using 
```
qsub -q oven jobscipt.sh 
```
## Output files

You can generate any intermediate files that are defined in output directive of any rule in the workflow. Likely the ones that are the most interesting are:
### 0_index_reference

### 1_mapreads

### 2_callvars

### 3_genotypeGVCF

### 4_filter

## Still to implement

While the workflow should be functional as it is now there are still several features that I want to implement or improve

* Improved plotting functions either using Snakemake reports or via direct integration of [multiQC](https://multiqc.info/)
* Better handling of crashed jobs, e.g. automated restart of jobs that run out of memory or walltime with increased resources
* Implementation of [updog](https://dcgerard.github.io/updog/) for improved polyploid genotyping





