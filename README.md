# Variant calling pipeline for mixed-ploidy datasets
## Introduction
This pipeline takes raw, paired-end Illumina reads of genomic datasets from plant populations with variable ploidy, aligns them against a (diploid) reference genome using BWA and does variant calling and filtering of results using GATK4. This variant calling workflow is similar to the [Fastq-to-vcf pipeline](https://github.com/vlkofly/Fastq-to-vcf) by Jakub Vlcek, however it implements a number of changes including:
* [GATK4](https://www.broadinstitute.org/news/broad-institute-releases-open-source-gatk4-software-genome-analysis-optimized-speed-and) instead of GATK3
* use of [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow system instead of a collection of Python/Bash scripts
* alternative filtering pipeline using bcftools instead of GATK

In the following part I will describe how to set up, configure and run the workflow on HPC computing clusters, both in general and specifically for users of the czech national computing infrastructure [MetaCentrum](https://metavo.metacentrum.cz/).

## General idea
This pipeline is intended to do variant calling on population genomic datasets with variable ploidy. To this there are four general steps:
1. Indexing of reference genome
2. Trimming of Illumina short reads (can now be deactivated bysetting the option **trim_reads** to 0), alignment and marking of PCR duplicates
3. Variant calling on the individual and population level
4. Filtering of variant calls (either using GATK/PicardTools or Bcftools)

### Snakemake
The workflow is implemented using [Snakemake](https://snakemake.readthedocs.io/en/stable/). You can find an introduction about the underlying concepts on the [Snakemake GitHub page](https://snakemake.github.io/). In general, individual tasks are defined as rules, based on which Snakemake determines how to generate desired output files from available input files and generates and runs required software. Snakemake is under continuous development and you should use the most recent version to run this workflow because slightly older versions may not work. I am currently running it under Snakemake 7.32.4.

## Setting up the workflow

The workflow needs a number of software packages. Depending on your cluster configuration, these may be installed in different ways (either by hand, via conda or by loading preconfigured modules for your cluster) and you can adjust the paths and commands for each of the packages individually.

### Installing Snakemake

Snakemake may be installed via [conda or Pip](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). Currently this workflow runs with Snakemake v.7, v.8 is not supported yet.
#### On MetaCentrum:

I installed Snakemake using a local Conda installation. First I dowloaded a minimal Conda installation (with the mamba package solver already preinstalled) [here](https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh) using
```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
```
I started the Conda installation using

```
sh Mambaforge-Linux-x86_64.sh
```
I followed the installation instructions and installed it in my home directory

To the last question, if conda should be activated by default I said no. You can activate your newly installed [conda environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) using
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
Install Snakemake v.7 using [mamba](https://github.com/mamba-org/mamba)
```
mamba install -c bioconda snakemake=7.32.4
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
With the activated conda environment, GATK4 can be run using the "gatk" wrapper script in the extracted GATK4 folder. From GATK4 v.4.4.0.0 Java v.17 or higher is required. On metaCentrum this can be loaded using
```
module load openjdk
```
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
* R

For estimating mean depth of bam files and for generating genome-wide estimates this workflow uses [PanDepth](https://github.com/HuiyangYu/PanDepth). PanDepth is very fast, but there is neither a Conda package nor a module on MetaCentrum. However, you can download [pre-build binaries](https://github.com/HuiyangYu/PanDepth/releases). After extracting them, you'll have to set the path to the binary by setting the *PANDEPTH* path in your prerun_script (1_mapreads.sh). 

2. Variant calling on the individual level
* GATK4
3. joint genotyping
* Bedtools
* GATK4
* Bcftools
4. Filtering of variant calls using bcftools
* bcftools
* Numpy
* matplotlib
* bedops
* Bedtools

For having a recent version of matplotlib on Metacentrum, it is recommended to create a Conda environment for the filtering step and install the dependencies using

```
mamba install -c bioconda bcftools numpy matplotlib bedops bedtools
```


### Installing the workflow

You can clone the workflow to your local directory using

```
git clone https://github.com/jgerchen/polyploid_variant_calling
```

#### Configuration

General [configuration of Snakemake](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html) is implemented using yaml files. There is one yaml file in the [config directory](config) for each of the four parts of the workflow, in which default parameters for some options are put. Importantly, you have to provide an additional custom yaml file ([example](examples/example.yaml)) for your specific Snakemake run in which you set a number of paths to directories and other config files, for which there are no default settings. You can also add any of the settings from the default config file to this file to override them.
The individual parameters that have to be set are:
##### Input files
* **input_fasta** Your reference genome in fasta format
* **sample_list** tab-separated file containing information about individual samples format is as follows:
 ```
sample_name<TAB>file_name<TAB>adapter_file<TAB>sample_ploidy
```
Here sample_name will determine the naming of other downstream files for each sample (wildcard {sample}) and file name is used to identify individual fastq files containing reads for this sample (file_name must be contained in the name of each fastq file). Files are identified using the [python glob module](https://docs.python.org/3/library/glob.html), e.g. it will also match bash wildcards contained in the name like * or ?. Multiple file_names can be matched by entering a comma-separated list (with no spaces). Furthermore, paired-end reads have to be matched correctly in pairs. Since there may be more complex situations than for sample names, I'm using the [RE module from python](https://docs.python.org/3/library/re.html) to match REs within the results from the glob module. Multiple PE pairs can be defined using the **paired_end_res** option (default: R1/R2 and fwd/rev), more explanations are in the comments of the [default config file](config/1_mapreads.yaml).
The column adapter_file determines the name of the adapter file used for trimmomatic for this sample (located in **adapter_dir**), sample_ploidy determines sample ploidy (2 for diploid, 4 for tetraploid etc.).

* **sub_intervals** defining genomic intervals over which GATK haplotypecaller and population level genotyping will be parallelized in a [scatter-gather](https://gatk.broadinstitute.org/hc/en-us/articles/360035532012-Parallelism-Multithreading-Scatter-Gather) fashion. Each line defines an interval, over which variant calling will run as a single job. The file is tab-seperated in two columns, the first line defines the name for output file for each interval, the second line is either the name of a single contig or [interval](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists) in the fasta file or, if multiple contigs or intervals should be run in a single job, they can be separated with commas (without spaces) and put in the same line (e.g. interval<TAB>contig1,contig2,contig3,...).

##### Pre-existing directories

* **fastq_dir** Directory in which raw illumina data in fastq.gz format are located. The data can be spread across multiple sub-directories and the script will find files with matching file names defined in **sample_list** recursively. Now you can also define multiple directories using a pythonic list, e.g.: fastq_dir: ["/path/to/dir1", "/path/to/dir2"] 
* **adapter_dir** directory in which sequencing adapters (used by trimmomatic) defined in **sample_list** are located
* **cluster_code_dir** scripts that will be run at the beginning of each rule (separate scripts for each step of the pipeline). This is used to load software modules or snakemake environments on metacentrum, functionality can be deactivated by setting load_cluster_code to 0

##### Output directories (will be automatically created by snakemake if non-existent)
* **fasta_dir** Reference genome will be copied here and associated indexes will be located here
* **fastqc_dir** results of fastqc analysis of reads (before and after trimming) will be put here. FastQC can be deactivated by setting the **run_fastqc** option to 0
* **fastq_trimmed_dir** Fastq files after trimming will be copied here
* **bam_dir** bam files will be copied here
* **gvcf_dir** gvcf files of individual samples (created by GATK HaplotypeCaller) will be put here
* **vcf_dir** unfiltered, population-level vcf files will be put here
* **vcf_filtered** filtered vcf files will be put here
* **log_dir** log files (direct output of stdout and/or stderr or rules and shell scripts) will be put here
* **depthmask_dir** location where the [depthmask script](https://github.com/jgerchen/polyploid_popgen/tree/main/depth_mask) will put the deptmask bed and additional output files
* **report_dir** location where additional plots and stats are saved, which can later be used to generate a [Snakemake report](https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html), although I did not implement extensive plotting functions for all of the rules yet
* **sample_read_file** This is a tab-separated file that collects information from samples, libraries, fastq-files, adapter files and ploidy. This file is generated based on your **sample_list** by running the rule **get_sample_reads**, however you may want to check it to see if forward and reverse reads are assigned correctly to libraries and you can edit it by hand if something went wrong. 

#### Pre-run scripts

In addition, there is one script for each of the four parts of the workflow which contain commands, which are run every time before any of the rules in this part of the workflow is run. The idea is that you can load any software modules or conda packages in the way that is appropriate for your computing environment. You can also set a different set of pre-run scripts for your individual Snakemake runs by changing the folder with the cluster_code_dir option.
If your computing environment does not require to run these scripts you can deactivate this functionality by setting the **load_cluster_code** option to 0.

#### On MetaCentrum:

[Here](examples/prerun_scripts) I provide pre-run scripts, which load the conda environments and metacentrum modules for running the scripts. If you want to use these, you'll have to adjust several paths to conda and gatk4 in 2_callvars.sh, 3_genotypeGVCF.sh and 4_filter_GATK.sh (if you want to use GATK for filtering).

## Running the pipeline
In priciple, you have to run Snakemake in the workflow directory (where the main Snakefile is located), giving the desired output file(s) as a parameter. In addition, you'll have to provide your main config file using the --configfile parameter and you'll have to define the number of parallel jobs using the -j parameter. Snakemake will then test if the output can be genrated given the rules and input files. If true, it will run the rules and generate output files. For most downstreaam output files that you're most likely interested in generating, you'll have to set the {species} [wildcard](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards), which will then be automatically used for naming all upstream files.


In general, it is useful to do a dry run, to test if the workflow (with your config file) can be resolved before you actually run it using

```
snakemake output_file -n --configfile config.yaml
```
You can also plot the [DAG](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html#step-4-indexing-read-alignments-and-visualizing-the-dag-of-jobs) of your workflow as a svg without running the actual workflow
```
snakemake output_file --configfile config.yaml --dag | dot -Tsvg > dag.svg          
```
For this you need [graphviz](https://graphviz.org/) installed, on metacentrum you can load it as a module using
```
module load graphviz 
```
While it should be possible to run this workflow on any HPC computing platform, I designed some aspects of it with the architecture of MetaCentrum in mind. That means that currently for each job all input files are copied to a temporary directory (defined by the **temp_dir** option, standard is the $SCRATCHDIR variable used to define the local filesystem of computing notes in MetaCentrum) and the results are copied back to the appropriate directory defined in the config file afterwards while all temporary files, which are not specifically designated as output files of rules, are deleted.
Since the primary merged, unfiltered VCF files (**config[vcf_dir]/{species}.merged.vcf.gz**) can get very large, they are not copied to the temporary directory but directly read by BCFtools in downstream rules. Copying of these files can activated by setting the option **copy_large_vcfs** to 1.

### Running it on MetaCentrum
After you installed the  [Metacentrum cluster profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles), you can run snakemake using 
```
snakemake  output_file --configfile config.yaml --profile snakemake_metacentrum -j 100
```
This assumes that your cluster profile is located in $HOME/.config/snakemake/snakemake_metacentrum and will submit up to 100 jobs in parallel.
For longer workflows, you may want to run Snakemake itself on the [oven node](https://wiki.metacentrum.cz/wiki/Oven_node). You can do this by putting your Snakemake command in a [jobscript](https://wiki.metacentrum.cz/wiki/Beginners_guide#Run_batch_jobs) and submit it using 
```
qsub -q oven jobscript.sh 
```
An example script can be found [here](examples/run_snakemake.sh)
## Output files

You can generate any intermediate files that are defined in output directive of any rule in the workflow. Here **config[option]** are paths defined in the config files and **{wildcard}** are wildcards which are defined either based on the name of ouput files or internally for example based on the names in the sample list.  Likely most interesting output files are:
### 0_index_reference
* **config[fasta_dir]/{species}.fasta** Renamed copy of the reference genome defined by **input\_fasta**
* **config[fasta_dir]/{species}.dict/.fasta.amb/.ann/.bwt/.fai/.pac/.sa** various indexes generated by BWA, samtools and PICARD tools

### 1_mapreads
* **config[fastq_trimmed_dir]/{sample}\_{lib}\_trimmed\_R1/R2/U1/U2.fastq.gz** Reads trimmed using trimmomatic with R1 and R2 for paired forward or reverse reads and U1 and U2 for unpaired forward and reverse reads
* **config[bam_dir]/{species}\_{sample}\_{lib}.bam** Aligned reads of a single PE library {lib} for a single individual {sample}
* **config[bam_dir]/{species}\_{sample}.merged.bam** Merged alignment of all libraries for a single individual {sample} containing all reads
* **config[bam_dir]/{species}\_{sample}.merged.dedup.bam** Merged alignment of all libraries for a single individual {sample} with PCR duplicates marked by PICARD tools MarkDuplicates
* **config[report_dir]/bam_depth/{species}\_{sample}.chr.stat.gz** gzipped table with raw sequencing depth calculated by PanDepth (including duplicates and low quality reads) per contig
* **config[report_dir]/bam_depth/{species}\_{sample}.chr.clean.stat.gz** gzipped table with filtered sequencing depth calculated by PanDepth (duplicates and reads with mapping quality < **pandepth_min_q** removed) per chromosome
* **config[report_dir]/bam_depth/{species}\_{sample}.win.stat.gz** gzipped table with filtered sequencing depth calculated by PanDepth (duplicates and reads with mapping quality < **pandepth_min_q** removed) in sliding windows **config[pandepth_window_size]**
* **config[report_dir]/bam_depth/{species}.bam.table.tsv** tab separated table collecting the results of flagstats and bamDepth for all samples

### 2_callvars
* **config[gvcf_dir]/{species}\_{sample}\_{sub}.gvcf.gz** GVCF file for single individual {sample} generated by GATK HaplotypeCaller for a specific interval {sub} defined in **sub_intervals**
* **config[gvcf_dir]/{species}\_{sub}\_GenomicsDB** Directory containing a [GenomicsDB](https://gatk.broadinstitute.org/hc/en-us/articles/360035891051-GenomicsDB) for all samples for a specific interval {sub}

### 3_genotypeGVCF
* **config[vcf_dir]/{species}\_{sub}.vcf.gz** Unfiltered VCF file for all samples for a specific interval {sub}, containing both variants and invariant sites
* **config[vcf_dir]/{species}.merged.vcf.gz** Unfiltered VCF file for all samples with all intervals merged

### 4_filter
#### GATK (deprecated, probably doesn't work anymore, use BCFtools!)
* **config[vcf_filtered]/{species}.bisel.vcf.gz** VCF file containing all biallelic SNPs (invariants, multi-allelic variants and INDELS and complex variants removed)
* **config[vcf_filtered]/{species}.bifilt.vcf.gz** Biallelic SNPs filtered for GATK best practice values by default (or individual values defined in the config files), with filtered Variants annotated in the FILTER column of the VCF file
* **config[vcf_filtered]/{species}.bipassed.vcf.gz** Biallelic SNPs with SNPs marked as filtered in the previous step removed
* **config[vcf_filtered]/{species}.novarsel.vcf.gz** VCF file containing invariant sites (bi- and multi-allelic variants and INDELS and complex variants removed)
* **config[vcf_filtered]/{species}.novarfilt.vcf.gz** Invariant sites filtered based on the QUAL score (default 15, can be changed by the **invariantQUAL_less** option) annotated in the FILTER column of the VCF file
* **config[vcf_filtered]/{species}.novarpass.vcf.gz** Invariants with sites filtered in the previous step removed
* **config[vcf_filtered]/{species}.merged.vcf.gz** Filtered biallelic SNPs and Invariants merged into a single VCF
* **config[vcf_filtered]/{species}.merged.filtered.vcf.gz** Biallelic SNPs and variants with additional exclusion of regions with excessive depth (defined in {species}\_dm.bed) and optional additonal bed file defining regions with excess heterozygosity (can be set using the **hetmask** option to a different value than the standard value "None")
* **config[vcf_filtered]/{species}.bipassed.dp.vcf.gz** Filtered bi-allic SNPs with genotypes marked that don't pass the minimum depth threshold (default 8, can be set using the **gen_min_depth** option)
* **config[vcf_filtered]/{species}.bipassed.dp.nc.vcf.gz** VCF file with SNPs filtered in the previous step set to missing
* **config[vcf_filtered]/{species}.bipassed.dp.nc.m.vcf.gz** VCF file with sites removed, which have a greater proportion of missing data than a pre-defined threshold (default 0.5, can be set using the **gen_max_missing** option)
* **config[vcf_filtered]/{species}.fourfold.filtered.vcf.gz** biallelic SNPs and Invariants merged into a single VCF, subset to contain only fourfold degenerate sites (defined by the **fourfold** option)
* **config[vcf_filtered]/{species}.bi.fourfold.dp.vcf.gz** Filtered bi-allic SNPs with genotypes marked that don't pass the minimum depth threshold (default 8, can be set using the **gen_min_depth** option), hard filtered for **fourfold** degenrate site and **depthmask** and optional **hetmask** 
* **config[vcf_filtered]/{species}.bi.fourfold.dp.nc.vcf.gz** previous file, but with marked genotypes set to no-call
* **config[vcf_filtered]/{species}.bi.fourfold.dp.nc.m.vcf.gz** previous file, but with sites removed, which have a greater proportion of missing data than a pre-defined threshold (default 0.5, can be set using the **gen_max_missing** option)

#### BCFtools
* **config[vcf_filtered]/{species}.bisel.bt.vcf.gz** VCF file containing all biallelic SNPs (invariants, multi-allelic variants and INDELS and complex variants removed)
* **config[vcf_filtered]/{species}.bipassed.bt.vcf.gz** Biallelic SNPs with SNPs marked as filtered in the previous step removed and sites filtered with [GATK best practices filters](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)
* **config[vcf_filtered]/{species}.novarpass.bt.vcf.gz** Invariants with sites filtered with either QUAL scores or based on sequencing depth at individual genotypes, followed by missing data filtering
* **config[vcf_filtered]/{species}.multivar.bt.vcf.gz** All multiallelic sites, no filtering implemented.
* **config[vcf_filtered]/{species}.bigt.dp.bt.vcf.gz** Biallelic SNPs, with sites with sequencing depth<**gen_min_depth** set to no-call
* **config[vcf_filtered]/{species}.bigt.dp.m.bt.vcf.gz** Biallelic SNPs, with sites removed with >**gen_max_missing** missing data removed
* **config[vcf_filtered]/{species}.merged.bt.vcf.gz** Filtered biallelic SNPs and Invariants merged into a single VCF
* **config[vcf_filtered]/{species}.merged.masked.bt.vcf.gz"** Biallelic SNPs and variants with additional exclusion of regions with excessive depth (defined in {species}\_dm.bed) and optional additonal bed file defining regions with excess heterozygosity (can be set using the **hetmask** option to a different value than the standard value "None")
* **config[vcf_filtered]/{species}.fourfold.bt.vcf.gz"** biallelic SNPs and Invariants merged into a single VCF, subset to contain only fourfold degenerate sites (defined by the **fourfold** option)
* **config[depthmask_dir]/{species}\_dm.bt.bed** Bedfile containing genomic regions with excessive depth, identified by the [depthmask script](https://github.com/jgerchen/polyploid_popgen/tree/main/depth_mask). Proportion of samples with excess depth required for site to be included in depth mask can be set with the **depthmask_prop_excess_depth** parameter

### Limitations of wildcards
I use underscores to separate multiple wildcards in filenames. Since Snakemakes tends to resolve wildcards greedily based on filenames, using underscores in your own wildcards can result in problems and undefined behaviour. This applies to the species wildcard as well as sample and contig names set in you sample and contig list.   

### Setting resources
Each rule has a standard value for [resources](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources). Memory (in mb) is defined by **mem_mb**, disk space (in mb) is defined by **disk_mb** and runtime is now defined in human readable format (1d for one day, 1h for 1hour 1m for one minute etc.) by **runtime**. Default resources are defined in the config files. Alternatively, these can be set manually using the --set-resources RULE:RESOURCE=VALUE option of Snakemake. So for example to set the memory requirement for rule GenotypeGenomicsDBSub to 100000 mb add the option  --set-resources GenotypeGenomicsDBSub:mem_mb=100000 to your snakemake command. Similarly, to set the number of threads for a specific rule you can use --set-threads RULE=THREADS. However, increasing the number of threads may potentially result in runtime improvements only for the rules **trimmomatic** and **map_reads**, but not for variant calling or genotyping or filtering.

### Failed jobs
You can set the number of times Snakemake is supposed to resubmit a failed job by setting the --retries option of Snakemake. If you don't set resources yourself using the --set-resources option (which overrides this functionality) you can set that the standard values for **mem_mb**, **disk_mb** and **runtime** are automatically increased with every rerun. The factor by which each resource is increased with every rerun can be set by the **repeat_mem_mb_factor** (for memory), **repeat_disk_mb_factor** (for disk space) and **repeat_runtime_factor** (for runtime) options in the Snakefile. If you set it to 0, the standard value will remain the same with every attempt, if you set it to 0.5 it will increase linearly by 50% with every attempt etc.

### Issues with GATK
Unfortunately GATK has many unsolved issues, which cause problems in this type of analysis, especially when dealing with invariant sites. Below are a number of relevant and unsolved issues that I noted:
* There are no more QUAL scores at most/all invariant sites. See [this report](https://gatk.broadinstitute.org/hc/en-us/community/posts/4412548615579-Not-quality-in-invariant-sites). This makes it impossible to filter invariant sites based on QUAL scores. Instead I filter sites based on sequencing depth at genotypes and then on missing data. In addition, for some invariant sites the QUAL score may be noted as "Infinity" or similar. However, there does not seem to be any correlation with this score and sequencing depth, so it does not indicate particularly high confidence sites.
* GATK tends to create non-sensical "long" reference alleles for some invariant sites, which can span several other following sites, which are properly called as either invariant or variant. See [this report](https://github.com/broadinstitute/gatk/issues/8030). So far, I leave these sites in the VCF files, if they are an issue with downstream analyses they should be filtered by hand.
* Starting with GATK4.2.3 genotypes without sequencing coverage are now called as homozygous reference (e.g. 0/0) instead of no-call (e.g. ./.) as before. According to GATK [this is not a bug but a feature](https://gatk.broadinstitute.org/hc/en-us/articles/6012243429531-GenotypeGVCFs-and-the-death-of-the-dot), unfortunately it took the GATK team more than two years to tell anybody about it. 


## Still to implement
While the workflow should be functional as it is now there are still several features that I want to implement or improve

* Improved plotting functions either using Snakemake reports or via direct integration of [multiQC](https://multiqc.info/)
* Implementation of [updog](https://dcgerard.github.io/updog/) for improved polyploid genotyping

