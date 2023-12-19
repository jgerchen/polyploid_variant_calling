#!/bin/bash
#PBS -N Snakemake
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=10gb
#PBS -l walltime=72:00:00 
#PBS -m n

### Change this to the path to the conda directory in your home directory and make sure that the name of your snakemake conda environment is actually Snakemake
source /storage/brno2/home/gerchenj/mambaforge/bin/activate Snakemake

### Change this to the folder where your workflow is located
cd /storage/brno2/home/gerchenj/polyploid_variant_calling/workflow

### This is the main snakemake command:
### It will generate two final output files alnus_AG001c.merged.dedup.bam and alnus_AG002c.merged.dedup.bam . You will have to adjust it to your own dataset.
### The first part of the output file must be the same as the bam_dir option in your config file, in this case /storage/brno2/home/gerchenj/snakemake_course/bam/
### the second part will be {species}_{sample}.merged.dedup.bam
### {species} is here alnus
### {sample}  is AG001c for the first file and AG002c for the second file, which are the names of specific samples, which have to be in the first column of your sample list file
### You will also have to change the paths to --configfile and --profile to point to the pathsin your own home directory
### You can leave the options -j (number of parallel jobs submitted at the same time) --max-status-checks-per-seconds (frequency at which Snakemake queries the job status) as they are 

#snakemake  /storage/brno2/home/gerchenj/snakemake_course/bam/alnus_AG00{1..3}c.merged.dedup.bam --configfile /storage/brno2/home/gerchenj/snakemake_course/alnus_dtol.yaml --profile snakemake_metacentrum -j 100 --max-status-checks-per-second 1 


#This is the snakemake command for variant calling. It will call variants for all individuals in the samples list. The option --rerun-triggers mtime ensures that Snakemake will only generate output files that do not exist yet, since it likes to rerun Rules where the code has been modified.

snakemake  /storage/brno2/home/gerchenj/snakemake_course/vcf/alnus.merged.vcf.gz --configfile /storage/brno2/home/gerchenj/snakemake_course/alnus_dtol.yaml --profile snakemake_metacentrum -j 100 --max-status-checks-per-second 1 --rerun-triggers mtime


# clean the SCRATCH directory
clean_scratch
