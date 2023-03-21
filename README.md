# Variant calling pipeline for mixed-ploidy datasets
## Introduction
This variant calling workflow is similar to the [Fastq-to-vcf pipeline](https://github.com/vlkofly/Fastq-to-vcf) by Jakub Vlcek, however it implements a number of changes including:
* [GATK4](https://www.broadinstitute.org/news/broad-institute-releases-open-source-gatk4-software-genome-analysis-optimized-speed-and)
* use of [Snakemake](https://snakemake.readthedocs.io/en/stable/)
