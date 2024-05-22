configfile: workflow.source_path("../../config/3_genotypeGVCF.yaml")
from humanfriendly import parse_timespan


def GenotypeGenomicsDBSub_mem_mb(wildcards, attempt):
	return int(config["GenotypeGenomicsDBSub_mem_mb"]+(config["GenotypeGenomicsDBSub_mem_mb"]*(attempt-1)*config["repeat_mem_mb_factor"]))
def GenotypeGenomicsDBSub_disk_mb(wildcards, attempt):
	return int(config["GenotypeGenomicsDBSub_disk_mb"]+(config["GenotypeGenomicsDBSub_disk_mb"]*(attempt-1)*config["repeat_disk_mb_factor"]))
def GenotypeGenomicsDBSub_runtime(wildcards, attempt):
	GenotypeGenomicsDBSub_runtime_seconds=parse_timespan(config["GenotypeGenomicsDBSub_runtime"])
	return str(GenotypeGenomicsDBSub_runtime_seconds+int((GenotypeGenomicsDBSub_runtime_seconds*(attempt-1))*config["repeat_runtime_factor"]))+"s"
	#	GenotypeGenomicsDBSub_runtime_cats=config["GenotypeGenomicsDBSub_runtime"].split(":")
#	return str(int(GenotypeGenomicsDBSub_runtime_cats[0])+int(int(GenotypeGenomicsDBSub_runtime_cats[0])*(attempt-1)*config["repeat_runtime_factor"]))+":"+GenotypeGenomicsDBSub_runtime_cats[1]+":"+GenotypeGenomicsDBSub_runtime_cats[2]
rule GenotypeGenomicsDBSub:
	input:
		genomicsDB=config["gvcf_dir"]+"/{species}_{sub}_GenomicsDB",
		ref_fasta=config["fasta_dir"]+"/{species}.fasta",
		ref_fasta_dict=config["fasta_dir"]+"/{species}.dict",
		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai"
	output:
		out_vcf=config["vcf_dir"]+"/{species}_{sub}.vcf.gz"
	threads: 1
	resources:
		mem_mb=GenotypeGenomicsDBSub_mem_mb,
		disk_mb=GenotypeGenomicsDBSub_disk_mb,
		runtime=GenotypeGenomicsDBSub_runtime
	params:
		sub_interval=lambda wildcards: interval_dict[wildcards.sub]
	log:
		config["log_dir"]+"/GenotypeGenomicsDBSub_{species}_{sub}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/GenotypeGenomicsDBSub_{wildcards.species}_{wildcards.sub}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/3_genotypeGVCF.sh
		fi
		cp -r {input.genomicsDB} $temp_folder
		cp {input.ref_fasta} $temp_folder
		cp {input.ref_fasta_dict} $temp_folder
		cp {input.ref_fasta_fai} $temp_folder
		cd $temp_folder
		mkdir tmp
		sub_interval={params.sub_interval}
		if [[ $sub_interval = *','* ]]
		then
			echo $sub_interval | sed -e $'s/,/\\n/g' > sub_intervals.list
			sub_interval=sub_intervals.list
		fi
		if [ {config[GATK_GenotypeGVCF_ignore_crash]} -eq 1  ]
		then
			! $GATK4 --java-options \"-Xmx{resources[mem_mb]}m\" GenotypeGVCFs  -R {wildcards.species}.fasta -V gendb://{wildcards.species}_{wildcards.sub}_GenomicsDB -L $sub_interval -O {wildcards.species}_{wildcards.sub}.vcf.gz --tmp-dir tmp --include-non-variant-sites  &>> {log}
		else
			$GATK4 --java-options \"-Xmx{resources[mem_mb]}m\" GenotypeGVCFs  -R {wildcards.species}.fasta -V gendb://{wildcards.species}_{wildcards.sub}_GenomicsDB -L $sub_interval -O {wildcards.species}_{wildcards.sub}.vcf.gz --tmp-dir tmp --include-non-variant-sites  &>> {log}
		fi
		cp {wildcards.species}_{wildcards.sub}.vcf.gz {output} 
		"""


def get_subVCF_stats_mem_mb(wildcards, attempt):
	return int(config["subVCF_stats_mem_mb"]+(config["subVCF_stats_mem_mb"]*(attempt-1)*config["repeat_mem_mb_factor"]))
def get_subVCF_stats_disk_mb(wildcards, attempt):
	return int(config["subVCF_stats_disk_mb"]+(config["subVCF_stats_disk_mb"]*(attempt-1)*config["repeat_disk_mb_factor"]))
def get_subVCF_stats_runtime(wildcards, attempt):
	subVCF_stats_runtime_seconds=parse_timespan(config["subVCF_stats_runtime"])
	return str(subVCF_stats_runtime_seconds+int((subVCF_stats_runtime_seconds*(attempt-1))*config["repeat_runtime_factor"]))+"s"

rule get_subVCF_stats:
	input:
		out_vcf=config["vcf_dir"]+"/{species}_{sub}.vcf.gz"
	output:
		vcf_stats_pickle=config["report_dir"]+"/MergeSubVCFs/{species}_{sub}.vcfstats.pickle"
	threads: 2
	resources:
		mem_mb=get_subVCF_stats_mem_mb,
		disk_mb=get_subVCF_stats_disk_mb,
		runtime=get_subVCF_stats_runtime
	params:
		bcftools_pickle_script=workflow.source_path("../scripts/parse_bcftools_stdout_pickle.py")
	log:
		config["log_dir"]+"/get_subVCF_stats_{species}_{sub}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/subVCF_stats_{wildcards.species}_{wildcards.sub}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_bcftools.sh
		fi
		cp {input} $temp_folder
		cp {params.bcftools_pickle_script} $temp_folder
		cd $temp_folder
		vcf_n_sites=$(bcftools view --threads 1 {wildcards.species}_{wildcards.sub}.vcf.gz | grep -c ^[^#])
		bcftools view --threads 1 {wildcards.species}_{wildcards.sub}.vcf.gz | python3 parse_bcftools_stdout_pickle.py --n_sites $vcf_n_sites --output {wildcards.species}_{wildcards.sub} --invariants --biallelic --multiallelic
		cp {wildcards.species}_{wildcards.sub}.pickle {output.vcf_stats_pickle} 
		"""

def MergeSubVCFsbcftools_mem_mb(wildcards, attempt):
	return int(config["MergeSubVCFsbcftools_mem_mb"]+(config["MergeSubVCFsbcftools_mem_mb"]*(attempt-1)*config["repeat_mem_mb_factor"]))
def MergeSubVCFsbcftools_disk_mb(wildcards, attempt):
	return int(config["MergeSubVCFsbcftools_disk_mb"]+(config["MergeSubVCFsbcftools_disk_mb"]*(attempt-1)*config["repeat_disk_mb_factor"]))
def MergeSubVCFsbcftools_runtime(wildcards, attempt):
	MergeSubVCFsbcftools_runtime_seconds=parse_timespan(config["MergeSubVCFsbcftools_runtime"])
	return str(MergeSubVCFsbcftools_runtime_seconds+int((MergeSubVCFsbcftools_runtime_seconds*(attempt-1))*config["repeat_runtime_factor"]))+"s"
#MergeSubVCFsbcftools_runtime_cats=config["MergeSubVCFsbcftools_runtime"].split(":")
#	return str(int(MergeSubVCFsbcftools_runtime_cats[0])+int(int(MergeSubVCFsbcftools_runtime_cats[0])*(attempt-1)*config["repeat_runtime_factor"]))+":"+MergeSubVCFsbcftools_runtime_cats[1]+":"+MergeSubVCFsbcftools_runtime_cats[2]
rule MergeSubVCFsbcftools:
	input:
		out_vcf=expand(config["vcf_dir"]+"/{{species}}_{sub}.vcf.gz", sub=interval_list),
		vcf_stats_pickles=expand(config["report_dir"]+"/MergeSubVCFs/{{species}}_{sub}.vcfstats.pickle", sub=interval_list),
		sub_interval_list=config["sub_intervals"]
	output:
		vcf_out=config["vcf_dir"]+"/{species}.merged.vcf.gz",
		vcf_out_index=config["vcf_dir"]+"/{species}.merged.vcf.gz.tbi",
		vcf_stats_table=config["report_dir"]+"/MergeSubVCFs/{species}.merged.tsv",
		vcf_stats_QUAL_biallelic=report(config["report_dir"]+"/MergeSubVCFs/{species}_biallelic_QUAL.pdf", category="MergeSubVCFs", subcategory="general", labels={"variant type":"biallelic", "statistic":"QUAL"}),
		vcf_stats_QUAL_categories_biallelic=report(config["report_dir"]+"/MergeSubVCFs/{species}_biallelic_QUAL_categories.pdf", category="MergeSubVCFs", subcategory="general", labels={"variant type":"biallelic", "statistic":"QUAL_categories"}),
		vcf_stats_INFO_biallelic=report(config["report_dir"]+"/MergeSubVCFs/{species}_biallelic_INFO.pdf", category="MergeSubVCFs", subcategory="INFO", labels={"variant type":"biallelic", "statistic":"INFO"}),
		vcf_stats_GT_counts_biallelic=report(config["report_dir"]+"/MergeSubVCFs/{species}_GT_counts_biallelic.pdf", category="MergeSubVCFs", subcategory="Genotype counts", labels={"variant type":"biallelic", "statistic":"GT counts"}),
		vcf_stats_GT_DP_biallelic=report(config["report_dir"]+"/MergeSubVCFs/{species}_GT_DP_biallelic.pdf", category="MergeSubVCFs", subcategory="Genotype stats", labels={"variant type":"biallelic", "statistic":"GT DP"}),
		vcf_stats_GT_GQ_biallelic=report(config["report_dir"]+"/MergeSubVCFs/{species}_GT_GQ_biallelic.pdf", category="MergeSubVCFs", subcategory="Genotype stats", labels={"variant type":"biallelic", "statistic":"GT GQ"}),
		vcf_stats_n_alt_alleles_multiallelic=report(config["report_dir"]+"/MergeSubVCFs/{species}_multiallelic_n_alt_alleles.pdf", category="MergeSubVCFs", subcategory="general", labels={"variant type":"multiallelic", "statistic":"N alt alleles"}),
		vcf_stats_QUAL_multiallelic=report(config["report_dir"]+"/MergeSubVCFs/{species}_multiallelic_QUAL.pdf", category="MergeSubVCFs", subcategory="general", labels={"variant type":"multiallelic", "statistic":"QUAL"}),
		vcf_stats_QUAL_categories_multiallelic=report(config["report_dir"]+"/MergeSubVCFs/{species}_multiallelic_QUAL_categories.pdf", category="MergeSubVCFs", subcategory="general", labels={"variant type":"multiallelic", "statistic":"QUAL_categories"}),
		vcf_stats_INFO_multiallelic=report(config["report_dir"]+"/MergeSubVCFs/{species}_multiallelic_INFO.pdf", category="MergeSubVCFs", subcategory="INFO", labels={"variant type":"multiallelic", "statistic":"INFO"}),
#		vcf_stats_GT_counts_multiallelic=report(config["report_dir"]+"/MergeSubVCFs/{species}_GT_counts_multiallelic.pdf", category="MergeSubVCFs", subcategory="Genotype counts", labels={"variant type":"multiallelic", "statistic":"GT counts"}),
		vcf_stats_GT_DP_multiallelic=report(config["report_dir"]+"/MergeSubVCFs/{species}_GT_DP_multiallelic.pdf", category="MergeSubVCFs", subcategory="Genotype stats", labels={"variant type":"multiallelic", "statistic":"GT DP"}),
		vcf_stats_GT_GQ_multiallelic=report(config["report_dir"]+"/MergeSubVCFs/{species}_GT_GQ_multiallelic.pdf", category="MergeSubVCFs", subcategory="Genotype stats", labels={"variant type":"multiallelic", "statistic":"GT GQ"}),
		vcf_stats_QUAL_invariant=report(config["report_dir"]+"/MergeSubVCFs/{species}_invariant_QUAL.pdf", category="MergeSubVCFs", subcategory="general", labels={"variant type":"invariant", "statistic":"QUAL"}),
		vcf_stats_QUAL_categories_invariant=report(config["report_dir"]+"/MergeSubVCFs/{species}_invariant_QUAL_categories.pdf", category="MergeSubVCFs", subcategory="general", labels={"variant type":"invariant", "statistic":"QUAL_categories"}),
		vcf_stats_INFO_invariant=report(config["report_dir"]+"/MergeSubVCFs/{species}_invariant_INFO.pdf", category="MergeSubVCFs", subcategory="INFO", labels={"variant type":"invariant", "statistic":"INFO"}),
		vcf_stats_GT_counts_invariant=report(config["report_dir"]+"/MergeSubVCFs/{species}_GT_counts_invariant.pdf", category="MergeSubVCFs", subcategory="Genotype counts", labels={"variant type":"invariant", "statistic":"GT counts"}),
		vcf_stats_GT_DP_invariant=report(config["report_dir"]+"/MergeSubVCFs/{species}_GT_DP_invariant.pdf", category="MergeSubVCFs", subcategory="Genotype stats", labels={"variant type":"invariant", "statistic":"GT DP"}),
		vcf_stats_GT_RGQ_invariant=report(config["report_dir"]+"/MergeSubVCFs/{species}_GT_RGQ_invariant.pdf", category="MergeSubVCFs", subcategory="Genotype stats", labels={"variant type":"invariant", "statistic":"GT RGQ"})
	threads: 3
	resources:
		mem_mb=MergeSubVCFsbcftools_mem_mb,
		disk_mb=MergeSubVCFsbcftools_disk_mb,
		runtime=MergeSubVCFsbcftools_runtime
	params:
		combine_stat_pickles=workflow.source_path("../scripts/combine_parse_bcftools_pickles.py")
	log:
		config["log_dir"]+"/MergeSubVCFsbcftools_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/MergeSubVCFs_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_bcftools.sh
		fi
		cp {input} $temp_folder
		cp {params.combine_stat_pickles} $temp_folder
		cd $temp_folder
		sub_intervals=$(awk -F/ '{{print $NF}}' <<< {input.sub_interval_list})
		#TODO: does this actually work??? No it doesn't! Fix it...
		awk '{{print "{wildcards.species}_"$1".vcf.gz"}}' $sub_intervals > input_files.list
		awk '{{print "{wildcards.species}_"$1".vcfstats.pickle"}}' $sub_intervals > input_pickles.list
		bcftools concat -f input_files.list -n -o {wildcards.species}.merged.bt.vcf.gz
		cp {wildcards.species}.merged.bt.vcf.gz {output.vcf_out}
		tabix {wildcards.species}.merged.bt.vcf.gz 

		python3 combine_parse_bcftools_pickles.py --histogram_bins 50 --output {wildcards.species} --biallelic --invariants --multiallelic --pickle_list input_pickles.list
		cp {wildcards.species}.merged.bt.vcf.gz.tbi {output.vcf_out_index}
		cp {wildcards.species}_table.tsv {output.vcf_stats_table}
		cp {wildcards.species}_QUAL_biallelic.pdf {output.vcf_stats_QUAL_biallelic}
		cp {wildcards.species}_QUAL_categories_biallelic.pdf {output.vcf_stats_QUAL_categories_biallelic}
		cp {wildcards.species}_INFO_biallelic.pdf {output.vcf_stats_INFO_biallelic}
		cp {wildcards.species}_GT_counts_biallelic.pdf {output.vcf_stats_GT_counts_biallelic}
		cp {wildcards.species}_GT_DP_biallelic.pdf {output.vcf_stats_GT_DP_biallelic}
		cp {wildcards.species}_GT_GQ_biallelic.pdf {output.vcf_stats_GT_GQ_biallelic}
		#multiallelic
		cp {wildcards.species}_n_alt_alleles_multiallelic.pdf {output.vcf_stats_n_alt_alleles_multiallelic}	
		cp {wildcards.species}_QUAL_multiallelic.pdf {output.vcf_stats_QUAL_multiallelic}
		cp {wildcards.species}_QUAL_categories_multiallelic.pdf {output.vcf_stats_QUAL_categories_multiallelic}
		cp {wildcards.species}_INFO_multiallelic.pdf {output.vcf_stats_INFO_multiallelic}
		cp {wildcards.species}_GT_DP_multiallelic.pdf {output.vcf_stats_GT_DP_multiallelic}
		cp {wildcards.species}_GT_GQ_multiallelic.pdf {output.vcf_stats_GT_GQ_multiallelic}
		#invariant
		cp {wildcards.species}_QUAL_invariant.pdf {output.vcf_stats_QUAL_invariant}
		cp {wildcards.species}_QUAL_categories_invariant.pdf {output.vcf_stats_QUAL_categories_invariant}
		cp {wildcards.species}_INFO_invariant.pdf {output.vcf_stats_INFO_invariant}
		cp {wildcards.species}_GT_DP_invariant.pdf {output.vcf_stats_GT_DP_invariant}
		cp {wildcards.species}_GT_RGQ_invariant.pdf {output.vcf_stats_GT_RGQ_invariant}
		cp {wildcards.species}_GT_counts_invariant.pdf {output.vcf_stats_GT_counts_invariant}
		"""
