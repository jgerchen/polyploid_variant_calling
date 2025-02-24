configfile: workflow.source_path("../../config/4_filter.yaml")
from humanfriendly import parse_timespan

def filter_bcftools_mem_mb(wildcards, attempt):
	return int(config["filter_mem_mb"]+(config["filter_mem_mb"]*(attempt-1)*config["repeat_mem_mb_factor"]))
def filter_bcftools_disk_mb(wildcards, attempt):
	return int(config["filter_disk_mb"]+(config["filter_disk_mb"]*(attempt-1)*config["repeat_disk_mb_factor"]))
def filter_bcftools_runtime(wildcards, attempt):
	bcftools_runtime_seconds=parse_timespan(config["filter_runtime"])
	return str(bcftools_runtime_seconds+int((bcftools_runtime_seconds*(attempt-1))*config["repeat_runtime_factor"]))+"s"
	#	filter_bcftools_runtime_cats=config["filter_runtime"].split(":")
	#return str(int(filter_bcftools_runtime_cats[0])+int(int(filter_bcftools_runtime_cats[0])*(attempt-1)*config["repeat_runtime_factor"]))+":"+filter_bcftools_runtime_cats[1]+":"+filter_bcftools_runtime_cats[2]

rule filter_bcftools_bisnp:
	input:
		vcf_input=config["vcf_dir"]+"/{species}.merged.vcf.gz",
		vcf_stats_table=config["report_dir"]+"/MergeSubVCFs/{species}.merged.tsv"
	output:
		vcf_filtered_mv_biallelic=config["vcf_filtered"]+"/{species}.mv_biallelic.bt.vcf.gz" if config["combine_filtered_mv_biallelic"]==True else [],
		vcf_filtered_mv_biallelic_index=config["vcf_filtered"]+"/{species}.mv_biallelic.bt.vcf.gz.tbi" if config["combine_filtered_mv_biallelic"]==True else [],
		#vcf_stats_table_ma_biallelic=config["report_dir"]+"/filtermultivariant/{species}.ma_biallelic.tsv" if config["combine_filtered_mv_biallelic"]==True else [],
		bisnp_sel=config["vcf_filtered"]+"/{species}.bisel.bt.vcf.gz",
		bisnp_sel_index=config["vcf_filtered"]+"/{species}.bisel.bt.vcf.gz.tbi",
		bisnp_filter=config["vcf_filtered"]+"/{species}.bifilter.bt.vcf.gz",
		bisnp_filter_index=config["vcf_filtered"]+"/{species}.bifilter.bt.vcf.gz.tbi",
		bisnp_passed=config["vcf_filtered"]+"/{species}.bipassed.bt.vcf.gz",
		bisnp_passed_index=config["vcf_filtered"]+"/{species}.bipassed.bt.vcf.gz.tbi",
		vcf_stats_table_bisel=config["report_dir"]+"/filterbisnp/{species}.bisel.merged.tsv",
		vcf_stats_QUAL_biallelic_bisel=report(config["report_dir"]+"/filterbisnp/{species}_biallelic_QUAL_bisel.pdf", category="filterbisnp", subcategory="general", labels={"variant type":"biallelic", "statistic":"QUAL", "filter":"all biallelic SNPs"}),
		vcf_stats_QUAL_categories_biallelic_bisel=report(config["report_dir"]+"/filterbisnp/{species}_biallelic_QUAL__categories_bisel.pdf", category="filterbisnp", subcategory="general", labels={"variant type":"biallelic", "statistic":"QUAL_categories", "filter":"all biallelic SNPs"}),
		vcf_stats_INFO_biallelic_bisel=report(config["report_dir"]+"/filterbisnp/{species}_biallelic_INFO_bisel.pdf", category="filterbisnp", subcategory="INFO", labels={"variant type":"biallelic", "statistic":"INFO", "filter":"all biallelic SNPs"}),
		vcf_stats_GT_counts_biallelic_bisel=report(config["report_dir"]+"/filterbisnp/{species}_GT_counts_biallelic_bisel.pdf", category="filterbisnp", subcategory="Genotype counts", labels={"variant type":"biallelic", "statistic":"GT counts", "filter":"all biallelic SNPs"}),
		vcf_stats_GT_DP_biallelic_bisel=report(config["report_dir"]+"/filterbisnp/{species}_GT_DP_biallelic_bisel.pdf", category="filterbisnp", subcategory="Genotype stats", labels={"variant type":"biallelic", "statistic":"GT DP", "filter":"all biallelic SNPs"}),
		vcf_stats_GT_GQ_biallelic_bisel=report(config["report_dir"]+"/filterbisnp/{species}_GT_GQ_biallelic_bisel.pdf", category="filterbisnp", subcategory="Genotype stats", labels={"variant type":"biallelic", "statistic":"GT GQ", "filter":"all biallelic SNPs"}),
		vcf_stats_table_bipassed=config["report_dir"]+"/filterbisnp/{species}.bipassed.merged.tsv",
		vcf_stats_QUAL_categories_biallelic_bipassed=report(config["report_dir"]+"/filterbisnp/{species}_biallelic_QUAL_categories_bipassed.pdf", category="filterbisnp", subcategory="general", labels={"variant type":"biallelic", "statistic":"QUAL_categories", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"])}),
		vcf_stats_QUAL_biallelic_bipassed=report(config["report_dir"]+"/filterbisnp/{species}_biallelic_QUAL_bipassed.pdf", category="filterbisnp", subcategory="general", labels={"variant type":"biallelic", "statistic":"QUAL", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"])}),
		vcf_stats_INFO_biallelic_bipassed=report(config["report_dir"]+"/filterbisnp/{species}_biallelic_INFO_bipassed.pdf", category="filterbisnp", subcategory="INFO", labels={"variant type":"biallelic", "statistic":"INFO", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"])}),
		vcf_stats_GT_counts_biallelic_bipassed=report(config["report_dir"]+"/filterbisnp/{species}_GT_counts_biallelic_bipassed.pdf", category="filterbisnp", subcategory="Genotype counts", labels={"variant type":"biallelic", "statistic":"GT counts", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"])}),
		vcf_stats_GT_DP_biallelic_bipassed=report(config["report_dir"]+"/filterbisnp/{species}_GT_DP_biallelic_bipassed.pdf", category="filterbisnp", subcategory="Genotype stats", labels={"variant type":"biallelic", "statistic":"GT DP", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"])}),
		vcf_stats_GT_GQ_biallelic_bipassed=report(config["report_dir"]+"/filterbisnp/{species}_GT_GQ_biallelic_bipassed.pdf", category="filterbisnp", subcategory="Genotype stats", labels={"variant type":"biallelic", "statistic":"GT GQ", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"])}),
		vcf_stats_QUAL_categories_biallelic_bifilter=report(config["report_dir"]+"/filterbisnp/{species}_biallelic_QUAL_categories_bifilter.pdf", category="filterbisnp", subcategory="general", labels={"variant type":"biallelic", "statistic":"QUAL_categories", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"])}),
		vcf_stats_table_bifilter=config["report_dir"]+"/filterbisnp/{species}.bifilter.merged.tsv",
		vcf_stats_QUAL_biallelic_bifilter=report(config["report_dir"]+"/filterbisnp/{species}_biallelic_QUAL_bifilter.pdf", category="filterbisnp", subcategory="general", labels={"variant type":"biallelic", "statistic":"QUAL", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"])}),
		vcf_stats_INFO_biallelic_bifilter=report(config["report_dir"]+"/filterbisnp/{species}_biallelic_INFO_bifilter.pdf", category="filterbisnp", subcategory="INFO", labels={"variant type":"biallelic", "statistic":"INFO", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"])}),
		vcf_stats_GT_counts_biallelic_bifilter=report(config["report_dir"]+"/filterbisnp/{species}_GT_counts_biallelic_bifilter.pdf", category="filterbisnp", subcategory="Genotype counts", labels={"variant type":"biallelic", "statistic":"GT counts", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"])}),
		vcf_stats_GT_DP_biallelic_bifilter=report(config["report_dir"]+"/filterbisnp/{species}_GT_DP_biallelic_bifilter.pdf", category="filterbisnp", subcategory="Genotype stats", labels={"variant type":"biallelic", "statistic":"GT DP", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"])}),
		vcf_stats_GT_GQ_biallelic_bifilter=report(config["report_dir"]+"/filterbisnp/{species}_GT_GQ_biallelic_bifilter.pdf", category="filterbisnp", subcategory="Genotype stats", labels={"variant type":"biallelic", "statistic":"GT GQ", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"])})
	threads: 3
	resources:
		mem_mb=filter_bcftools_mem_mb,
		disk_mb=filter_bcftools_disk_mb,
		runtime=filter_bcftools_runtime
	params:
		bcftools_parse_script=workflow.source_path("../scripts/parse_bcftools_stdout.py")
	log:
		config["log_dir"]+"/filter_bisnp_bcftools_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/filter_bcftools_bisnp_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_bcftools.sh
		fi
		cd $temp_folder
		cp {params.bcftools_parse_script} $temp_folder
		cp {input.vcf_stats_table} $temp_folder
		n_sites=$(grep $'biallelic\tgeneral\tsite_count' {wildcards.species}.merged.tsv | cut -f 7)
		if [ {config[copy_large_vcfs]} -eq 1 ]
		then
			cp {input} $temp_folder
		fi
		if [  {config[combine_filtered_mv_biallelic]} -eq 1 ]
		then

			if [ {config[copy_large_vcfs]} -eq 1 ]
			then
				bcftools view --threads 1 -m3 -M3 -C 1:minor {wildcards.species}.merged.vcf.gz | python3 replace_minor.py | bgzip > {wildcards.species}.mvbiallelic.bt.vcf.gz
			else
				bcftools view --threads 1 -m3 -M3 -C 1:minor {input.vcf_input} | python3 replace_minor.py | bgzip > {wildcards.species}.mvbiallelic.bt.vcf.gz

			fi

			cp {wildcards.species}.mvbiallelic.bt.vcf.gz {output.vcf_filtered_mv_biallelic} 
			#TODO: add stat plots to output
			tabix {wildcards.species}.mvbiallelic.bt.vcf.gz
			cp {wildcards.species}.mvbiallelic.bt.vcf.gz.tbi {output.vcf_filtered_mv_biallelic_index} 
		fi


		if [ {config[copy_large_vcfs]} -eq 1 ]
		then
			bcftools view --threads 1 -m2 -M2 -v snps {wildcards.species}.merged.vcf.gz | tee >(python3 parse_bcftools_stdout.py --n_sites $n_sites --histogram_bins 50 --output {wildcards.species}.bisel --biallelic) | bgzip > {wildcards.species}.bisel.bt.vcf.gz 

		else
			bcftools view --threads 1 -m2 -M2 -v snps {input.vcf_input} | tee >(python3 parse_bcftools_stdout.py --n_sites $n_sites --histogram_bins 50 --output {wildcards.species}.bisel --biallelic) | bgzip > {wildcards.species}.bisel.bt.vcf.gz 
		fi

		cp {wildcards.species}.bisel_table.tsv {output.vcf_stats_table_bisel}
		cp {wildcards.species}.bisel_QUAL_categories_biallelic.pdf {output.vcf_stats_QUAL_categories_biallelic_bisel}
		cp {wildcards.species}.bisel_QUAL_biallelic.pdf {output.vcf_stats_QUAL_biallelic_bisel}
		cp {wildcards.species}.bisel_GT_counts_biallelic.pdf {output.vcf_stats_GT_counts_biallelic_bisel}
		cp {wildcards.species}.bisel_GT_DP_biallelic.pdf {output.vcf_stats_GT_DP_biallelic_bisel}
		cp {wildcards.species}.bisel_GT_GQ_biallelic.pdf {output.vcf_stats_GT_GQ_biallelic_bisel}
		cp {wildcards.species}.bisel_INFO_biallelic.pdf {output.vcf_stats_INFO_biallelic_bisel}

		cp {wildcards.species}.bisel.bt.vcf.gz {output.bisnp_sel} &>> {log}
		tabix {wildcards.species}.bisel.bt.vcf.gz &>> {log}
		cp {wildcards.species}.bisel.bt.vcf.gz.tbi {output.bisnp_sel_index} 
		n_sites_bisel=$(grep $'biallelic\tgeneral\tsite_count' {wildcards.species}.bisel_table.tsv | cut -f 7 )

		bcftools filter --threads 1 -m+ -s'MQ' -e'MQ<{config[MQ_less]}' {wildcards.species}.bisel.bt.vcf.gz | bcftools filter -m+ -s'QD' -e'QD<{config[QD_less]}' | bcftools filter -m+ -s'FS' -e'FS>{config[FS_more]}' | bcftools filter -m+ -s'MQRankSum' -e'MQRankSum<{config[MQRS_less]}' |  bcftools filter -m+ -s'ReadPosRankSum' -e'ReadPosRankSum<{config[RPRS_less]}' | bcftools filter -m+ -s'SOR' -e'SOR>{config[SOR_more]}' | tee >(python3 parse_bcftools_stdout.py --n_sites $n_sites_bisel --histogram_bins 50 --output {wildcards.species}.bifilter --biallelic --plot_filter) | bgzip > {wildcards.species}.bifilter.bt.vcf.gz

		cp {wildcards.species}.bifilter_table.tsv {output.vcf_stats_table_bifilter}
		cp {wildcards.species}.bifilter_QUAL_biallelic.pdf {output.vcf_stats_QUAL_biallelic_bifilter}
		cp {wildcards.species}.bifilter_QUAL_categories_biallelic.pdf {output.vcf_stats_QUAL_categories_biallelic_bifilter}
		cp {wildcards.species}.bifilter_GT_counts_biallelic.pdf {output.vcf_stats_GT_counts_biallelic_bifilter}
		cp {wildcards.species}.bifilter_GT_DP_biallelic.pdf {output.vcf_stats_GT_DP_biallelic_bifilter}
		cp {wildcards.species}.bifilter_GT_GQ_biallelic.pdf {output.vcf_stats_GT_GQ_biallelic_bifilter}
		cp {wildcards.species}.bifilter_INFO_biallelic.pdf {output.vcf_stats_INFO_biallelic_bifilter}
		cp {wildcards.species}.bifilter.bt.vcf.gz {output.bisnp_filter} 
		tabix {wildcards.species}.bifilter.bt.vcf.gz
		cp {wildcards.species}.bifilter.bt.vcf.gz.tbi {output.bisnp_filter_index} 

		bcftools view -f.,PASS {wildcards.species}.bifilter.bt.vcf.gz | tee >(python3 parse_bcftools_stdout.py --n_sites $n_sites_bisel --histogram_bins 50 --output {wildcards.species}.bipassed --biallelic)  | bgzip > {wildcards.species}.bipassed.bt.vcf.gz
		cp {wildcards.species}.bipassed_table.tsv {output.vcf_stats_table_bipassed}
		cp {wildcards.species}.bipassed_QUAL_biallelic.pdf {output.vcf_stats_QUAL_biallelic_bipassed}
		cp {wildcards.species}.bipassed_QUAL_categories_biallelic.pdf {output.vcf_stats_QUAL_categories_biallelic_bipassed}
		cp {wildcards.species}.bipassed_GT_counts_biallelic.pdf {output.vcf_stats_GT_counts_biallelic_bipassed}
		cp {wildcards.species}.bipassed_GT_DP_biallelic.pdf {output.vcf_stats_GT_DP_biallelic_bipassed}
		cp {wildcards.species}.bipassed_GT_GQ_biallelic.pdf {output.vcf_stats_GT_GQ_biallelic_bipassed}
		cp {wildcards.species}.bipassed_INFO_biallelic.pdf {output.vcf_stats_INFO_biallelic_bipassed}
		cp {wildcards.species}.bipassed.bt.vcf.gz {output.bisnp_passed} 
		tabix {wildcards.species}.bipassed.bt.vcf.gz
		cp {wildcards.species}.bipassed.bt.vcf.gz.tbi {output.bisnp_passed_index} 
		"""
#TODO filter multivariants
rule filter_bcftools_multivariants:
	input:
		config["vcf_dir"]+"/{species}.merged.vcf.gz",
		vcf_stats_table=config["report_dir"]+"/MergeSubVCFs/{species}.merged.tsv"
	output:
		#all multivariant sites
		multivar_all=config["vcf_filtered"]+"/{species}.multivar.bt.vcf.gz",
		#vcf_filtered_mvbiallelic=config["vcf_filtered"]+"/{species}.mv_biallelic.bt.vcf.gz" if config["combine_filtered_mv_biallelic"]==True else [],
		#vcf_filtered_mvbiallelic_index=config["vcf_filtered"]+"/{species}.mv_biallelic.bt.vcf.gz.tbi" if config["combine_filtered_mv_biallelic"]==True else [],
		#vcf_stats_table_mvbiallelic=config["report_dir"]+"/filtermultivariant/{species}.mv_biallelic.tsv" if config["combine_filtered_mv_biallelic"]==True else []
		#sites where a multivariate site has a low-frequency complex variant, covering a SNP
		#multivar_bisel=
	threads: 2
	resources:
		mem_mb=filter_bcftools_mem_mb,
		disk_mb=filter_bcftools_disk_mb,
		runtime=filter_bcftools_runtime
	params:
		bcftools_parse_script=workflow.source_path("../scripts/parse_bcftools_stdout.py")
	log:
		config["log_dir"]+"/filter_bisnp_bcftools_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/filter_bcftools_multivariant_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_bcftools.sh
		fi
		cd $temp_folder
		cp {params.bcftools_parse_script} $temp_folder
		if [ {config[copy_large_vcfs]} -eq 1 ]
		then
			cp {input} $temp_folder
			bcftools view --threads 1 -m3 -o {wildcards.species}.multivar.bt.vcf.gz -O z {wildcards.species}.merged.vcf.gz  &>> {log}
		else
			bcftools view --threads 1 -m3 -o {wildcards.species}.multivar.bt.vcf.gz -O z {input} &>> {log}
		"""

rule filter_bcftools_invariants:
	input:
		vcf_input=config["vcf_dir"]+"/{species}.merged.vcf.gz",
		vcf_stats_table=config["report_dir"]+"/MergeSubVCFs/{species}.merged.tsv"
	output:
		#novar_sel=config["vcf_filtered"]+"/{species}.novarsel.bt.vcf.gz",# selected novariants
		#novar_sel_index=config["vcf_filtered"]+"/{species}.novarsel.bt.vcf.gz.tbi",# selected novariants
		#novar_filt=config["vcf_filtered"]+"/{species}.novarfilt.bt.vcf.gz" ,# filtered novariants
		#novar_filt_index=config["vcf_filtered"]+"/{species}.novarfilt.bt.vcf.gz.tbi" ,# filtered novariants
		novar_passed=config["vcf_filtered"]+"/{species}.novarpass.bt.vcf.gz",  # sites that have not passed filters removed
		novar_passed_index=config["vcf_filtered"]+"/{species}.novarpass.bt.vcf.gz.tbi",  # sites that have not passed filters removed
		#vcf_stats_table_novarsel=config["report_dir"]+"/filterinvariants/{species}.novarsel.merged.tsv",
		#vcf_stats_general_invariant_novarsel=report(config["report_dir"]+"/filterinvariants/{species}_invariant_novarsel_general.pdf", category="filterinvariants", subcategory="general", labels={"variant type":"invariant", "statistic":"multiple", "filter":"all invariant sites"}),
		#vcf_stats_QUAL_invariant_novarsel=report(config["report_dir"]+"/filterinvariants/{species}_invariant_novarsel_QUAL.pdf", category="filterinvariants", subcategory="general", labels={"variant type":"invariant", "statistic":"QUAL", "filter":"all invariant sites"}),
		#vcf_stats_INFO_invariant_novarsel=report(config["report_dir"]+"/filterinvariants/{species}_invariant_novarsel_INFO.pdf", category="filterinvariants", subcategory="INFO", labels={"variant type":"invariant", "statistic":"INFO", "filter":"all invariant sites"}),
		#vcf_stats_GT_DP_invariant_novarsel=report(config["report_dir"]+"/filterinvariants/{species}_GT_DP_invariant_novarsel.pdf", category="filterinvariants", subcategory="Genotype stats", labels={"variant type":"invariant", "statistic":"GT DP", "filter":"all invariant sites"}),
		#vcf_stats_GT_RGQ_invariant_novarsel=report(config["report_dir"]+"/filterinvariants/{species}_GT_RGQ_invariant_novarsel.pdf", category="filterinvariants", subcategory="Genotype stats", labels={"variant type":"invariant", "statistic":"GT RGQ", "filter":"all invariant sites"}),
		vcf_stats_table_novarpassed=config["report_dir"]+"/filterinvariants/{species}.novarpassed.merged.tsv",
		vcf_stats_QUAL_invariant_novarpassed=report(config["report_dir"]+"/filterinvariants/{species}_invariant_novarpassed_QUAL.pdf", category="filterinvariants", subcategory="general", labels={"variant type":"invariant", "statistic":"QUAL", "filter":"QUAL<%s, GT DP<%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"], config["invariant_max_missing"]) if config["filter_qual"]==1 else "GT DP>%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"])}),
		vcf_stats_QUAL_categories_invariant_novarpassed=report(config["report_dir"]+"/filterinvariants/{species}_invariant_novarpassed_QUAL_categories.pdf", category="filterinvariants", subcategory="general", labels={"variant type":"invariant", "statistic":"QUAL_categories", "filter":"QUAL<%s, GT DP<%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"], config["invariant_max_missing"]) if config["filter_qual"]==1 else "GT DP>%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"])}),
		vcf_stats_INFO_invariant_novarpassed=report(config["report_dir"]+"/filterinvariants/{species}_invariant_novarpassed_INFO.pdf", category="filterinvariants", subcategory="INFO", labels={"variant type":"invariant", "statistic":"INFO", "filter":"QUAL<%s, GT DP<%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"], config["invariant_max_missing"]) if config["filter_qual"]==1 else "GT DP>%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"])}),
		vcf_stats_GT_count_invariant_novarpassed=report(config["report_dir"]+"/filterinvariants/{species}_GT_count_invariant_novarpassed.pdf", category="filterinvariants", subcategory="Genotype stats", labels={"variant type":"invariant", "statistic":"GT count", "filter":"QUAL<%s, GT DP<%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"], config["invariant_max_missing"]) if config["filter_qual"]==1 else "GT DP>%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"])}),
		vcf_stats_GT_DP_invariant_novarpassed=report(config["report_dir"]+"/filterinvariants/{species}_GT_DP_invariant_novarpassed.pdf", category="filterinvariants", subcategory="Genotype stats", labels={"variant type":"invariant", "statistic":"GT DP", "filter":"QUAL<%s, GT DP<%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"], config["invariant_max_missing"]) if config["filter_qual"]==1 else "GT DP>%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"])}),
		vcf_stats_GT_RGQ_invariant_novarpassed=report(config["report_dir"]+"/filterinvariants/{species}_GT_RGQ_invariant_novarpassed.pdf", category="filterinvariants", subcategory="Genotype stats", labels={"variant type":"invariant", "statistic":"GT RGQ", "filter":"QUAL<%s, GT DP<%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"], config["invariant_max_missing"]) if config["filter_qual"]==1 else "GT DP>%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"])})
	threads: 4
	resources:
		mem_mb=filter_bcftools_mem_mb,
		disk_mb=filter_bcftools_disk_mb,
		runtime=filter_bcftools_runtime
	params:
		bcftools_parse_script=workflow.source_path("../scripts/parse_bcftools_stdout.py")
	log:
		config["log_dir"]+"/filter_bcftools_invariants_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/filter_bcftools_invariants_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_bcftools.sh
		fi
		cp {params.bcftools_parse_script} $temp_folder
		cp {input.vcf_stats_table} $temp_folder 
		cd $temp_folder
		n_sites=$(grep $'invariant\tgeneral\tsite_count' {wildcards.species}.merged.tsv | cut -f 7)
		if [ {config[copy_large_vcfs]} -eq 1 ]
		then
			cp {input.vcf_input} $temp_folder
		else
			ln -s {input.vcf_input} {wildcards.species}.merged.vcf.gz
		fi
		#TODO: test if this works!
		if [ {config[filter_qual]} -eq 1 ]
		then
			bcftools view --threads 1 -C 0 {wildcards.species}.merged.vcf.gz | bcftools filter -m+ -s+ -e'QUAL<{config[invariantQUAL_less]} | QUAL=\".\"' | bcftools view -f.,PASS| bcftools filter --threads 1 -i 'FMT/DP>{config[invariant_min_depth]}' --set-GTs . | bcftools view -i 'F_MISSING<{config[invariant_max_missing]}'| tee >(python3 parse_bcftools_stdout.py --n_sites $n_sites --histogram_bins 50 --output {wildcards.species}.novarpassed --invariants) | bgzip > {wildcards.species}.novarpassed.bt.vcf.gz
		else
			bcftools view --threads 1 -C 0 {wildcards.species}.merged.vcf.gz | bcftools filter -i 'FMT/DP>{config[invariant_min_depth]}' --set-GTs . | bcftools view -i 'F_MISSING<{config[invariant_max_missing]}' | tee >(python3 parse_bcftools_stdout.py --n_sites $n_sites --histogram_bins 50 --output {wildcards.species}.novarpassed --invariants) | bgzip > {wildcards.species}.novarpassed.bt.vcf.gz
		fi
		cp {wildcards.species}.novarpassed_table.tsv {output.vcf_stats_table_novarpassed}
		cp {wildcards.species}.novarpassed_QUAL_invariant.pdf {output.vcf_stats_QUAL_invariant_novarpassed}
		cp {wildcards.species}.novarpassed_QUAL_categories_invariant.pdf {output.vcf_stats_QUAL_categories_invariant_novarpassed}
		cp {wildcards.species}.novarpassed_GT_counts_invariant.pdf {output.vcf_stats_GT_count_invariant_novarpassed}
		cp {wildcards.species}.novarpassed_GT_RGQ_invariant.pdf {output.vcf_stats_GT_RGQ_invariant_novarpassed}
		cp {wildcards.species}.novarpassed_GT_DP_invariant.pdf {output.vcf_stats_GT_DP_invariant_novarpassed}
		cp {wildcards.species}.novarpassed_INFO_invariant.pdf {output.vcf_stats_INFO_invariant_novarpassed}
		cp {wildcards.species}.novarpassed.bt.vcf.gz {output.novar_passed}
		tabix {wildcards.species}.novarpassed.bt.vcf.gz &>> {log}
		cp {wildcards.species}.novarpassed.bt.vcf.gz.tbi {output.novar_passed_index} 
		"""

rule bcftools_filter_gt_snps:
	input:
		bisnp_passed=config["vcf_filtered"]+"/{species}.bipassed.bt.vcf.gz",
		bisnp_passed_index=config["vcf_filtered"]+"/{species}.bipassed.bt.vcf.gz.tbi",
		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai",
		vcf_stats_table_bipassed=config["report_dir"]+"/filterbisnp/{species}.bipassed.merged.tsv"
	output:
		bigt_dp=config["vcf_filtered"]+"/{species}.bigt.dp.bt.vcf.gz",
		bigt_dp_index=config["vcf_filtered"]+"/{species}.bigt.dp.bt.vcf.gz.tbi",
		bigt_dp_m=config["vcf_filtered"]+"/{species}.bigt.dp.m.bt.vcf.gz",
		bigt_dp_m_index=config["vcf_filtered"]+"/{species}.bigt.dp.m.bt.vcf.gz.tbi",
		vcf_stats_table_bigt_dp=config["report_dir"]+"/filterbigt/{species}.bigt.dp.merged.tsv",
		vcf_stats_QUAL_bigt_dp=report(config["report_dir"]+"/filterbigt/{species}_bigt_dp_QUAL.pdf", category="filterbigt", subcategory="general", labels={"variant type":"biallelic", "statistic":"QUAL", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"])}),
		vcf_stats_QUAL_categories_bigt_dp=report(config["report_dir"]+"/filterbigt/{species}_bigt_dp_QUAL_categories.pdf", category="filterbigt", subcategory="general", labels={"variant type":"biallelic", "statistic":"QUAL_categories", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"])}),
		vcf_stats_INFO_bigt_dp=report(config["report_dir"]+"/filterbigt/{species}_bigt_dp_INFO.pdf", category="filterbigt", subcategory="INFO", labels={"variant type":"biallelic", "statistic":"INFO", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"])}),
		vcf_stats_GT_counts_bigt_dp=report(config["report_dir"]+"/filterbigt/{species}_GT_counts_bigt_dp.pdf", category="filterbigt", subcategory="Genotype counts", labels={"variant type":"biallelic", "statistic":"GT counts", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"])}),
		vcf_stats_GT_DP_bigt_dp=report(config["report_dir"]+"/filterbigt/{species}_GT_DP_bigt_dp.pdf", category="filterbigt", subcategory="Genotype stats", labels={"variant type":"biallelic", "statistic":"GT DP", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"])}),
		vcf_stats_GT_GQ_bigt_dp=report(config["report_dir"]+"/filterbigt/{species}_GT_GQ_bigt_dp.pdf", category="filterbigt", subcategory="Genotype stats", labels={"variant type":"biallelic", "statistic":"GT GQ", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"])}),
		vcf_stats_table_bigt_dp_m=config["report_dir"]+"/filterbigt/{species}.bigt.dp.m.merged.tsv",
		vcf_stats_QUAL_bigt_dp_m=report(config["report_dir"]+"/filterbigt/{species}_bigt_dp_m_QUAL.pdf", category="filterbigt", subcategory="general", labels={"variant type":"biallelic", "statistic":"QUAL", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])}),
		vcf_stats_QUAL_categories_bigt_dp_m=report(config["report_dir"]+"/filterbigt/{species}_bigt_dp_m_QUAL_categories.pdf", category="filterbigt", subcategory="general", labels={"variant type":"biallelic", "statistic":"QUAL_categories", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])}),
		vcf_stats_INFO_bigt_dp_m=report(config["report_dir"]+"/filterbigt/{species}_bigt_dp_m_INFO.pdf", category="filterbigt", subcategory="INFO", labels={"variant type":"biallelic", "statistic":"INFO", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])}),
		vcf_stats_GT_counts_bigt_dp_m=report(config["report_dir"]+"/filterbigt/{species}_GT_counts_bigt_dp_m.pdf", category="filterbigt", subcategory="Genotype counts", labels={"variant type":"biallelic", "statistic":"GT counts", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])}),
		vcf_stats_GT_DP_bigt_dp_m=report(config["report_dir"]+"/filterbigt/{species}_GT_DP_bigt_dp_m.pdf", category="filterbigt", subcategory="Genotype stats", labels={"variant type":"biallelic", "statistic":"GT DP", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])}),
		vcf_stats_GT_GQ_bigt_dp_m=report(config["report_dir"]+"/filterbigt/{species}_GT_GQ_bigt_dp_m.pdf", category="filterbigt", subcategory="Genotype stats", labels={"variant type":"biallelic", "statistic":"GT GQ", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])})

	threads: 2
	resources:
		mem_mb=filter_bcftools_mem_mb,
		disk_mb=filter_bcftools_disk_mb,
		runtime=filter_bcftools_runtime
	params:
		bcftools_parse_script=workflow.source_path("../scripts/parse_bcftools_stdout.py")
	log:
		config["log_dir"]+"/bcftools_filter_gt_snps_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/bcftools_filter_gt_snps_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_bcftools.sh
		fi
		cp {params.bcftools_parse_script} $temp_folder
		cp {input} $temp_folder
		cd $temp_folder
		awk '{{print$1\"\\t\"$2}}' {wildcards.species}.fasta.fai > ref.sizes

		n_sites_bipassed=$(grep $'biallelic\tgeneral\tsite_count' {wildcards.species}.bipassed.merged.tsv | cut -f 7)
		bcftools filter --threads 1 -i 'FMT/DP>{config[gen_min_depth]}' --set-GTs . {wildcards.species}.bipassed.bt.vcf.gz | tee >(python3 parse_bcftools_stdout.py --n_sites $n_sites_bipassed --histogram_bins 50 --output {wildcards.species}.bigt_dp --biallelic) | bgzip > {wildcards.species}.bigt.dp.bt.vcf.gz
		cp {wildcards.species}.bigt_dp_table.tsv {output.vcf_stats_table_bigt_dp}
		cp {wildcards.species}.bigt_dp_QUAL_biallelic.pdf {output.vcf_stats_QUAL_bigt_dp}
		cp {wildcards.species}.bigt_dp_QUAL_categories_biallelic.pdf {output.vcf_stats_QUAL_categories_bigt_dp}
		cp {wildcards.species}.bigt_dp_GT_counts_biallelic.pdf {output.vcf_stats_GT_counts_bigt_dp}
		cp {wildcards.species}.bigt_dp_GT_DP_biallelic.pdf {output.vcf_stats_GT_DP_bigt_dp}
		cp {wildcards.species}.bigt_dp_GT_GQ_biallelic.pdf {output.vcf_stats_GT_GQ_bigt_dp}
		cp {wildcards.species}.bigt_dp_INFO_biallelic.pdf {output.vcf_stats_INFO_bigt_dp}
		cp {wildcards.species}.bigt.dp.bt.vcf.gz {output.bigt_dp}
		tabix {wildcards.species}.bigt.dp.bt.vcf.gz  &>> {log}
		cp {wildcards.species}.bigt.dp.bt.vcf.gz.tbi {output.bigt_dp_index}

		#Fix header
		bcftools view -h {wildcards.species}.bigt.dp.bt.vcf.gz > {wildcards.species}.bigt.dp.bt.header
		sed -i s/ID=PL,Number=G/ID=PL,Number=./g {wildcards.species}.bigt.dp.bt.header
		bcftools reheader -h {wildcards.species}.bigt.dp.bt.header -o {wildcards.species}.bigt.dp.bt.rh.vcf.gz {wildcards.species}.bigt.dp.bt.vcf.gz
		bcftools +fill-AN-AC {wildcards.species}.bigt.dp.bt.rh.vcf.gz | bcftools view -a | bcftools view -m2 -v snps -i 'F_MISSING<{config[gen_max_missing]}' | tee >(python3 parse_bcftools_stdout.py --n_sites $n_sites_bipassed --histogram_bins 50 --output {wildcards.species}.bigt_dp_m --biallelic) | bgzip > {wildcards.species}.bigt.dp.m.bt.vcf.gz 
		cp {wildcards.species}.bigt_dp_m_table.tsv {output.vcf_stats_table_bigt_dp_m}
		cp {wildcards.species}.bigt_dp_m_QUAL_biallelic.pdf {output.vcf_stats_QUAL_bigt_dp_m}
		cp {wildcards.species}.bigt_dp_m_QUAL_categories_biallelic.pdf {output.vcf_stats_QUAL_categories_bigt_dp_m}
		cp {wildcards.species}.bigt_dp_m_GT_counts_biallelic.pdf {output.vcf_stats_GT_counts_bigt_dp_m}
		cp {wildcards.species}.bigt_dp_m_GT_DP_biallelic.pdf {output.vcf_stats_GT_DP_bigt_dp_m}
		cp {wildcards.species}.bigt_dp_m_GT_GQ_biallelic.pdf {output.vcf_stats_GT_GQ_bigt_dp_m}
		cp {wildcards.species}.bigt_dp_m_INFO_biallelic.pdf {output.vcf_stats_INFO_bigt_dp_m}
		cp {wildcards.species}.bigt.dp.m.bt.vcf.gz {output.bigt_dp_m}
		tabix {wildcards.species}.bigt.dp.m.bt.vcf.gz &>> {log}
		cp {wildcards.species}.bigt.dp.m.bt.vcf.gz.tbi {output.bigt_dp_m_index}
		"""
#change dependencies: Depthmask should be based also on pre-filtered SNPs!

rule combine_bcftools:
	input:
		novar_passed=config["vcf_filtered"]+"/{species}.novarpass.bt.vcf.gz",  # sites that have not passed filters removed
		novar_passed_index=config["vcf_filtered"]+"/{species}.novarpass.bt.vcf.gz.tbi",  # sites that have not passed filters removed
		bisnp_passed=config["vcf_filtered"]+"/{species}.bigt.dp.m.bt.vcf.gz",
		bisnp_passed_index=config["vcf_filtered"]+"/{species}.bigt.dp.m.bt.vcf.gz.tbi"
	output:
		merged=config["vcf_filtered"]+"/{species}.merged.bt.vcf.gz",
		merged_index=config["vcf_filtered"]+"/{species}.merged.bt.vcf.gz.tbi"
	threads: 2
	resources:
		mem_mb=filter_bcftools_mem_mb,
		disk_mb=filter_bcftools_disk_mb,
		runtime=filter_bcftools_runtime
	log:
		config["log_dir"]+"/combine_bcftools_variants_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/combine_bcftools_variants_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_bcftools.sh
		fi
		cp {input} $temp_folder
		cd $temp_folder
		bcftools concat --threads 1 -a -o {wildcards.species}.merged.bt.vcf.gz -O z {wildcards.species}.novarpass.bt.vcf.gz {wildcards.species}.bigt.dp.m.bt.vcf.gz &>> {log}
		cp {wildcards.species}.merged.bt.vcf.gz {output.merged}	
		tabix {wildcards.species}.merged.bt.vcf.gz &>> {log}
		cp {wildcards.species}.merged.bt.vcf.gz.tbi {output.merged_index} 
		"""


#filter sites based on depth and hetmask if available
rule bcftools_mask:
	input:
		merged=config["vcf_filtered"]+"/{species}.merged.bt.vcf.gz",
		merged_index=config["vcf_filtered"]+"/{species}.merged.bt.vcf.gz.tbi",
		hetmask=config["hetmask"] if config["hetmask"]!="None" else [],
		depthmask=config["depthmask_dir"]+"/{species}_dm.bt.bed",
		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai",
		vcf_stats_table_bigt_dp_m=config["report_dir"]+"/filterbigt/{species}.bigt.dp.m.merged.tsv",
		vcf_stats_table_novarpassed=config["report_dir"]+"/filterinvariants/{species}.novarpassed.merged.tsv"
	output:
		merged_masked=config["vcf_filtered"]+"/{species}.merged.masked.bt.vcf.gz",
		merged_masked_index=config["vcf_filtered"]+"/{species}.merged.masked.bt.vcf.gz.tbi",
		vcf_stats_table_mask=config["report_dir"]+"/filterbigt/{species}.bigt.merged.bt.tsv",
		vcf_stats_QUAL_bigt_mask=report(config["report_dir"]+"/mask/{species}_bigt_mask_QUAL.pdf", category="mask", subcategory="general", labels={"variant type":"biallelic", "statistic":"QUAL", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])}),
		vcf_stats_QUAL_categories_bigt_mask=report(config["report_dir"]+"/mask/{species}_bigt_mask_QUAL_categories.pdf", category="mask", subcategory="general", labels={"variant type":"biallelic", "statistic":"QUAL_categories", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])}),
		vcf_stats_INFO_bigt_mask=report(config["report_dir"]+"/mask/{species}_bigt_mask_INFO.pdf", category="mask", subcategory="INFO", labels={"variant type":"biallelic", "statistic":"INFO", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])}),
		vcf_stats_GT_counts_bigt_mask=report(config["report_dir"]+"/mask/{species}_GT_counts_bigt_mask.pdf", category="mask", subcategory="Genotype counts", labels={"variant type":"biallelic", "statistic":"GT counts", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])}),
		vcf_stats_GT_DP_bigt_mask=report(config["report_dir"]+"/mask/{species}_GT_DP_bigt_mask.pdf", category="mask", subcategory="Genotype stats", labels={"variant type":"biallelic", "statistic":"GT DP", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])}),
		vcf_stats_GT_GQ_bigt_mask=report(config["report_dir"]+"/mask/{species}_GT_GQ_bigt_mask.pdf", category="mask", subcategory="Genotype stats", labels={"variant type":"biallelic", "statistic":"GT GQ", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])}),
		
		vcf_stats_QUAL_invariant_mask=report(config["report_dir"]+"/mask/{species}_invariant_mask_QUAL.pdf", category="mask", subcategory="general", labels={"variant type":"invariant", "statistic":"QUAL", "filter":"QUAL<%s, GT DP<%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"], config["invariant_max_missing"]) if config["filter_qual"]==1 else "GT DP>%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"])}),
		vcf_stats_QUAL_categories_invariant_mask=report(config["report_dir"]+"/mask/{species}_invariant_mask_QUAL_categories.pdf", category="mask", subcategory="general", labels={"variant type":"invariant", "statistic":"QUAL_categories", "filter":"QUAL<%s, GT DP<%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"], config["invariant_max_missing"]) if config["filter_qual"]==1 else "GT DP>%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"])}),
		vcf_stats_INFO_invariant_mask=report(config["report_dir"]+"/mask/{species}_invariant_mask_INFO.pdf", category="mask", subcategory="INFO", labels={"variant type":"invariant", "statistic":"INFO", "filter":"QUAL<%s, GT DP<%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"], config["invariant_max_missing"]) if config["filter_qual"]==1 else "GT DP>%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"])}),
		vcf_stats_GT_counts_invariant_mask=report(config["report_dir"]+"/mask/{species}_GT_counts_invariant_mask.pdf", category="mask", subcategory="Genotype stats", labels={"variant type":"invariant", "statistic":"GT counts", "filter":"QUAL<%s, GT DP<%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"], config["invariant_max_missing"]) if config["filter_qual"]==1 else "GT DP>%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"])}),
		vcf_stats_GT_DP_invariant_mask=report(config["report_dir"]+"/mask/{species}_GT_DP_invariant_mask.pdf", category="mask", subcategory="Genotype stats", labels={"variant type":"invariant", "statistic":"GT DP", "filter":"QUAL<%s, GT DP<%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"], config["invariant_max_missing"]) if config["filter_qual"]==1 else "GT DP>%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"])}),
		vcf_stats_GT_RGQ_invariant_mask=report(config["report_dir"]+"/mask/{species}_GT_RGQ_invariant_mask.pdf", category="mask", subcategory="Genotype stats", labels={"variant type":"invariant", "statistic":"GT RGQ", "filter":"QUAL<%s, GT DP<%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"], config["invariant_max_missing"]) if config["filter_qual"]==1 else "GT DP>%s, MaxMissing=%s" % (config["invariantQUAL_less"], config["invariant_min_depth"])})

	threads: 2
	resources:
		mem_mb=filter_bcftools_mem_mb,
		disk_mb=filter_bcftools_disk_mb,
		runtime=filter_bcftools_runtime
	params:
		bcftools_parse_script=workflow.source_path("../scripts/parse_bcftools_stdout.py")
	log:
		config["log_dir"]+"/bcftools_mask_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/bcftools_mask_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_bcftools.sh
		fi
		cp {params.bcftools_parse_script} $temp_folder
		cp {input} $temp_folder
		cd $temp_folder
		#get depthmask file name
		depthmask=$(awk -F/ '{{print $NF}}' <<< {input.depthmask})
		if [ "{config[hetmask]}" != "None" ]
		then
			hetmask=$(awk -F/ '{{print $NF}}' <<< {config[hetmask]})
			bedtools merge -i $depthmask -i $hetmask > exclude.bed
		else
			mv $depthmask exclude.bed
		fi
		
		#awk '{{print$1\"\\t\"$2}}' {wildcards.species}.fasta.fai > ref.sizes

		#get n_sizes

		n_sites_bigt=$(grep $'biallelic\tgeneral\tsite_count' {wildcards.species}.bigt.dp.m.merged.tsv | cut -f 7)
		n_sites_invariant=$(grep $'invariant\tgeneral\tsite_count' {wildcards.species}.novarpassed.merged.tsv | cut -f 7)
		n_sites_total=$((n_sites_bigt+n_sites_invariant))
		#bedtools complement -i exclude.bed -g ref.sizes > include.bed
		bcftools filter --threads 1 -M exclude.bed -s MASK {wildcards.species}.merged.bt.vcf.gz | bcftools view -f .,PASS | tee >(python3 parse_bcftools_stdout.py --n_sites $n_sites_total --histogram_bins 50 --output {wildcards.species}.mask --biallelic --invariants) | bgzip > {wildcards.species}.merged.masked.bt.vcf.gz 
		

		cp {wildcards.species}.mask_table.tsv {output.vcf_stats_table_mask}
		cp {wildcards.species}.mask_QUAL_biallelic.pdf {output.vcf_stats_QUAL_bigt_mask}
		cp {wildcards.species}.mask_QUAL_categories_biallelic.pdf {output.vcf_stats_QUAL_categories_bigt_mask}
		cp {wildcards.species}.mask_GT_counts_biallelic.pdf {output.vcf_stats_GT_counts_bigt_mask}
		cp {wildcards.species}.mask_GT_DP_biallelic.pdf {output.vcf_stats_GT_DP_bigt_mask}
		cp {wildcards.species}.mask_GT_GQ_biallelic.pdf {output.vcf_stats_GT_GQ_bigt_mask}
		cp {wildcards.species}.mask_INFO_biallelic.pdf {output.vcf_stats_INFO_bigt_mask}

		cp {wildcards.species}.mask_QUAL_invariant.pdf {output.vcf_stats_QUAL_invariant_mask}
		cp {wildcards.species}.mask_QUAL_categories_invariant.pdf {output.vcf_stats_QUAL_categories_invariant_mask}
		cp {wildcards.species}.mask_GT_counts_invariant.pdf {output.vcf_stats_GT_counts_invariant_mask}
		cp {wildcards.species}.mask_GT_RGQ_invariant.pdf {output.vcf_stats_GT_RGQ_invariant_mask}
		cp {wildcards.species}.mask_INFO_invariant.pdf {output.vcf_stats_INFO_invariant_mask}
		cp {wildcards.species}.mask_GT_DP_invariant.pdf {output.vcf_stats_GT_DP_invariant_mask}
		cp {wildcards.species}.merged.masked.bt.vcf.gz {output.merged_masked}
		tabix {wildcards.species}.merged.masked.bt.vcf.gz &>> {log}
		cp {wildcards.species}.merged.masked.bt.vcf.gz.tbi {output.merged_masked_index}
		"""


rule bcftools_filter_fourfold:
	input:
		merged_masked=config["vcf_filtered"]+"/{species}.merged.masked.bt.vcf.gz",
		merged_masked_index=config["vcf_filtered"]+"/{species}.merged.masked.bt.vcf.gz.tbi",
		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai",
		fourfold_sites=config["fourfold"],
		vcf_stats_table_mask=config["report_dir"]+"/filterbigt/{species}.bigt.merged.bt.tsv"
	output:
		fourfold=config["vcf_filtered"]+"/{species}.fourfold.bt.vcf.gz",
		fourfold_index=config["vcf_filtered"]+"/{species}.fourfold.bt.vcf.gz.tbi",
		vcf_stats_table_fourfold=config["report_dir"]+"/filterbigt/{species}.bigt.dp.m.merged.tsv",
		vcf_stats_QUAL_bigt_fourfold=report(config["report_dir"]+"/fourfold/{species}_bigt_fourfold_QUAL.pdf", category="fourfold", subcategory="general", labels={"variant type":"biallelic", "statistic":"QUAL", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s, fourfold" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])}),
		vcf_stats_QUAL_categories_bigt_fourfold=report(config["report_dir"]+"/fourfold/{species}_bigt_fourfold_QUAL_categories.pdf", category="fourfold", subcategory="general", labels={"variant type":"biallelic", "statistic":"QUAL_categories", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s, fourfold" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])}),
		vcf_stats_INFO_bigt_fourfold=report(config["report_dir"]+"/fourfold/{species}_bigt_fourfold_INFO.pdf", category="fourfold", subcategory="INFO", labels={"variant type":"biallelic", "statistic":"INFO", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s, fourfold" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])}),
		vcf_stats_GT_counts_bigt_fourfold=report(config["report_dir"]+"/fourfold/{species}_GT_counts_bigt_fourfold.pdf", category="fourfold", subcategory="Genotype counts", labels={"variant type":"biallelic", "statistic":"GT counts", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s, fourfold" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])}),
		vcf_stats_GT_DP_bigt_fourfold=report(config["report_dir"]+"/fourfold/{species}_GT_DP_bigt_fourfold.pdf", category="fourfold", subcategory="Genotype stats", labels={"variant type":"biallelic", "statistic":"GT DP", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s, fourfold" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])}),
		vcf_stats_GT_GQ_bigt_fourfold=report(config["report_dir"]+"/fourfold/{species}_GT_GQ_bigt_fourfold.pdf", category="fourfold", subcategory="Genotype stats", labels={"variant type":"biallelic", "statistic":"GT GQ", "filter":"MQ<%s, QD<%s, FS>%s, MQRankSum<%s, ReadPosRankSum<%s, SOR>%s, GT DP>%s, MaxMissing=%s, fourfold" % (config["MQ_less"],config["QD_less"], config["FS_more"], config["MQRS_less"], config["RPRS_less"], config["SOR_more"], config["gen_min_depth"], config["gen_max_missing"])}),
		
		vcf_stats_QUAL_invariant_fourfold=report(config["report_dir"]+"/fourfold/{species}_invariant_fourfold_QUAL.pdf", category="fourfold", subcategory="general", labels={"variant type":"invariant", "statistic":"QUAL", "filter":"QUAL<%s, GT DP<%s, MaxMissing=%s, fourfold" % (config["invariantQUAL_less"], config["invariant_min_depth"], config["invariant_max_missing"]) if config["filter_qual"]==1 else "GT DP>%s, MaxMissing=%s, fourfold" % (config["invariantQUAL_less"], config["invariant_min_depth"])}),
		vcf_stats_QUAL_categories_invariant_fourfold=report(config["report_dir"]+"/fourfold/{species}_invariant_fourfold_QUAL_categories.pdf", category="fourfold", subcategory="general", labels={"variant type":"invariant", "statistic":"QUAL_categories", "filter":"QUAL<%s, GT DP<%s, MaxMissing=%s, fourfold" % (config["invariantQUAL_less"], config["invariant_min_depth"], config["invariant_max_missing"]) if config["filter_qual"]==1 else "GT DP>%s, MaxMissing=%s, fourfold" % (config["invariantQUAL_less"], config["invariant_min_depth"])}),
		vcf_stats_INFO_invariant_fourfold=report(config["report_dir"]+"/fourfold/{species}_invariant_fourfold_INFO.pdf", category="fourfold", subcategory="INFO", labels={"variant type":"invariant", "statistic":"INFO", "filter":"QUAL<%s, GT DP<%s, MaxMissing=%s, fourfold" % (config["invariantQUAL_less"], config["invariant_min_depth"], config["invariant_max_missing"]) if config["filter_qual"]==1 else "GT DP>%s, MaxMissing=%s, fourfold" % (config["invariantQUAL_less"], config["invariant_min_depth"])}),
		vcf_stats_GT_DP_invariant_fourfold=report(config["report_dir"]+"/fourfold/{species}_GT_DP_invariant_fourfold.pdf", category="fourfold", subcategory="Genotype stats", labels={"variant type":"invariant", "statistic":"GT DP", "filter":"QUAL<%s, GT DP<%s, MaxMissing=%s, fourfold" % (config["invariantQUAL_less"], config["invariant_min_depth"], config["invariant_max_missing"]) if config["filter_qual"]==1 else "GT DP>%s, MaxMissing=%s, fourfold" % (config["invariantQUAL_less"], config["invariant_min_depth"])}),
		vcf_stats_GT_RGQ_invariant_fourfold=report(config["report_dir"]+"/fourfold/{species}_GT_RGQ_invariant_fourfold.pdf", category="fourfold", subcategory="Genotype stats", labels={"variant type":"invariant", "statistic":"GT RGQ", "filter":"QUAL<%s, GT DP<%s, MaxMissing=%s, fourfold" % (config["invariantQUAL_less"], config["invariant_min_depth"], config["invariant_max_missing"]) if config["filter_qual"]==1 else "GT DP>%s, MaxMissing=%s, fourfold" % (config["invariantQUAL_less"], config["invariant_min_depth"])})

	threads: 1
	resources:
		mem_mb=filter_bcftools_mem_mb,
		disk_mb=filter_bcftools_disk_mb,
		runtime=filter_bcftools_runtime
	params:
		bcftools_parse_script=workflow.source_path("../scripts/parse_bcftools_stdout.py")
	log:
		config["log_dir"]+"/bcftools_filter_fourfold_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/bcftools_filter_fourfold_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_bcftools.sh
		fi
		cp {input} $temp_folder
		cp {params.bcftools_parse_script} $temp_folder
		cd $temp_folder
		n_sites_bi=$(grep $'biallelic\tgeneral\tsite_count' {wildcards.species}.bigt.merged.bt.tsv | cut -f 7)
		n_sites_inv=$(grep $'invariant\tgeneral\tsite_count' {wildcards.species}.bigt.merged.bt.tsv | cut -f 7)
		n_sites_total=$((n_sites_bi+n_sites_inv))
		fourfold_sites=$(awk -F/ '{{print $NF}}' <<< {input.fourfold_sites})
		bcftools view --threads 1 -R $fourfold_sites {wildcards.species}.merged.masked.bt.vcf.gz | tee >(python3 parse_bcftools_stdout.py --n_sites $n_sites_total --histogram_bins 50 --output {wildcards.species}.fourfold --biallelic --invariants) | bgzip >{wildcards.species}.fourfold.bt.vcf.gz 
		cp {wildcards.species}.fourfold_table.tsv {output.vcf_stats_table_fourfold}
		cp {wildcards.species}.fourfold_QUAL_biallelic.pdf {output.vcf_stats_QUAL_bigt_fourfold}
		cp {wildcards.species}.fourfold_QUAL_categories_biallelic.pdf {output.vcf_stats_QUAL_categories_bigt_fourfold}
		cp {wildcards.species}.fourfold_GT_counts_biallelic.pdf {output.vcf_stats_GT_counts_bigt_fourfold}
		cp {wildcards.species}.fourfold_GT_DP_biallelic.pdf {output.vcf_stats_GT_DP_bigt_fourfold}
		cp {wildcards.species}.fourfold_GT_GQ_biallelic.pdf {output.vcf_stats_GT_GQ_bigt_fourfold}
		cp {wildcards.species}.fourfold_INFO_biallelic.pdf {output.vcf_stats_INFO_bigt_fourfold}

		cp {wildcards.species}.fourfold_QUAL_invariant.pdf {output.vcf_stats_QUAL_invariant_fourfold}
		cp {wildcards.species}.fourfold_QUAL_categories_invariant.pdf {output.vcf_stats_QUAL_categories_invariant_fourfold}
		cp {wildcards.species}.fourfold_GT_RGQ_invariant.pdf {output.vcf_stats_GT_RGQ_invariant_fourfold}
		cp {wildcards.species}.fourfold_INFO_invariant.pdf {output.vcf_stats_INFO_invariant_fourfold}
		cp {wildcards.species}.fourfold_GT_DP_invariant.pdf {output.vcf_stats_GT_DP_invariant_fourfold}

		cp {wildcards.species}.fourfold.bt.vcf.gz {output.fourfold}
		tabix {wildcards.species}.fourfold.bt.vcf.gz
		cp {wildcards.species}.fourfold.bt.vcf.gz.tbi {output.fourfold_index}
		"""
def make_depth_mask_mem_mb_bt(wildcards, attempt):
	return int(config["make_depth_mask_mem_mb"]+(config["make_depth_mask_mem_mb"]*(attempt-1)*config["repeat_mem_mb_factor"]))
def make_depth_mask_disk_mb_bt(wildcards, attempt):
	return int(config["make_depth_mask_disk_mb"]+(config["make_depth_mask_disk_mb"]*(attempt-1)*config["repeat_disk_mb_factor"]))
def make_depth_mask_runtime_bt(wildcards, attempt):
	make_depth_mask_runtime_seconds=parse_timespan(config["make_depth_mask_runtime"])
	return str(make_depth_mask_runtime_seconds+int((make_depth_mask_runtime_seconds*(attempt-1))*config["repeat_runtime_factor"]))+"s"
#make_depth_mask_runtime_cats=config["make_depth_mask_runtime"].split(":")
#	return str(int(make_depth_mask_runtime_cats[0])+int(int(make_depth_mask_runtime_cats[0])*(attempt-1)*config["repeat_runtime_factor"]))+":"+make_depth_mask_runtime_cats[1]+":"+make_depth_mask_runtime_cats[2]
#make depth mask->check manually if it makes sense
rule make_depth_mask_bt:
	input:
		merged=config["vcf_filtered"]+"/{species}.merged.bt.vcf.gz"
	output:
		#depthmask=config["depthmask"],
		depth_out=config["depthmask_dir"]+"/{species}_dm_out.bt.tsv",
		depth_counts=config["depthmask_dir"]+"/{species}_dm_counts.bt.tsv",
		depth_hist=config["depthmask_dir"]+"/{species}_dm_hist.bt.tsv",
		depth_bed=config["depthmask_dir"]+"/{species}_dm.bt.bed"
	resources:
		mem_mb=make_depth_mask_mem_mb_bt,
		disk_mb=make_depth_mask_disk_mb_bt,
		runtime=make_depth_mask_runtime_bt
	params:
		depth_mask_script=workflow.source_path("../scripts/make_depth_mask.py")
	log:
		config["log_dir"]+"/make_depth_mask_{species}.bt.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/make_depth_mask_{wildcards.species}_bt
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_bcftools.sh
		fi
		cp {input} $temp_folder
		cp {params.depth_mask_script} $temp_folder
		cd $temp_folder
		n_loci=$(zcat {wildcards.species}.merged.bt.vcf.gz | grep -c \"^[^#]\")
		echo "Ran n_loci!"
		echo $n_loci
		ls -lh
		vcf_samples=$(bcftools query -l {wildcards.species}.merged.bt.vcf.gz)
		echo $vcf_samples
		vcf_n_samples=$(echo $vcf_samples | awk '{{print NF}}')
		echo $vcf_n_samples
		n_excess_depth=$(echo '{config[depthmask_prop_excess_depth]} * '$vcf_n_samples | bc)
		echo $n_excess_depth
		min_depth=$(({config[gen_min_depth]}+1))
		python3 make_depth_mask.py -v {wildcards.species}.merged.bt.vcf.gz -o out_file.tsv -c counts_out.tsv -p hist_out.tsv -n $n_excess_depth -l $n_loci -m $min_depth
		cp out_file.tsv {output.depth_out}
		cp counts_out.tsv {output.depth_counts}
		cp hist_out.tsv {output.depth_hist}
		awk '{{print $1\"\\t\"($2 - 1)\"\\t\"$2}}' out_file.tsv > out_file.bed
		bedops -m out_file.bed > out_file_merged.bed
		cp out_file_merged.bed {output.depth_bed}
		"""
