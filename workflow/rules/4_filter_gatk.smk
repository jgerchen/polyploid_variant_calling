configfile: "../config/4_filter.yaml"

def filter_gatk_mem_mb(wildcards, attempt):
	return int(config["filter_mem_mb"]+(config["filter_mem_mb"]*(attempt-1)*config["repeat_mem_mb_factor"]))
def filter_gatk_disk_mb(wildcards, attempt):
	return int(config["filter_disk_mb"]+(config["filter_disk_mb"]*(attempt-1)*config["repeat_disk_mb_factor"]))
def filter_gatk_runtime(wildcards, attempt):
	filter_gatk_runtime_cats=config["filter_runtime"].split(":")
	return str(int(filter_gatk_runtime_cats[0])+int(int(filter_gatk_runtime_cats[0])*(attempt-1)*config["repeat_runtime_factor"]))+":"+filter_gatk_runtime_cats[1]+":"+filter_gatk_runtime_cats[2]
#rule MergeSubVCFs:
#	input:
#		out_vcf=expand(config["vcf_dir"]+"/{{species}}_{sub}.vcf.gz", sub=interval_list)
#	output:
#		vcf_out=config["vcf_dir"]+"/{species}.merged.vcf.gz",
#		vcf_out_index=config["vcf_dir"]+"/{species}.merged.vcf.gz.tbi"
#	threads: 1
#	resources:
#		mem_mb=16000,
#		disk_mb=20000,
#		runtime="4:00:00"
#	log:
#		config["log_dir"]+"/MergeSubVCFs_{species}.log"
#	shell:
#		"""
#		temp_folder={config[temp_dir]}/MergeSubVCFs_{wildcards.species}
#		mkdir -p $temp_folder
#		trap 'rm -rf $temp_folder' TERM EXIT
#		if [ {config[load_cluster_code]} -eq 1 ]
#		then
#			source {config[cluster_code_dir]}/4_filter_GATK.sh
#		fi
#		cp {input} $temp_folder
#		cd $temp_folder
#		ls *.vcf.gz > input_files.list
#		java -jar $PICARD MergeVcfs I=input_files.list O={wildcards.species}.merged.vcf.gz &>> {log}
#		cp {wildcards.species}.merged.vcf.gz {output.vcf_out}
#		$GATK4 IndexFeatureFile -I {wildcards.species}.merged.vcf.gz  &>> {log}
#		cp {wildcards.species}.merged.vcf.gz.tbi {output.vcf_out_index}
#		"""
#
rule filter_GATK_bisnp:
	input:
		vcf_input=config["vcf_dir"]+"/{species}.merged.vcf.gz",
		vcf_input_index=config["vcf_dir"]+"/{species}.merged.vcf.gz.tbi",
		ref_fasta=config["fasta_dir"]+"/{species}.fasta",
		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai",
		ref_fasta_dict=config["fasta_dir"]+"/{species}.dict"
	output:
		bisnp_sel=config["vcf_filtered"]+"/{species}.bisel.vcf.gz",
		bisnp_sel_index=config["vcf_filtered"]+"/{species}.bisel.vcf.gz.tbi",
		bisnp_filt=config["vcf_filtered"]+"/{species}.bifilt.vcf.gz",
		bisnp_filt_index=config["vcf_filtered"]+"/{species}.bifilt.vcf.gz.tbi",
		bisnp_passed=config["vcf_filtered"]+"/{species}.bipassed.vcf.gz",
		bisnp_passed_index=config["vcf_filtered"]+"/{species}.bipassed.vcf.gz.tbi"
	threads: 1
	resources:
		mem_mb=filter_gatk_mem_mb,
		disk_mb=filter_gatk_disk_mb,
		runtime=filter_gatk_runtime
	log:
		config["log_dir"]+"/filter_GATK_bisnp_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/filter_GATK_bisnp_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_GATK.sh
		fi
		cp {input} $temp_folder
		cd $temp_folder
		$GATK4 SelectVariants -R {wildcards.species}.fasta -V {wildcards.species}.merged.vcf.gz -O {wildcards.species}.bisel.vcf.gz  --restrict-alleles-to BIALLELIC -select-type SNP &>> {log}
		cp {wildcards.species}.bisel.vcf.gz {output.bisnp_sel} &>> {log}
		$GATK4 IndexFeatureFile -I {wildcards.species}.bisel.vcf.gz 
		cp {wildcards.species}.bisel.vcf.gz.tbi {output.bisnp_sel_index} 
		$GATK4 VariantFiltration -R {wildcards.species}.fasta -V {wildcards.species}.bisel.vcf.gz -O {wildcards.species}.bifilt.vcf.gz --filter-name \"QD\" --filter-expression \"QD < {config[QD_less]}\" --filter-name \"FS\" --filter-expression \"FS > {config[FS_more]}\"  --filter-name \"MQ\" --filter-expression \"MQ < {config[MQ_less]}\" --filter-name \"MQRS\" --filter-expression \"MQRankSum < {config[MQRS_less]}\" --filter-name \"RPRS\" --filter-expression \"ReadPosRankSum < {config[RPRS_less]}\" --filter-name \"SOR\" --filter-expression \"SOR > {config[SOR_more]}\" &>> {log}
		cp {wildcards.species}.bifilt.vcf.gz {output.bisnp_filt}
		$GATK4 IndexFeatureFile -I {wildcards.species}.bifilt.vcf.gz &>> {log}
		cp {wildcards.species}.bifilt.vcf.gz.tbi {output.bisnp_filt_index} 
		$GATK4 SelectVariants -R {wildcards.species}.fasta -V {wildcards.species}.bifilt.vcf.gz -O {wildcards.species}.bipassed.vcf.gz --exclude-filtered &>> {log}
		cp {wildcards.species}.bipassed.vcf.gz {output.bisnp_passed} 
		$GATK4 IndexFeatureFile -I {wildcards.species}.bipassed.vcf.gz &>> {log}
		cp {wildcards.species}.bipassed.vcf.gz.tbi {output.bisnp_passed_index} 
		"""

rule filter_GATK_invariants:
	input:
		config["vcf_dir"]+"/{species}.merged.vcf.gz",
		config["vcf_dir"]+"/{species}.merged.vcf.gz.tbi",
		ref_fasta=config["fasta_dir"]+"/{species}.fasta",
		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai",
		ref_fasta_dict=config["fasta_dir"]+"/{species}.dict"
	output:
		novar_sel=config["vcf_filtered"]+"/{species}.novarsel.vcf.gz",# selected novariants
		novar_sel_index=config["vcf_filtered"]+"/{species}.novarsel.vcf.gz.tbi",# selected novariants
		novar_filt=config["vcf_filtered"]+"/{species}.novarfilt.vcf.gz" ,# filtered novariants
		novar_filt_index=config["vcf_filtered"]+"/{species}.novarfilt.vcf.gz.tbi" ,# filtered novariants
		novar_passed=config["vcf_filtered"]+"/{species}.novarpass.vcf.gz",  # sites that have not passed filters removed
		novar_passed_index=config["vcf_filtered"]+"/{species}.novarpass.vcf.gz.tbi"  # sites that have not passed filters removed
	threads: 1
	resources:
		mem_mb=filter_gatk_mem_mb,
		disk_mb=filter_gatk_disk_mb,
		runtime=filter_gatk_runtime
	log:
		config["log_dir"]+"/filter_GATK_invariants_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/filter_GATK_invariants_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_GATK.sh
		fi
		cp {input} $temp_folder
		cd $temp_folder
		$GATK4 SelectVariants -R {wildcards.species}.fasta -V {wildcards.species}.merged.vcf.gz -O {wildcards.species}.novarsel.vcf.gz -select-type NO_VARIATION -xl-select-type MNP &>> {log}
		cp {wildcards.species}.novarsel.vcf.gz {output.novar_sel}
		$GATK4 IndexFeatureFile -I {wildcards.species}.novarsel.vcf.gz &>> {log}
		cp {wildcards.species}.novarsel.vcf.gz.tbi {output.novar_sel_index} 
		$GATK4 VariantFiltration -R {wildcards.species}.fasta -V {wildcards.species}.novarsel.vcf.gz -O {wildcards.species}.novarfilt.vcf.gz --filter-name \"QUAL\" --filter-expression \"QUAL < {config[invariantQUAL_less]}\" &>> {log}
		cp {wildcards.species}.novarfilt.vcf.gz {output.novar_filt}
		$GATK4 IndexFeatureFile -I {wildcards.species}.novarfilt.vcf.gz &>> {log}
		cp {wildcards.species}.novarfilt.vcf.gz.tbi {output.novar_filt_index} 
		$GATK4 SelectVariants -R {wildcards.species}.fasta -V {wildcards.species}.novarfilt.vcf.gz -O {wildcards.species}.novarpassed.vcf.gz --exclude-filtered &>> {log}
		cp {wildcards.species}.novarpassed.vcf.gz {output.novar_passed}
		$GATK4 IndexFeatureFile -I {wildcards.species}.novarpassed.vcf.gz &>> {log}
		cp {wildcards.species}.novarpassed.vcf.gz.tbi {output.novar_passed_index} 
		"""
rule combine_GATK_variants:
	input:
		novar_passed=config["vcf_filtered"]+"/{species}.novarpass.vcf.gz",  # sites that have not passed filters removed
		novar_passed_index=config["vcf_filtered"]+"/{species}.novarpass.vcf.gz.tbi",  # sites that have not passed filters removed
		bisnp_passed=config["vcf_filtered"]+"/{species}.bipassed.vcf.gz",
		bisnp_passed_index=config["vcf_filtered"]+"/{species}.bipassed.vcf.gz.tbi",
		ref_fasta=config["fasta_dir"]+"/{species}.fasta",
		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai",
		ref_fasta_dict=config["fasta_dir"]+"/{species}.dict"
	output:
		merged=config["vcf_filtered"]+"/{species}.merged.vcf.gz",
		merged_index=config["vcf_filtered"]+"/{species}.merged.vcf.gz.tbi"
	threads: 1
	resources:
		mem_mb=filter_gatk_mem_mb,
		disk_mb=filter_gatk_disk_mb,
		runtime=filter_gatk_runtime
	log:
		config["log_dir"]+"/combine_GATK_variants_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/combine_GATK_variants_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_GATK.sh
		fi
		cp {input} $temp_folder
		cd $temp_folder
		#do PICARD MergeVcfs
		java -jar $PICARD MergeVcfs I={wildcards.species}.novarpass.vcf.gz I={wildcards.species}.bipassed.vcf.gz O={wildcards.species}.merged.vcf.gz &>> {log}
		cp {wildcards.species}.merged.vcf.gz {output.merged}	
		$GATK4 IndexFeatureFile -I {wildcards.species}.merged.vcf.gz &>> {log}
		cp {wildcards.species}.merged.vcf.gz.tbi {output.merged_index} 
		"""


#filter sites based on depth and hetmask if available
rule GATK_mask:
	input:
		merged=config["vcf_filtered"]+"/{species}.merged.vcf.gz",
		merged_index=config["vcf_filtered"]+"/{species}.merged.vcf.gz.tbi",
		hetmask=config["hetmask"] if config["hetmask"]!="None" else [],
		depthmask=config["depthmask_dir"]+"/{species}_dm.bed",
		ref_fasta=config["fasta_dir"]+"/{species}.fasta",
		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai",
		ref_fasta_dict=config["fasta_dir"]+"/{species}.dict"
	output:
		merged_filtered=config["vcf_filtered"]+"/{species}.merged.filtered.vcf.gz",
		merged_filtered_index=config["vcf_filtered"]+"/{species}.mergedfiltered.vcf.gz.tbi"
	threads: 1
	resources:
		mem_mb=filter_gatk_mem_mb,
		disk_mb=filter_gatk_disk_mb,
		runtime=filter_gatk_runtime
	log:
		config["log_dir"]+"/GATK_mask_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/GATK_mask_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_GATK.sh
		fi
		cp {input} $temp_folder
		cd $temp_folder
		#get depthmask file name
		depthmask=$(awk -F/ '{{print $NF}}' <<< {input.depthmask})
		if [ {config[hetmask]} != "None" ]
		then
			hetmask=$(awk -F/ '{{print $NF}}' <<< {input.hetmask})
			$GATK4 VariantFiltration -R {wildcards.species}.fasta -V {wildcards.species}.merged.filtered.vcf.gz -O {wildcards.species}.merged.filtered.vcf.gz -XL $depthmask -XL $hetmask &>> {log}
		else
			$GATK4 VariantFiltration -R {wildcards.species}.fasta -V {wildcards.species}.merged.vcf.gz -O {wildcards.species}.merged.filtered.vcf.gz -XL $depthmask &>> {log}
		fi
		$GATK4 IndexFeatureFile -I {wildcards.species}.merged.filtered.vcf.gz &>> {log}
		cp {wildcards.species}.merged.filtered.vcf.gz {output.merged_filtered}
		cp {wildcards.species}.merged.filtered.vcf.gz.tbi {output.merged_filtered_index}
		"""
rule GATK_filter_gt_snps:
	input:
		bisnp_passed=config["vcf_filtered"]+"/{species}.bipassed.vcf.gz",
		bisnp_passed_index=config["vcf_filtered"]+"/{species}.bipassed.vcf.gz.tbi",
		ref_fasta=config["fasta_dir"]+"/{species}.fasta",
		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai",
		hetmask=config["hetmask"] if config["hetmask"]!="None" else [],
		#depthmask=config["depthmask"],
		depthmask=config["depthmask_dir"]+"/{species}_dm.bed",
		ref_fasta_dict=config["fasta_dir"]+"/{species}.dict"
	output:
		bisnp_passed_dp=config["vcf_filtered"]+"/{species}.bipassed.dp.vcf.gz",
		bisnp_passed_dp_index=config["vcf_filtered"]+"/{species}.bipassed.dp.vcf.gz.tbi",
		bisnp_passed_dp_nc=config["vcf_filtered"]+"/{species}.bipassed.dp.nc.vcf.gz",
		bisnp_passed_dp_nc_index=config["vcf_filtered"]+"/{species}.bipassed.dp.nc.vcf.gz.tbi",
		bisnp_passed_dp_nc_m=config["vcf_filtered"]+"/{species}.bipassed.dp.nc.m.vcf.gz",
		bisnp_passed_dp_nc_m_index=config["vcf_filtered"]+"/{species}.bipassed.dp.nc.m.vcf.gz.tbi"
	threads: 1
	resources:
		mem_mb=filter_gatk_mem_mb,
		disk_mb=filter_gatk_disk_mb,
		runtime=filter_gatk_runtime
	log:
		config["log_dir"]+"/GATK_filter_gt_snps_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/GATK_filter_gt_snps_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_GATK.sh
		fi
		cp {input} $temp_folder
		cd $temp_folder
		#get depthmask file name
		depthmask=$(awk -F/ '{{print $NF}}' <<< {input.depthmask})
		if [ {config[hetmask]} != "None" ]
		then
			hetmask=$(awk -F/ '{{print $NF}}' <<< {input.hetmask})
			$GATK4 VariantFiltration -R {wildcards.species}.fasta -V {wildcards.species}.bipassed.vcf.gz -O {wildcards.species}.bipassed.filtered.vcf.gz -XL $depthmask -XL $hetmask &>> {log}
		else
			$GATK4 VariantFiltration -R {wildcards.species}.fasta -V {wildcards.species}.bipassed.vcf.gz -O {wildcards.species}.bipassed.filtered.vcf.gz -XL $depthmask &>> {log}
		fi
		$GATK4 IndexFeatureFile -I {wildcards.species}.bipassed.filtered.vcf.gz &>> {log}
		
		$GATK4 VariantFiltration -R {wildcards.species}.fasta -V {wildcards.species}.bipassed.filtered.vcf.gz -O {wildcards.species}.bipassed.dp.vcf.gz --genotype-filter-expression \"DP < {config[gen_min_depth]}\" --genotype-filter-name \"DP\" &>> {log}
		cp {wildcards.species}.bipassed.dp.vcf.gz {output.bisnp_passed_dp}
		$GATK4 IndexFeatureFile -I {wildcards.species}.bipassed.dp.vcf.gz &>> {log}
		cp {wildcards.species}.bipassed.dp.vcf.gz.tbi {output.bisnp_passed_dp_index}
		$GATK4 VariantFiltration -R {wildcards.species}.fasta -V {wildcards.species}.bipassed.dp.vcf.gz -O {wildcards.species}.bipassed.dpnc.vcf.gz --set-filtered-genotype-to-no-call &>> {log}
		cp {wildcards.species}.bipassed.dpnc.vcf.gz {output.bisnp_passed_dp_nc}
		$GATK4 IndexFeatureFile -I {wildcards.species}.bipassed.dpnc.vcf.gz &>> {log}
		cp {wildcards.species}.bipassed.dpnc.vcf.gz.tbi {output.bisnp_passed_dp_nc_index}

		$GATK4 SelectVariants -R {wildcards.species}.fasta -V {wildcards.species}.bipassed.dpnc.vcf.gz  -O  {wildcards.species}.bipassed.dpncm.vcf.gz --max-nocall-fraction {config[gen_max_missing]}  &>> {log}
		cp {wildcards.species}.bipassed.dpncm.vcf.gz {output.bisnp_passed_dp_nc_m}
		$GATK4 IndexFeatureFile -I {wildcards.species}.bipassed.dpncm.vcf.gz &>> {log}
		cp {wildcards.species}.bipassed.dpncm.vcf.gz.tbi {output.bisnp_passed_dp_nc_m_index}
		"""



rule GATK_filter_gt_fourfold:
	input:
		merged_filtered=config["vcf_filtered"]+"/{species}.merged.filtered.vcf.gz",
		bisnp_passed=config["vcf_filtered"]+"/{species}.bipassed.vcf.gz",
		bisnp_passed_index=config["vcf_filtered"]+"/{species}.bipassed.vcf.gz.tbi",
		ref_fasta=config["fasta_dir"]+"/{species}.fasta",
		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai",
		hetmask=config["hetmask"] if config["hetmask"]!="None" else [],
		#depthmask=config["depthmask"],
		depthmask=config["depthmask_dir"]+"/{species}_dm.bed",
		ref_fasta_dict=config["fasta_dir"]+"/{species}.dict",
		fourfold_sites=config["fourfold"] if config["fourfold"]!="None" else []
	output:
		fourfold_filtered=config["vcf_filtered"]+"/{species}.fourfold.filtered.vcf.gz",
		fourfold_filtered_index=config["vcf_filtered"]+"/{species}.fourfold.filtered.vcf.gz.tbi",
		fourfold_bi_dp=config["vcf_filtered"]+"/{species}.bi.fourfold.dp.vcf.gz",
		fourfold_bi_dp_index=config["vcf_filtered"]+"/{species}.bi.fourfold.dp.vcf.gz.tbi",
		fourfold_bi_dp_nc=config["vcf_filtered"]+"/{species}.bi.fourfold.dp.nc.vcf.gz",
		fourfold_bi_dp_nc_index=config["vcf_filtered"]+"/{species}.bi.fourfold.dp.nc.vcf.gz.tbi",
		fourfold_bi_dp_nc_m=config["vcf_filtered"]+"/{species}.bi.fourfold.dp.nc.m.vcf.gz",
		fourfold_bi_dp_nc_m_index=config["vcf_filtered"]+"/{species}.bi.fourfold.dp.nc.m.vcf.gz.tbi"
	threads: 1
	resources:
		mem_mb=filter_gatk_mem_mb,
		disk_mb=filter_gatk_disk_mb,
		runtime=filter_gatk_runtime
	log:
		config["log_dir"]+"/GATK_filter_gt_fourfold_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/GATK_filter_gt_fourfold_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_GATK.sh
		fi
		cp {input} $temp_folder
		cd $temp_folder
		fourfold=$(awk -F/ '{{print $NF}}' <<< {input.fourfold})
		$GATK4 VariantFiltration -R {wildcards.species}.fasta -V {wildcards.species}.merged.filtered.vcf.gz -O {wildcards.species}.fourfold.filtered.vcf.gz -L $fourfold &>> {log}
		cp {wildcards.species}.fourfold.filtered.vcf.gz {output.fourfold_filtered} 
		$GATK4 IndexFeatureFile -I {wildcards.species}.fourfold.filtered.vcf.gz &>> {log}
		cp {wildcards.species}.fourfold.filtered.vcf.gz.tbi {output.fourfold_filtered_index} 
		depthmask=$(awk -F/ '{{print $NF}}' <<< {input.depthmask})
		if [ {config[hetmask]} != "None" ]
		then
			hetmask=$(awk -F/ '{{print $NF}}' <<< {input.hetmask})
			$GATK4 VariantFiltration -R {wildcards.species}.fasta -V {wildcards.species}.bipassed.vcf.gz -O {wildcards.species}.bipassed.filtered.vcf.gz -XL $depthmask -XL $hetmask -L $fourfold &>> {log}
		else
			$GATK4 VariantFiltration -R {wildcards.species}.fasta -V {wildcards.species}.bipassed.vcf.gz -O {wildcards.species}.bi.fourfold.filtered.vcf.gz -XL $depthmask -L $fourfold &>> {log}
		fi
		$GATK4 IndexFeatureFile -I {wildcards.species}.bi.fourfold.filtered.vcf.gz &>> {log}
		$GATK4 VariantFiltration -R {wildcards.species}.fasta -V {wildcards.species}.bi.fourfold.filtered.vcf.gz -O {wildcards.species}.bi.fourfold.dp.vcf.gz --genotype-filter-expression \"DP < {config[gen_min_depth]}\" --genotype-filter-name \"DP\" &>> {log}
		cp {wildcards.species}.bi.fourfold.dp.vcf.gz {output.fourfold_bi_dp}
		$GATK4 IndexFeatureFile -I {wildcards.species}.bi.fourfold.dp.vcf.gz &>> {log}
		cp {wildcards.species}.bipassed.dp.vcf.gz.tbi {output.fourfold_bi_dp_index}
		$GATK4 VariantFiltration -R {wildcards.species}.fasta -V {wildcards.species}.bi.fourfold.dp.vcf.gz -O {wildcards.species}.bi.fourfold.dpnc.vcf.gz --set-filtered-genotype-to-no-call &>> {log}
		cp {wildcards.species}.bi.fourfold.dpnc.vcf.gz {output.fourfold_bi_dp_nc}
		$GATK4 IndexFeatureFile -I {wildcards.species}.bi.fourfold.dpnc.vcf.gz &>> {log}
		cp {wildcards.species}.bi.fourfold.dpnc.vcf.gz.tbi {config.bisnp_passed_dp_nc_index}
		$GATK4 SelectVariants -R {wildcards.species}.fasta -V {wildcards.species}.bi.fourfold.dpnc.vcf.gz  -O  {wildcards.species}.bi.fourfold.dpncm.vcf.gz --max-nocall-fraction {config[gen_max_missing]}  &>> {log}
		cp {wildcards.species}.bi.fourfold.dpncm.vcf.gz {output.fourfold_bi_dp_nc_m}
		$GATK4 IndexFeatureFile -I {wildcards.species}.bi.fourfold.dpncm.vcf.gz &>> {log}
		cp {wildcards.species}.bi.fourfold.dpncm.vcf.gz.tbi {config.bisnp_passed_dp_nc_m_index}
		"""
