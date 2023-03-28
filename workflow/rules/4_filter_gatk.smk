configfile: "../config/4_filter.yaml"

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
		mem_mb=32000,
		disk_mb=20000,
		runtime="24:00:00"
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
		#TODO: replace GATK standard filter values with values from config file
		$GATK4 VariantFiltration -R {wildcards.species}.fasta -V {wildcards.species}.bisel.vcf.gz -O {wildcards.species}.bifilt.vcf.gz --filter-name 'QD' --filter-expression 'QD < 2.0' --filter-name 'FS' --filter-expression 'FS > 60.0'  --filter-name 'MQ' --filter-expression 'MQ < 40.0' --filter-name 'MQRS' --filter-expression 'MQRankSum < -12.5' --filter-name 'RPRS' --filter-expression 'ReadPosRankSum < -8.0' --filter-name 'SOR' --filter-expression 'SOR > 3.0' &>> {log}
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
		mem_mb=32000,
		disk_mb=20000,
		runtime="24:00:00"
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
		$GATK4 VariantFiltration -R {wildcards.species}.fasta -V {wildcards.species}.novarsel.vcf.gz -O {wildcards.species}.novarfilt.vcf.gz --filter-name 'QUAL' --filter-expression 'QUAL < 15.0' &>> {log}
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
		mem_mb=16000,
		disk_mb=20000,
		runtime="24:00:00"
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
		depthmask=config["depthmask"],
		ref_fasta=config["fasta_dir"]+"/{species}.fasta",
		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai",
		ref_fasta_dict=config["fasta_dir"]+"/{species}.dict"
	output:
		merged_filtered=config["vcf_filtered"]+"/{species}.merged.filtered.vcf.gz",
		merged_filtered_index=config["vcf_filtered"]+"/{species}.mergedfiltered.vcf.gz.tbi"
	threads: 1
	resources:
		mem_mb=16000,
		disk_mb=20000,
		runtime="24:00:00"
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
		depthmask=config["depthmask"],
		ref_fasta_dict=config["fasta_dir"]+"/{species}.dict"
	output:
		bisnp_passed_dp=config["vcf_filtered"]+"/{species}.bipassed.dp"+str(config["gen_min_depth"])+".vcf.gz",
		bisnp_passed_dp_index=config["vcf_filtered"]+"/{species}.bipassed.dp"+str(config["gen_min_depth"])+".vcf.gz.tbi",
		bisnp_passed_dp_nc=config["vcf_filtered"]+"/{species}.bipassed.dp"+str(config["gen_min_depth"])+"nc.vcf.gz",
		bisnp_passed_dp_nc_index=config["vcf_filtered"]+"/{species}.bipassed.dp"+str(config["gen_min_depth"])+"nc.vcf.gz.tbi",
		bisnp_passed_dp_nc_m=config["vcf_filtered"]+"/{species}.bipassed.dp"+str(config["gen_min_depth"])+"nc.m"+str(config["gen_max_missing"])+".vcf.gz",
		bisnp_passed_dp_nc_m_index=config["vcf_filtered"]+"/{species}.bipassed.dp"+str(config["gen_min_depth"])+"nc.m"+str(config["gen_max_missing"])+".vcf.gz.tbi"
	threads: 1
	resources:
		mem_mb=32000,
		disk_mb=40000,
		runtime="24:00:00"
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
		
		$GATK4 VariantFiltration -R {wildcards.species}.fasta -V {wildcards.species}.bipassed.filtered.vcf.gz -O {wildcards.species}.bipassed.dp.vcf.gz --genotype-filter-expression \"DP < {config[gen_min_depth]}\" --genotype-filter-name 'DP' &>> {log}
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
		depthmask=config["depthmask"],
		ref_fasta_dict=config["fasta_dir"]+"/{species}.dict",
		fourfold_sites=config["fourfold"] if config["fourfold"]!="None" else []
	output:
		fourfold_filtered=config["vcf_filtered"]+"/{species}.fourfold.filtered.vcf.gz",
		fourfold_filtered_index=config["vcf_filtered"]+"/{species}.fourfold.filtered.vcf.gz.tbi",
		fourfold_bi_dp=config["vcf_filtered"]+"/{species}.bi.fourfold.dp"+str(config["gen_min_depth"])+".vcf.gz",
		fourfold_bi_dp_index=config["vcf_filtered"]+"/{species}.bi.fourfold.dp"+str(config["gen_min_depth"])+".vcf.gz.tbi",
		fourfold_bi_dp_nc=config["vcf_filtered"]+"/{species}.bi.fourfold.dp"+str(config["gen_min_depth"])+"nc.vcf.gz",
		fourfold_bi_dp_nc_index=config["vcf_filtered"]+"/{species}.bi.fourfold.dp"+str(config["gen_min_depth"])+"nc.vcf.gz.tbi",
		fourfold_bi_dp_nc_m=config["vcf_filtered"]+"/{species}.bi.fourfold.dp"+str(config["gen_min_depth"])+"nc.m"+str(config["gen_max_missing"])+".vcf.gz",
		fourfold_bi_dp_nc_m_index=config["vcf_filtered"]+"/{species}.bi.fourfold.dp"+str(config["gen_min_depth"])+"nc.m"+str(config["gen_max_missing"])+".vcf.gz.tbi"
	threads: 1
	resources:
		mem_mb=32000,
		disk_mb=40000,
		runtime="24:00:00"
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
		$GATK4 VariantFiltration -R {wildcards.species}.fasta -V {wildcards.species}.bi.fourfold.filtered.vcf.gz -O {wildcards.species}.bi.fourfold.dp.vcf.gz --genotype-filter-expression \"DP < {config[gen_min_depth]}\" --genotype-filter-name 'DP' &>> {log}
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

#make depth mask->check manually if it makes sense
rule make_depth_mask:
	input:
		merged=config["vcf_filtered"]+"/{species}.merged.vcf.gz"
	output:
		#depthmask=config["depthmask"],
		depth_out=config["depthmask_dir"]+"/{species}_dm_out.tsv",
		depth_counts=config["depthmask_dir"]+"/{species}_dm_counts.tsv",
		depth_hist=config["depthmask_dir"]+"/{species}_dm_hist.tsv",
		depth_bed=config["depthmask_dir"]+"/{species}_dm.bed"
	resources:
		mem_mb=32000,
		disk_mb=40000,
		runtime="24:00:00"
	log:
		config["log_dir"]+"/make_depth_mask_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/make_depth_mask_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_GATK.sh
		fi
		cp {input} $temp_folder
		cp scripts/make_depth_mask.py $temp_folder
		cd $temp_folder
		n_loci=$(zcat {wildcards.species}.merged.vcf.gz | grep -c \"^[^#]\")
		python3 make_depth_mask.py -v {wildcards.species}.merged.vcf.gz -o out_file.tsv -c counts_out.tsv -p hist_out.tsv -n {config[depthmask_n]} -l $n_loci
		cp out_file.tsv {output.depth_out}
		cp counts_out.tsv {output.depth_counts}
		cp hist_out.tsv {output.depth_hist}
		awk '{{print $1\"\\t\"($2 - 1)\"\\t\"$2}}' out_file.tsv > out_file.bed
		bedops -m out_file.bed > out_file_merged.bed
		cp out_file_merged.bed {output.depth_bed}
		"""
