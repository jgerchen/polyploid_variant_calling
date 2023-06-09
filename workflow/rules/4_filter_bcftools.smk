configfile: "../config/4_filter.yaml"

def filter_bcftools_mem_mb(wildcards, attempt):
	return int(config["filter_mem_mb"]+(config["filter_mem_mb"]*(attempt-1)*config["repeat_mem_mb_factor"]))
def filter_bcftools_disk_mb(wildcards, attempt):
	return int(config["filter_disk_mb"]+(config["filter_disk_mb"]*(attempt-1)*config["repeat_disk_mb_factor"]))
def filter_bcftools_runtime(wildcards, attempt):
	filter_bcftools_runtime_cats=config["filter_runtime"].split(":")
	return str(int(filter_bcftools_runtime_cats[0])+int(int(filter_bcftools_runtime_cats[0])*(attempt-1)*config["repeat_runtime_factor"]))+":"+filter_bcftools_runtime_cats[1]+":"+filter_bcftools_runtime_cats[2]

rule filter_bcftools_bisnp:
	input:
		vcf_input=config["vcf_dir"]+"/{species}.merged.vcf.gz",
	output:
		bisnp_sel=config["vcf_filtered"]+"/{species}.bisel.bt.vcf.gz",
		bisnp_sel_index=config["vcf_filtered"]+"/{species}.bisel.bt.vcf.gz.tbi",
		bisnp_filt=config["vcf_filtered"]+"/{species}.bifilt.bt.vcf.gz",
		bisnp_filt_index=config["vcf_filtered"]+"/{species}.bifilt.bt.vcf.gz.tbi",
		bisnp_passed=config["vcf_filtered"]+"/{species}.bipassed.bt.vcf.gz",
		bisnp_passed_index=config["vcf_filtered"]+"/{species}.bipassed.bt.vcf.gz.tbi"
	threads: 2
	resources:
		mem_mb=filter_bcftools_mem_mb,
		disk_mb=filter_bcftools_disk_mb,
		runtime=filter_bcftools_runtime
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
		cp {input} $temp_folder
		cd $temp_folder
		bcftools view --threads 1 -m2 -M2 -v snps -o {wildcards.species}.bisel.bt.vcf.gz -O z {wildcards.species}.merged.vcf.gz  &>> {log}
 		
		cp {wildcards.species}.bisel.bt.vcf.gz {output.bisnp_sel} &>> {log}
		tabix {wildcards.species}.bisel.bt.vcf.gz &>> {log}
		cp {wildcards.species}.bisel.bt.vcf.gz.tbi {output.bisnp_sel_index} 
		bcftools filter --threads 1 -m+ -s+ -e'MQ<{config[MQ_less]} | QD<{config[QD_less]} | FS>{config[FS_more]} | MQRankSum<{config[MQRS_less]} | ReadPosRankSum<{config[RPRS_less]} | SOR>{config[SOR_more]} | QUAL=\".\"' -o {wildcards.species}.bifilt.bt.vcf.gz -O z MedicagoChen.bisel.bt.vcf.gz &>> {log}
		cp {wildcards.species}.bifilt.bt.vcf.gz {output.bisnp_filt} &>> {log}
		tabix {wildcards.species}.bifilt.bt.vcf.gz  &>> {log}
		cp {wildcards.species}.bifilt.bt.vcf.gz.tbi {output.bisnp_filt_index} 
		bcftools view --threads 1 -f.,PASS -o {wildcards.species}.bipassed.bt.vcf.gz -O z {wildcards.species}.bifilt.bt.vcf.gz &>> {log}
		cp {wildcards.species}.bipassed.bt.vcf.gz {output.bisnp_passed} 
		tabix {wildcards.species}.bipassed.bt.vcf.gz &>> {log}
		cp {wildcards.species}.bipassed.bt.vcf.gz.tbi {output.bisnp_passed_index} 
		"""

rule filter_bcftools_invariants:
	input:
		config["vcf_dir"]+"/{species}.merged.vcf.gz",
	output:
		novar_sel=config["vcf_filtered"]+"/{species}.novarsel.bt.vcf.gz",# selected novariants
		novar_sel_index=config["vcf_filtered"]+"/{species}.novarsel.bt.vcf.gz.tbi",# selected novariants
		novar_filt=config["vcf_filtered"]+"/{species}.novarfilt.bt.vcf.gz" ,# filtered novariants
		novar_filt_index=config["vcf_filtered"]+"/{species}.novarfilt.bt.vcf.gz.tbi" ,# filtered novariants
		novar_passed=config["vcf_filtered"]+"/{species}.novarpass.bt.vcf.gz",  # sites that have not passed filters removed
		novar_passed_index=config["vcf_filtered"]+"/{species}.novarpass.bt.vcf.gz.tbi"  # sites that have not passed filters removed
	threads: 2
	resources:
		mem_mb=filter_bcftools_mem_mb,
		disk_mb=filter_bcftools_disk_mb,
		runtime=filter_bcftools_runtime
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
		cp {input} $temp_folder
		cd $temp_folder
		bcftools view --threads 1 -C 0 -o {wildcards.species}.novarsel.bt.vcf.gz -O z {wildcards.species}.merged.vcf.gz &>> {log}
		cp {wildcards.species}.novarsel.bt.vcf.gz {output.novar_sel}
		tabix {wildcards.species}.novarsel.bt.vcf.gz &>> {log}
		cp {wildcards.species}.novarsel.bt.vcf.gz.tbi {output.novar_sel_index} 
		bcftools filter --threads 1 -m+ -s+ -e'QUAL<{config[invariantQUAL_less]} | QUAL=\".\"' -o {wildcards.species}.novarfilt.bt.vcf.gz -O z MedicagoChen.novarsel.bt.vcf.gz &>> {log}
		cp {wildcards.species}.novarfilt.bt.vcf.gz {output.novar_filt}
		tabix {wildcards.species}.novarfilt.bt.vcf.gz &>> {log}
		cp {wildcards.species}.novarfilt.bt.vcf.gz.tbi {output.novar_filt_index} 
		bcftools view --threads 1 -f.,PASS -o {wildcards.species}.novarpassed.bt.vcf.gz -O z {wildcards.species}.novarfilt.bt.vcf.gz &>> {log}
		cp {wildcards.species}.novarpassed.bt.vcf.gz {output.novar_passed}
		tabix {wildcards.species}.novarpassed.bt.vcf.gz &>> {log}
		cp {wildcards.species}.novarpassed.bt.vcf.gz.tbi {output.novar_passed_index} 
		"""
rule combine_bcftools:
	input:
		novar_passed=config["vcf_filtered"]+"/{species}.novarpass.bt.vcf.gz",  # sites that have not passed filters removed
		novar_passed_index=config["vcf_filtered"]+"/{species}.novarpass.bt.vcf.gz.tbi",  # sites that have not passed filters removed
		bisnp_passed=config["vcf_filtered"]+"/{species}.bipassed.bt.vcf.gz",
		bisnp_passed_index=config["vcf_filtered"]+"/{species}.bipassed.bt.vcf.gz.tbi"
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
		bcftools concat --threads 1 -a -o {wildcards.species}.merged.bt.vcf.gz -O z {wildcards.species}.novarpass.bt.vcf.gz {wildcards.species}.bipassed.bt.vcf.gz &>> {log}
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
		depthmask=config["depthmask_dir"]+"/{species}_dm.bed",
		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai"
	output:
		merged_filtered=config["vcf_filtered"]+"/{species}.merged.filtered.bt.vcf.gz",
		merged_filtered_index=config["vcf_filtered"]+"/{species}.merged.filtered.bt.vcf.gz.tbi"
	threads: 1
	resources:
		mem_mb=filter_bcftools_mem_mb,
		disk_mb=filter_bcftools_disk_mb,
		runtime=filter_bcftools_runtime
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
		cp {input} $temp_folder
		cd $temp_folder
		#get depthmask file name
		depthmask=$(awk -F/ '{{print $NF}}' <<< {input.depthmask})
		if [ {config[hetmask]} != "None" ]
		then
			hetmask=$(awk -F/ '{{print $NF}}' <<< {input.hetmask})
			bedtools merge -i $depthmask -i $hetmask > exclude.bed
		else
			mv $depthmask exclude.bed
		fi
		
		awk '{{print$1\"\\t\"$2}}' {wildcards.species}.fasta.fai > ref.sizes
		bedtools complement -i exclude.bed -g ref.sizes > include.bed
		bcftools view --threads 1 -R include.bed -o {wildcards.species}.merged.filtered.bt.vcf.gz -O z {wildcards.species}.merged.bt.vcf.gz &>> {log}
		cp {wildcards.species}.merged.filtered.bt.vcf.gz {output.merged_filtered}
		tabix {wildcards.species}.merged.filtered.bt.vcf.gz &>> {log}
		cp {wildcards.species}.merged.filtered.bt.vcf.gz.tbi {output.merged_filtered_index}
		"""
rule bcftools_filter_gt_snps:
	input:
		bisnp_passed=config["vcf_filtered"]+"/{species}.bipassed.bt.vcf.gz",
		bisnp_passed_index=config["vcf_filtered"]+"/{species}.bipassed.bt.vcf.gz.tbi",
		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai",
		hetmask=config["hetmask"] if config["hetmask"]!="None" else [],
		#depthmask=config["depthmask"],
		depthmask=config["depthmask_dir"]+"/{species}_dm.bed",
	output:
		bisnp_passed_dp=config["vcf_filtered"]+"/{species}.bipassed.dp.bt.vcf.gz",
		bisnp_passed_dp_index=config["vcf_filtered"]+"/{species}.bipassed.dp.bt.vcf.gz.tbi",
		bisnp_passed_dp_m=config["vcf_filtered"]+"/{species}.bipassed.dp.m.bt.vcf.gz",
		bisnp_passed_dp_m_index=config["vcf_filtered"]+"/{species}.bipassed.dp.m.bt.vcf.gz.tbi"
	threads: 2
	resources:
		mem_mb=filter_bcftools_mem_mb,
		disk_mb=filter_bcftools_disk_mb,
		runtime=filter_bcftools_runtime
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
		cp {input} $temp_folder
		cd $temp_folder
		#get depthmask file name
		depthmask=$(awk -F/ '{{print $NF}}' <<< {input.depthmask})
		if [ {config[hetmask]} != "None" ]
		then
			hetmask=$(awk -F/ '{{print $NF}}' <<< {input.hetmask})
			bedtools merge -i $depthmask -i $hetmask > exclude.bed
		else
			mv $depthmask exclude.bed
		fi
		awk '{{print$1\"\\t\"$2}}' {wildcards.species}.fasta.fai > ref.sizes
		bedtools complement -i exclude.bed -g ref.sizes > include.bed
	#Todo: IMPLEMENT ALSO IN gatk with depth/het mask!
		bcftools filter --threads 1 -i 'FMT/DP>{config[gen_min_depth]}' --set-GTs . -o {wildcards.species}.bipassed.dp.bt.vcf.gz -O z {wildcards.species}.bipassed.bt.vcf.gz &>> {log}
		cp {wildcards.species}.bipassed.dp.bt.vcf.gz {output.bisnp_passed_dp}
		tabix {wildcards.species}.bipassed.dp.bt.vcf.gz  &>> {log}
		cp {wildcards.species}.bipassed.dp.bt.vcf.gz.tbi {output.bisnp_passed_dp_index}
		bcftools view --threads 1 -R include.bed -o {wildcards.species}.bipassed.dp.m.bt.vcf.gz -O z -i 'F_MISSING<{config[gen_max_missing]}' {wildcards.species}.bipassed.dp.bt.vcf.gz &>> {log}
		cp {wildcards.species}.bipassed.dp.m.bt.vcf.gz {output.bisnp_passed_dp_m}
		tabix {wildcards.species}.bipassed.dp.m.bt.vcf.gz &>> {log}
		cp {wildcards.species}.bipassed.dp.m.bt.vcf.gz.tbi {output.bisnp_passed_dp_m_index}
		"""



rule bcftools_filter_gt_fourfold:
	input:
		bisnp_passed=config["vcf_filtered"]+"/{species}.bipassed.bt.vcf.gz",
		bisnp_passed_index=config["vcf_filtered"]+"/{species}.bipassed.bt.vcf.gz.tbi",
		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai",
		hetmask=config["hetmask"] if config["hetmask"]!="None" else [],
		#depthmask=config["depthmask"],
		depthmask=config["depthmask_dir"]+"/{species}_dm.bed",
		fourfold_sites=config["fourfold"]
	output:
		fourfold_passed_dp=config["vcf_filtered"]+"/{species}.fourfold.dp.bt.vcf.gz",
		fourfold_passed_dp_index=config["vcf_filtered"]+"/{species}.fourfold.dp.bt.vcf.gz.tbi",
		fourfold_passed_dp_m=config["vcf_filtered"]+"/{species}.fourfold.dp.m.bt.vcf.gz",
		fourfold_passed_dp_m_index=config["vcf_filtered"]+"/{species}.fourfold.dp.m.bt.vcf.gz.tbi"
	threads: 1
	resources:
		mem_mb=filter_bcftools_mem_mb,
		disk_mb=filter_bcftools_disk_mb,
		runtime=filter_bcftools_runtime
	log:
		config["log_dir"]+"/bcftools_filter_gt_fourfold_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/bcftools_filter_fourfold_snps_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/4_filter_bcftools.sh
		fi
		cp {input} $temp_folder
		cd $temp_folder
		#get depthmask file name
		depthmask=$(awk -F/ '{{print $NF}}' <<< {input.depthmask})
		if [ {config[hetmask]} != "None" ]
		then
			hetmask=$(awk -F/ '{{print $NF}}' <<< {input.hetmask})
			bedtools merge -i $depthmask -i $hetmask > exclude.bed
		else
			mv $depthmask exclude.bed
		fi
		awk '{{print$1\"\\t\"$2}}' {wildcards.species}.fasta.fai > ref.sizes
		bedtools complement -i exclude.bed -g ref.sizes > include.bed
		bedtools intersect -a {input.fourfold_sites} -b include.bed > include_fourfold.bed

		bcftools filter -i 'FMT/DP>{config[gen_min_depth]}' --set-GTs . -o {wildcards.species}.bipassed.dp.bt.vcf.gz -O z {species}.fourfold.bt.vcf.gz &>> {log}
		cp {wildcards.species}.fourfold.dp.bt.vcf.gz {output.fourfold_passed_dp}
		tabix {wildcards.species}.fourfold.dp.bt.vcf.gz  &>> {log}
		cp {wildcards.species}.fourfold.dp.bt.vcf.gz.tbi {output.fourfold_passed_dp_index}
		bcftools view --threads 1 -R include_fourfold.bed -o {wildcards.species}.fourfold.dp.m.bt.vcf.gz -O z -i 'F_MISSING<{config[gen_max_missing]}' {wildcards.species}.fourfold.dp.bt.vcf.gz &>> {log}
		cp {wildcards.species}.fourfold.dp.m.bt.vcf.gz {output.fourfold_passed_dp_m}
		tabix {wildcards.species}.fourfold.dp.m.bt.vcf.gz &>> {log}
		cp {wildcards.species}.fourfold.dp.m.bt.vcf.gz.tbi {output.fourfold_passed_dp_m_index}
		"""
