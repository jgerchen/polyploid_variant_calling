configfile: "../config/4_filter.yaml"
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
	output:
		bisnp_sel=config["vcf_filtered"]+"/{species}.bisel.bt.vcf.gz",
		bisnp_sel_index=config["vcf_filtered"]+"/{species}.bisel.bt.vcf.gz.tbi",
		#bisnp_filt=config["vcf_filtered"]+"/{species}.bifilt.bt.vcf.gz",
		#bisnp_filt_index=config["vcf_filtered"]+"/{species}.bifilt.bt.vcf.gz.tbi",
		bisnp_passed=config["vcf_filtered"]+"/{species}.bipassed.bt.vcf.gz",
		bisnp_passed_index=config["vcf_filtered"]+"/{species}.bipassed.bt.vcf.gz.tbi",
		plot_MQ=report(config["report_dir"]+"/filter_bcftools_bisnp/{species}_biSNP_MQ.pdf", category="filter_bcftools_bisnp", labels={"Filtered":"Biallelic SNPs","Statistic":"MQ"}),
		plot_QD=report(config["report_dir"]+"/filter_bcftools_bisnp/{species}_biSNP_QD.pdf", category="filter_bcftools_bisnp", labels={"Filtered":"Biallelic SNPs","Statistic":"QD"}),
		plot_FS=report(config["report_dir"]+"/filter_bcftools_bisnp/{species}_biSNP_FS.pdf", category="filter_bcftools_bisnp", labels={"Filtered":"Biallelic SNPs","Statistic":"FS"}),
		plot_MQRS=report(config["report_dir"]+"/filter_bcftools_bisnp/{species}_biSNP_MQRankSum.pdf", category="filter_bcftools_bisnp", labels={"Filtered":"Biallelic SNPs","Statistic":"MQRankSum"}),
		plot_RPRS=report(config["report_dir"]+"/filter_bcftools_bisnp/{species}_biSNP_ReadPosRankSum.pdf", category="filter_bcftools_bisnp", labels={"Filtered":"Biallelic SNPs","Statistic":"ReadPosRankSum"}),
		plot_SOR=report(config["report_dir"]+"/filter_bcftools_bisnp/{species}_biSNP_SOR.pdf", category="filter_bcftools_bisnp", labels={"Filtered":"Biallelic SNPs","Statistic":"SOR"})

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
		cp scripts/plot_SNPs_GATK_best_practices.R $temp_folder 
		cd $temp_folder

		if [ {config[copy_large_vcfs]} -eq 1 ]
		then
			cp {input} $temp_folder
			bcftools view --threads 1 -m2 -M2 -v snps -o {wildcards.species}.bisel.bt.vcf.gz -O z {wildcards.species}.merged.vcf.gz  &>> {log}
		else
			bcftools view --threads 1 -m2 -M2 -v snps -o {wildcards.species}.bisel.bt.vcf.gz -O z {input}  &>> {log}
		fi
		cp {wildcards.species}.bisel.bt.vcf.gz {output.bisnp_sel} &>> {log}
		tabix {wildcards.species}.bisel.bt.vcf.gz &>> {log}
		cp {wildcards.species}.bisel.bt.vcf.gz.tbi {output.bisnp_sel_index} 
		#run bcftools stats for unfiltered SNPs
		#TODO: run stats plot script
		bcftools stats -d 1,10000,10 {wildcards.species}.bisel.bt.vcf.gz >  {wildcards.species}.bisel.bt.stats

		grep "^SN" {wildcards.species}.bisel.bt.stats > {wildcards.species}.SN.stats
		grep "^AF" {wildcards.species}.bisel.bt.stats > {wildcards.species}.AF.stats
		grep "^QUAL" {wildcards.species}.bisel.bt.stats > {wildcards.species}.QUALY.stats
		grep "^DP" {wildcards.species}.bisel.bt.stats > {wildcards.species}.QUALY.stats

		#TODO: make script?
		#Rscript plot_bcftools_stats.R   bcftools_stats_plot.pdf 
		#run bcftools query for GATK best practices and plot
		bcftools query -f '%CHROM\t%POS\t%QUAL\t%INFO/MQ\t%INFO/QD\t%INFO/FS\t%INFO/MQRankSum\t%INFO/ReadPosRankSum\t%INFO/SOR\n' {wildcards.species}.bisel.bt.vcf.gz > {wildcards.species}.bisel.query 
		Rscript plot_SNPs_GATK_best_practices.R {wildcards.species}.bisel.query {config[MQ_less]} MQ_plot.pdf {config[QD_less]} QD_plot.pdf {config[FS_more]} FS_plot.pdf {config[MQRS_less]} MQRS_plot.pdf  {config[RPRS_less]} RPRS_plot.pdf {config[SOR_more]} SOR_plot.pdf 
		cp MQ_plot.pdf {output.plot_MQ}
		cp QD_plot.pdf {output.plot_QD}
		cp FS_plot.pdf {output.plot_FS}
		cp MQRS_plot.pdf {output.plot_MQRS}
		cp RPRS_plot.pdf {output.plot_RPRS}
		cp SOR_plot.pdf {output.plot_SOR}

		bcftools filter -m+ -s+ -e'MQ<{config[MQ_less]} | QD<{config[QD_less]} | FS>{config[FS_more]} | MQRankSum<{config[MQRS_less]} | ReadPosRankSum<{config[RPRS_less]} | SOR>{config[SOR_more]} | QUAL=\".\"'  {wildcards.species}.bisel.bt.vcf.gz |  bcftools view --threads 1 -f.,PASS -o {wildcards.species}.bipassed.bt.vcf.gz -O z &>> {log}
		#run bcftools stats for filtered SNPs again

		cp {wildcards.species}.bipassed.bt.vcf.gz {output.bisnp_passed} 
		tabix {wildcards.species}.bipassed.bt.vcf.gz &>> {log}
		cp {wildcards.species}.bipassed.bt.vcf.gz.tbi {output.bisnp_passed_index} 
		"""
#TODO filter multivariants
rule filter_bcftools_multivariants:
	input:
		config["vcf_dir"]+"/{species}.merged.vcf.gz"
	output:
		#all multivariant sites
		multivar_all=config["vcf_filtered"]+"/{species}.multivar.bt.vcf.gz",
		#sites where a multivariate site has a low-frequency complex variant, covering a SNP
		#multivar_bisel=
	threads: 2
	resources:
		mem_mb=filter_bcftools_mem_mb,
		disk_mb=filter_bcftools_disk_mb,
		runtime=filter_bcftools_runtime
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
		if [ {config[copy_large_vcfs]} -eq 1 ]
		then
			cp {input} $temp_folder
			bcftools view --threads 1 -m3 -o {wildcards.species}.multivar.bt.vcf.gz -O z {wildcards.species}.merged.vcf.gz  &>> {log}
		else
			bcftools view --threads 1 -m3 -o {wildcards.species}.multivar.bt.vcf.gz -O z {input} &>> {log}
		#select multi-allelic sites (either type)
		#split multiallelic sites, retain only bi-allelic SNP and set other types of variants to no-call
		#bcftools norm 	
		"""

rule filter_bcftools_invariants:
	input:
		config["vcf_dir"]+"/{species}.merged.vcf.gz",
	output:
		novar_sel=config["vcf_filtered"]+"/{species}.novarsel.bt.vcf.gz",# selected novariants
		novar_sel_index=config["vcf_filtered"]+"/{species}.novarsel.bt.vcf.gz.tbi",# selected novariants
		#novar_filt=config["vcf_filtered"]+"/{species}.novarfilt.bt.vcf.gz" ,# filtered novariants
		#novar_filt_index=config["vcf_filtered"]+"/{species}.novarfilt.bt.vcf.gz.tbi" ,# filtered novariants
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
		cd $temp_folder
		if [ {config[copy_large_vcfs]} -eq 1 ]
		then
			cp {input} $temp_folder
			bcftools view --threads 1 -C 0 -o {wildcards.species}.novarsel.bt.vcf.gz -O z {wildcards.species}.merged.vcf.gz &>> {log}
		else
			bcftools view --threads 1 -C 0 -o {wildcards.species}.novarsel.bt.vcf.gz -O z {input} &>> {log}
		fi
		bcftools stats {wildcards.species}.novarsel.bt.vcf.gz >  {wildcards.species}.novarsel.bt.stats
		#TODO: make stats plot script!
		#bcftools query -f '%CHROM\t%POS\t%ALT\t%QUAL\t%INFO/DP\n' {wildcards.species}.novarsel.bt.vcf.gz > {wildcards.species}.novarsel.query 

		cp {wildcards.species}.novarsel.bt.vcf.gz {output.novar_sel}
		tabix {wildcards.species}.novarsel.bt.vcf.gz &>> {log}
		cp {wildcards.species}.novarsel.bt.vcf.gz.tbi {output.novar_sel_index}
		#TODO: test if this works!
		if [ {config[filter_qual]} -eq 1 ]
		then
			bcftools filter --threads 1 -m+ -s+ -e'QUAL<{config[invariantQUAL_less]} | QUAL=\".\"' {wildcards.species}.novarsel.bt.vcf.gz | bcftools view --threads 1 -f.,PASS| bcftools filter --threads 1 -i 'FMT/DP>{config[invariant_min_depth]}' --set-GTs . | bcftools view--threads 1 -i 'F_MISSING<{config[invariant_max_missing]}' -o {wildcards.species}.novarpassed.bt.vcf.gz -O z &>> {log}
		else
			bcftools filter --threads 1 -i 'FMT/DP>{config[invariant_min_depth]}' --set-GTs . {wildcards.species}.novarsel.bt.vcf.gz | bcftools view--threads 1 -i 'F_MISSING<{config[invariant_max_missing]}' -o {wildcards.species}.novarpassed.bt.vcf.gz -O z &>> {log}
		fi
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
		bcftools view --threads 1 -R include.bed -o {wildcards.species}.bipassed.dp.m.bt.vcf.gz -O z -i 'F_MISSING<{config[gen_max_missing]} & AN>1' {wildcards.species}.bipassed.dp.bt.vcf.gz &>> {log}
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
		cp scripts/make_depth_mask.py $temp_folder
		cd $temp_folder
		n_loci=$(zcat {wildcards.species}.merged.bt.vcf.gz | grep -c \"^[^#]\")
		python3 make_depth_mask.py -v {wildcards.species}.merged.bt.vcf.gz -o out_file.tsv -c counts_out.tsv -p hist_out.tsv -n {config[depthmask_n]} -l $n_loci
		cp out_file.tsv {output.depth_out}
		cp counts_out.tsv {output.depth_counts}
		cp hist_out.tsv {output.depth_hist}
		awk '{{print $1\"\\t\"($2 - 1)\"\\t\"$2}}' out_file.tsv > out_file.bed
		bedops -m out_file.bed > out_file_merged.bed
		cp out_file_merged.bed {output.depth_bed}
		"""
