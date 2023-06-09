configfile: "../config/3_genotypeGVCF.yaml"


def GenotypeGenomicsDBSub_mem_mb(wildcards, attempt):
	return int(config["GenotypeGenomicsDBSub_mem_mb"]+(config["GenotypeGenomicsDBSub_mem_mb"]*(attempt-1)*config["repeat_mem_mb_factor"]))
def GenotypeGenomicsDBSub_disk_mb(wildcards, attempt):
	return int(config["GenotypeGenomicsDBSub_disk_mb"]+(config["GenotypeGenomicsDBSub_disk_mb"]*(attempt-1)*config["repeat_disk_mb_factor"]))
def GenotypeGenomicsDBSub_runtime(wildcards, attempt):
	GenotypeGenomicsDBSub_runtime_cats=config["GenotypeGenomicsDBSub_runtime"].split(":")
	return str(int(GenotypeGenomicsDBSub_runtime_cats[0])+int(int(GenotypeGenomicsDBSub_runtime_cats[0])*(attempt-1)*config["repeat_runtime_factor"]))+":"+GenotypeGenomicsDBSub_runtime_cats[1]+":"+GenotypeGenomicsDBSub_runtime_cats[2]
rule GenotypeGenomicsDBSub:
	input:
		genomicsDB=config["gvcf_dir"]+"/{species}_{sub}_GenomicsDB",
		ref_fasta=config["fasta_dir"]+"/{species}.fasta",
		ref_fasta_dict=config["fasta_dir"]+"/{species}.dict",
		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai",
		sub_interval_list=config["sub_intervals"]
	output:
		out_vcf=config["vcf_dir"]+"/{species}_{sub}.vcf.gz"
	threads: 1
	resources:
		mem_mb=GenotypeGenomicsDBSub_mem_mb,
		disk_mb=GenotypeGenomicsDBSub_disk_mb,
		runtime=GenotypeGenomicsDBSub_runtime
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
		cp {input.sub_interval_list} $temp_folder
		cd $temp_folder
		mkdir tmp
		sub_interval_list=$(awk -F/ '{{print $NF}}' <<< {input.sub_interval_list})
		sub_interval=$(awk '{{if ($1==\"{wildcards.sub}\") print $2 }}' $sub_interval_list)
		$GATK4 --java-options \"-Xmx{resources[mem_mb]}m\" GenotypeGVCFs  -R {wildcards.species}.fasta -V gendb://{wildcards.species}_{wildcards.sub}_GenomicsDB -L $sub_interval -O {wildcards.species}_{wildcards.sub}.vcf.gz --tmp-dir tmp --include-non-variant-sites  &>> {log}
		cp {wildcards.species}_{wildcards.sub}.vcf.gz {output} 
		"""
def make_depth_mask_mem_mb(wildcards, attempt):
	return int(config["make_depth_mask_mem_mb"]+(config["make_depth_mask_mem_mb"]*(attempt-1)*config["repeat_mem_mb_factor"]))
def make_depth_mask_disk_mb(wildcards, attempt):
	return int(config["make_depth_mask_disk_mb"]+(config["make_depth_mask_disk_mb"]*(attempt-1)*config["repeat_disk_mb_factor"]))
def make_depth_mask_runtime(wildcards, attempt):
	make_depth_mask_runtime_cats=config["make_depth_mask_runtime"].split(":")
	return str(int(make_depth_mask_runtime_cats[0])+int(int(make_depth_mask_runtime_cats[0])*(attempt-1)*config["repeat_runtime_factor"]))+":"+make_depth_mask_runtime_cats[1]+":"+make_depth_mask_runtime_cats[2]
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
		mem_mb=make_depth_mask_mem_mb,
		disk_mb=make_depth_mask_disk_mb,
		runtime=make_depth_mask_runtime
	log:
		config["log_dir"]+"/make_depth_mask_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/make_depth_mask_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/3_genotypeGVCF.sh
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

def MergeSubVCFsbcftools_mem_mb(wildcards, attempt):
	return int(config["MergeSubVCFsbcftools_mem_mb"]+(config["MergeSubVCFsbcftools_mem_mb"]*(attempt-1)*config["repeat_mem_mb_factor"]))
def MergeSubVCFsbcftools_disk_mb(wildcards, attempt):
	return int(config["MergeSubVCFsbcftools_disk_mb"]+(config["MergeSubVCFsbcftools_disk_mb"]*(attempt-1)*config["repeat_disk_mb_factor"]))
def MergeSubVCFsbcftools_runtime(wildcards, attempt):
	MergeSubVCFsbcftools_runtime_cats=config["MergeSubVCFsbcftools_runtime"].split(":")
	return str(int(MergeSubVCFsbcftools_runtime_cats[0])+int(int(MergeSubVCFsbcftools_runtime_cats[0])*(attempt-1)*config["repeat_runtime_factor"]))+":"+MergeSubVCFsbcftools_runtime_cats[1]+":"+MergeSubVCFsbcftools_runtime_cats[2]
rule MergeSubVCFsbcftools:
	input:
		out_vcf=expand(config["vcf_dir"]+"/{{species}}_{sub}.vcf.gz", sub=interval_list)
	output:
		vcf_out=config["vcf_dir"]+"/{species}.merged.vcf.gz",
		vcf_out_index=config["vcf_dir"]+"/{species}.merged.vcf.gz.tbi"
	threads: 1
	resources:
		mem_mb=MergeSubVCFsbcftools_mem_mb,
		disk_mb=MergeSubVCFsbcftools_disk_mb,
		runtime=MergeSubVCFsbcftools_runtime
	log:
		config["log_dir"]+"/MergeSubVCFsbcftools_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/MergeSubVCFs_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/3_genotypeGVCF.sh
		fi
		cp {input} $temp_folder
		cd $temp_folder
		ls *.vcf.gz > input_files.list
		bcftools concat -f input_files.list -n -O z -o {wildcards.species}.merged.bt.vcf.gz &>> {log}
		cp {wildcards.species}.merged.bt.vcf.gz {output.vcf_out}
		tabix {wildcards.species}.merged.bt.vcf.gz 
		cp {wildcards.species}.merged.bt.vcf.gz.tbi {output.vcf_out_index}
		"""
