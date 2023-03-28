configfile: "../config/3_genotypeGVCF.yaml"

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
		mem_mb=96000,
		disk_mb=50000,
		runtime="24:00:00"
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
		$GATK4 --java-options \"-Xmx90G\" GenotypeGVCFs  -R {wildcards.species}.fasta -V gendb://{wildcards.species}_{wildcards.sub}_GenomicsDB -L $sub_interval -O {wildcards.species}_{wildcards.sub}.vcf.gz --tmp-dir tmp --include-non-variant-sites  &>> {log}
		cp {wildcards.species}_{wildcards.sub}.vcf.gz {output} 
		"""

#def get_intervals(wildcards):
#	return [config["gvcf_dir"]+"/"+wildcards.species+"_"+u+".gvcf.gz" for u in interval_list]


rule MergeSubVCFs:
	input:
		out_vcf=expand(config["vcf_dir"]+"/{{species}}_{sub}.vcf.gz", sub=interval_list)
	output:
		vcf_out=config["vcf_dir"]+"/{species}.merged.vcf.gz",
		vcf_out_index=config["vcf_dir"]+"/{species}.merged.vcf.gz.tbi"
	threads: 1
	resources:
		mem_mb=16000,
		disk_mb=20000,
		runtime="4:00:00"
	log:
		config["log_dir"]+"/MergeSubVCFs_{species}.log"
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
		java -jar $PICARD MergeVcfs I=input_files.list O={wildcards.species}.merged.vcf.gz &>> {log}
		cp {wildcards.species}.merged.vcf.gz {output.vcf_out}
		$GATK4 IndexFeatureFile -I {wildcards.species}.merged.vcf.gz  &>> {log}
		cp {wildcards.species}.merged.vcf.gz.tbi {output.vcf_out_index}
		"""




