configfile: workflow.source_path("../../config/0_indexreference.yaml")


rule index_reference:
	input:
		config["input_fasta"]
	output:
		config["fasta_dir"]+"/{species}.fasta",
		config["fasta_dir"]+"/{species}.fasta.amb",
		config["fasta_dir"]+"/{species}.fasta.ann",
		config["fasta_dir"]+"/{species}.fasta.bwt",
		config["fasta_dir"]+"/{species}.fasta.fai",
		config["fasta_dir"]+"/{species}.fasta.pac",
		config["fasta_dir"]+"/{species}.fasta.sa",
		config["fasta_dir"]+"/{species}.dict"
	
	threads: 1
	resources:
		mem_mb=config["index_mem_mb"],
		disk_mb=config["index_disk_mb"],
		runtime=config["index_runtime"]
	log:
		config["log_dir"]+"/{species}_indexref.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/index_reference
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/0_indexreference.sh
		fi
		cp {input} $temp_folder
		cd $temp_folder
		ref=$(awk -F/ '{{print $NF}}' <<< {input})
		ref_ending=$(echo -n $ref | tail -c 3)
		if [ $ref_ending = ".gz" ]
		then
			gunzip $ref
			ref=${{ref%.gz}}
		fi
		mv -n $ref {wildcards.species}.fasta
		picard CreateSequenceDictionary R={wildcards.species}.fasta O={wildcards.species}.dict &>> {log}
		samtools faidx {wildcards.species}.fasta &>> {log}
		bwa index {wildcards.species}.fasta &>> {log}
		cp * {config[fasta_dir]}
		"""
