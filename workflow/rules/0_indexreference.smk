configfile: "../config/0_indexreference.yaml"

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
		mem_mb=8000,
		disk_mb=4000,
		runtime="4:00:00"
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
		cp {input} $temp_folder/{wildcards.species}.fasta
		cd $temp_folder
		ref={wildcards.species}.fasta
		java -jar $PICARD CreateSequenceDictionary R=$ref O=${{ref/fasta/dict}} &>> {log}
		samtools faidx $ref &>> {log}
		bwa index $ref &>> {log}
		cp * {config[fasta_dir]}
		"""
