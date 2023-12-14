from humanfriendly import parse_timespan
configfile: "../config/2_callvars.yaml"


with open(config["sub_intervals"]) as interval_file:
	interval_list=[i.strip().split()[0] for i in interval_file]
###DO not use too slow.
#rule haplotypecaller:
#	input:
#		bam=config["bam_dir"]+"/{species}_{sample}.merged.dedup.bam",
#		bam_index=config["bam_dir"]+"/{species}_{sample}.merged.dedup.bam.bai",
#		ploidy_input=config["sample_ploidies"],
#		ref_fasta=config["fasta_dir"]+"/{species}.fasta",
#		ref_fasta_dict=config["fasta_dir"]+"/{species}.dict",
#		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai"
#	output:
#		gvcf_out=config["gvcf_dir"]+"/{species}_{sample}.gvcf.gz"
#	threads: 4
#	resources:
#		mem_mb=32000,
#		disk_mb=10000,
#		runtime="24:00:00"
#	log:
#		config["log_dir"]+"/haplotypecaller_{species}_{sample}.log"
#	shell:
#		"""
#		temp_folder={config[temp_dir]}/haplotypecaller_{wildcards.species}_{wildcards.sample}
#		mkdir -p $temp_folder
#		trap 'rm -rf $temp_folder' TERM EXIT
#		if [ {config[load_cluster_code]} -eq 1 ]
#		then
#			source {config[cluster_code_dir]}/2_callvars.sh
#		fi
#		cp {input} $temp_folder
#		cd $temp_folder
#		ploidy_file=$(awk -F/ '{{print $NF}}' <<< {input.ploidy_input})
#		s_ploidy=$(awk '{{if ($1==\"{wildcards.sample}\") print $2 }}' $ploidy_file )
#		echo "Sample ploidy is "$s_ploidy >> {log}
#		$GATK4 HaplotypeCaller -I {wildcards.species}_{wildcards.sample}.merged.dedup.bam -R {wildcards.species}.fasta -O {wildcards.species}_{wildcards.sample}.g.vcf.gz -ERC GVCF --min-base-quality-score {config[hapcaller_minbaseq]} --minimum-mapping-quality {config[hapcaller_minmapq]} -ploidy $s_ploidy -stand-call-conf 20 --pcr-indel-model NONE --max-genotype-count 350 &>> {log}
#		#removed rf BadMate
#		cp {wildcards.species}_{wildcards.sample}.g.vcf.gz {output.gvcf_out}
#		"""
#
##Does not work due to JAVA crap
#rule haplotypecaller_spark:
#	input:
#		bam=config["bam_dir"]+"/{species}_{sample}.merged.dedup.bam",
#		bam_index=config["bam_dir"]+"/{species}_{sample}.merged.dedup.bam.bai",
#		ploidy_input=config["sample_ploidies"],
#		ref_fasta=config["fasta_dir"]+"/{species}.fasta",
#		ref_fasta_dict=config["fasta_dir"]+"/{species}.dict",
#		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai"
#	output:
#		gvcf_out=config["gvcf_dir"]+"/{species}_{sample}.spark.gvcf.gz"
#	threads: 10
#	resources:
#		mem_mb=64000,
#		disk_mb=10000,
#		runtime="24:00:00"
#	log:
#		config["log_dir"]+"/haplotypecallerSpark_{species}_{sample}.log"
#	shell:
#		"""
#		temp_folder={config[temp_dir]}/haplotypecallerSpark_{wildcards.species}_{wildcards.sample}
#		mkdir -p $temp_folder
#		trap 'rm -rf $temp_folder' TERM EXIT
#		if [ {config[load_cluster_code]} -eq 1 ]
#		then
#			source {config[cluster_code_dir]}/2_callvars.sh
#		fi
#		cp {input} $temp_folder
#		cd $temp_folder
#		ploidy_file=$(awk -F/ '{{print $NF}}' <<< {input.ploidy_input})
#		s_ploidy=$(awk '{{if ($1==\"{wildcards.sample}\") print $2 }}' $ploidy_file )
#		echo "Sample ploidy is "$s_ploidy >> {log}
#		$GATK4 HaplotypeCallerSpark -I {wildcards.species}_{wildcards.sample}.merged.dedup.bam -R {wildcards.species}.fasta -O {wildcards.species}_{wildcards.sample}.g.vcf.gz -ERC GVCF --min-base-quality-score {config[hapcaller_minbaseq]} --minimum-mapping-quality {config[hapcaller_minmapq]} -ploidy $s_ploidy -stand-call-conf 20 --pcr-indel-model NONE --max-genotype-count 350 &>> {log}
#		#removed rf BadMate
#		cp {wildcards.species}_{wildcards.sample}.g.vcf.gz {output.gvcf_out}
#		"""
#use this one
def hapcallerSub_mem_mb(wildcards, attempt):
	return int(config["hapcallerSub_mem_mb"]+(config["hapcallerSub_mem_mb"]*(attempt-1)*config["repeat_mem_mb_factor"]))
def hapcallerSub_disk_mb(wildcards, attempt):
	return int(config["hapcallerSub_disk_mb"]+(config["hapcallerSub_disk_mb"]*(attempt-1)*config["repeat_disk_mb_factor"]))
def hapcallerSub_runtime(wildcards, attempt):
	hapcallerSub_runtime_seconds=parse_timespan(config["hapcallerSub_runtime"])
	return str(hapcallerSub_runtime_seconds+int((hapcallerSub_runtime_seconds*(attempt-1))*config["repeat_runtime_factor"]))+"s"
	#hapcallerSub_runtime_cats=config["hapcallerSub_runtime"].split(":")
	#return str(int(hapcallerSub_runtime_cats[0])+int(int(hapcallerSub_runtime_cats[0])*(attempt-1)*config["repeat_runtime_factor"]))+":"+hapcallerSub_runtime_cats[1]+":"+hapcallerSub_runtime_cats[2]

#TODO: extract read information from GATK haplotypeCaller output
rule hapcallerSub:
	input:
		bam=config["bam_dir"]+"/{species}_{sample}.merged.dedup.bam",
		bam_index=config["bam_dir"]+"/{species}_{sample}.merged.dedup.bam.bai",
		ploidy_input=ancient(config["sample_ploidies"]),
		sub_interval_list=config["sub_intervals"],
		ref_fasta=config["fasta_dir"]+"/{species}.fasta",
		ref_fasta_dict=config["fasta_dir"]+"/{species}.dict",
		ref_fasta_fai=config["fasta_dir"]+"/{species}.fasta.fai"
	output:
		gvcf_out=config["gvcf_dir"]+"/{species}_{sample}_{sub}.gvcf.gz"
	threads: 2
	resources:
		mem_mb=hapcallerSub_mem_mb,
		disk_mb=hapcallerSub_disk_mb,
		runtime=hapcallerSub_runtime
	log:
		config["log_dir"]+"/hapcallerSub_{species}_{sample}_{sub}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/hapcallerSub_{wildcards.species}_{wildcards.sample}_{wildcards.sub}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/2_callvars.sh
		fi
		cp {input} $temp_folder
		cd $temp_folder
		ploidy_file=$(awk -F/ '{{print $NF}}' <<< {input.ploidy_input})
		s_ploidy=$(awk '{{if ($1==\"{wildcards.sample}\") print $2 }}' $ploidy_file )
		echo "Sample ploidy is "$s_ploidy >> {log}
		sub_interval_list=$(awk -F/ '{{print $NF}}' <<< {input.sub_interval_list})
		sub_interval=$(awk '{{if ($1==\"{wildcards.sub}\") print $2 }}' $sub_interval_list )

		$GATK4 HaplotypeCaller -I {wildcards.species}_{wildcards.sample}.merged.dedup.bam -R {wildcards.species}.fasta -O {wildcards.species}_{wildcards.sample}_{wildcards.sub}.g.vcf.gz -ERC GVCF --min-base-quality-score {config[hapcaller_minbaseq]} --minimum-mapping-quality {config[hapcaller_minmapq]} -ploidy $s_ploidy -stand-call-conf 20 --pcr-indel-model NONE --max-genotype-count 350 -L $sub_interval &>> {log}
		#removed rf BadMate
		cp {wildcards.species}_{wildcards.sample}_{wildcards.sub}.g.vcf.gz {output.gvcf_out} 		
		"""


#def get_sample_gvcfs(wildcards):
#	return [config["gvcf_dir"]+"/"+wildcards.species+"_"+s+"_"wildcards.sub".gvcf.gz" for s in sample_dict]

def GenomicsDBimportSub_mem_mb(wildcards, attempt):
	return int(config["GenomicsDBimportSub_mem_mb"]+(config["GenomicsDBimportSub_mem_mb"]*(attempt-1)*config["repeat_mem_mb_factor"]))
def GenomicsDBimportSub_disk_mb(wildcards, attempt):
	return int(config["GenomicsDBimportSub_disk_mb"]+(config["GenomicsDBimportSub_disk_mb"]*(attempt-1)*config["repeat_disk_mb_factor"]))
def GenomicsDBimportSub_runtime(wildcards, attempt):
	GenomicsDBImportSub_runtime_seconds=parse_timespan(config["GenomicsDBimportSub_runtime"])
	return str(GenomicsDBImportSub_runtime_seconds+int((GenomicsDBImportSub_runtime_seconds*(attempt-1))*config["repeat_runtime_factor"]))+"s"
#GenomicsDBimportSub_runtime_cats=config["GenomicsDBimportSub_runtime"].split(":")
#	return str(int(GenomicsDBimportSub_runtime_cats[0])+int(int(GenomicsDBimportSub_runtime_cats[0])*(attempt-1)*config["repeat_runtime_factor"]))+":"+GenomicsDBimportSub_runtime_cats[1]+":"+GenomicsDBimportSub_runtime_cats[2]

rule GenomicsDBimportSub:
	input:
		gvfs=expand(config["gvcf_dir"]+"/{{species}}_{sample}_{{sub}}.gvcf.gz", sample=list(sample_dict.keys())),
		sub_interval_list=config["sub_intervals"]
	output:
		directory(config["gvcf_dir"]+"/{species}_{sub}_GenomicsDB")
	threads: 4
	resources:
		mem_mb=GenomicsDBimportSub_mem_mb,
		disk_mb=GenomicsDBimportSub_disk_mb,
		runtime=GenomicsDBimportSub_runtime
	log:
		config["log_dir"]+"/GenomicsDBimportSub_{species}_{sub}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/GenomicsDBimportSub_{wildcards.species}_{wildcards.sub}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/2_callvars.sh
		fi
		cp {input} $temp_folder
		cd $temp_folder
		mkdir tmp
		cp {config[sub_intervals]} .
		sub_interval_list=$(awk -F/ '{{print $NF}}' <<< {input.sub_interval_list})
		sub_interval=$(awk '{{if ($1==\"{wildcards.sub}\") print $2 }}' $sub_interval_list )
		touch cohort.sample_map
		for in_gvzf in *.gvcf.gz
		do
			$GATK4 IndexFeatureFile -I $in_gvzf &>> {log}
			echo $in_gvzf | awk -F_ '{{print $2\"\\t\"$0}}' >> cohort.sample_map 
		done
		
		$GATK4 GenomicsDBImport --genomicsdb-workspace-path GDB_database --batch-size 50 -L $sub_interval --sample-name-map cohort.sample_map --tmp-dir tmp --reader-threads 4 &>> {log}
		cp -rf GDB_database {output}
		"""

#def get_intervals(wildcards):
#	return [config["gvcf_dir"]+"/"+wildcards.species+"_"+u+".gvcf.gz" for u in interval_list]

#rule collect_intervals:
#	input:	
#		expand(config["gvcf_dir"]+"/{{species}}_{sub}_GenomicsDB", sub=interval_list)
#	output:
#		"GenomicsDBintervals_{species}.dummy"
#	shell:
#		"""
#		touch {output}
#		"""
