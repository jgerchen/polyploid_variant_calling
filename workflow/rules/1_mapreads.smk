import glob
import re
from humanfriendly import parse_timespan

configfile: "../config/1_mapreads.yaml"

#fastqname\tnewname\tadapter
#files will be named: newname_1.bam etc. and eventually newname.merged.bam
sample_dict={}
with open(config["sample_list"]) as sample_list:
	for line in sample_list:
		sample_cats=line.strip().split("\t")
		sample_libs=[]
		for s_cat_multiple in sample_cats[1].split(","):
			if type(config["fastq_dir"])==str:
				sample_glob=glob.glob(config["fastq_dir"]+"/**/*"+s_cat_multiple+"*.f*q*", recursive=True)
				if len(sample_glob)>0:
					sample_libs=list(set(sample_libs).union(set(sample_glob)))
			elif type(config["fastq_dir"]) in {list, tuple}:
				for fastq_dir  in config["fastq_dir"]:
					sample_glob=glob.glob(fastq_dir+"/**/*"+s_cat_multiple+"*.f*q*", recursive=True)
					if len(sample_glob)>0:
						sample_libs=list(set(sample_libs).union(set(sample_glob)))
		assert len(sample_libs)>=2, "There have to be at least 2 libraries per sample, however sample %s only has %s." % (sample_cats[0], len(sample_libs))
		r1_re=re.compile(".*R1.*")
		samples_r1=list(filter(r1_re.match, sample_libs))
		assert len(samples_r1)>0, "No library from sample %s matches regular expression for first read pair" % sample_cats[0]
		sample_libs_dict={str(i+1):(samples_r1[i] ,samples_r1[i].replace("R1", "R2")) for i in range(len(samples_r1))}
		sample_dict.update({sample_cats[0]:(sample_libs_dict, config["adapter_dir"]+"/"+sample_cats[2])})	
def get_fwd_reads(wildcards):
	return sample_dict[wildcards.sample][0][wildcards.lib][0]
def get_rev_reads(wildcards):
	return sample_dict[wildcards.sample][0][wildcards.lib][1]
def get_adapter_file(wildcards):
	return sample_dict[wildcards.sample][1]

##rule fastqc:
#	input:
#		fwd_reads=get_fwd_reads,
#		rev_reads=get_rev_reads
#	output:
#		out_fwd=directory(config["fastqc_dir"]+"/{sample}_{lib}")
#	threads: 2
#	resources:
#		mem_mb=8000,
#		disk_mb=4000,
#		runtime="4:00:00"
#	log:
#		config["log_dir"]+"/fastqc_{sample}_{lib}.log"
#	shell:
#		"""
#		temp_folder={config[temp_dir]}/map_reads_{wildcards.species}_{wildcards.sample}
#		mkdir -p $temp_folder
#		trap 'rm -rf $temp_folder' TERM EXIT
#		if [ {config[load_cluster_code]} -eq 1 ]
#		then
#			source {config[cluster_code_dir]}/1_mapreads.sh
#		fi
#		cp {input} $temp_folder
#		cd $temp_folder
#		mkdir fastqc_out
#		fastq_files=$(ls *.f*q*)
#		fastqc -o fastqc_out -t 2 -f fastq  $fastq_files
#		cp -r fastqc_files {output}
#		"""

def trimmomatic_mem_mb(wildcards, attempt):
	return int(config["trimmomatic_mem_mb"]+(config["trimmomatic_mem_mb"]*(attempt-1)*config["repeat_mem_mb_factor"]))
def trimmomatic_disk_mb(wildcards, attempt):
	return int(config["trimmomatic_disk_mb"]+(config["trimmomatic_disk_mb"]*(attempt-1)*config["repeat_disk_mb_factor"]))
def trimmomatic_runtime(wildcards, attempt):
	trimmomatic_runtime_seconds=parse_timespan(config["trimmomatic_runtime"])
	return str(trimmomatic_runtime_seconds+int((trimmomatic_runtime_seconds*(attempt-1))*config["repeat_runtime_factor"]))+"s"
	#trimmomatic_runtime_cats=config["trimmomatic_runtime"].split(":")
	#return str(int(trimmomatic_runtime_cats[0])+int(int(trimmomatic_runtime_cats[0])*(attempt-1)*config["repeat_runtime_factor"]))+":"+trimmomatic_runtime_cats[1]+":"+trimmomatic_runtime_cats[2]

rule trimmomatic:
	input:
		fwd_reads=get_fwd_reads,
		rev_reads=get_rev_reads,
		adapter_file=get_adapter_file
	output:
		fwd_reads_trimmed=config["fastq_trimmed_dir"]+"/{sample}_{lib}_trimmed_R1.fastq.gz",
		fwd_reads_unpaired=config["fastq_trimmed_dir"]+"/{sample}_{lib}_trimmed_U1.fastq.gz",
		rev_reads_trimmed=config["fastq_trimmed_dir"]+"/{sample}_{lib}_trimmed_R2.fastq.gz",
		rev_reads_unpaired=config["fastq_trimmed_dir"]+"/{sample}_{lib}_trimmed_U2.fastq.gz",
		pre_fwd=report(config["report_dir"]+"/trimmomatic/fastQC_{sample}_{lib}/pre_fwd.html", category="trimmomatic", subcategory="fastQC before trimming", labels={"sample":"{sample}", "library":"{lib}", "read-pair":"R1"}) if config["run_fastqc"]==1 else [],
		pre_rev=report(config["report_dir"]+"/trimmomatic/fastQC_{sample}_{lib}/pre_rev.html", category="trimmomatic", subcategory="fastQC before trimming", labels={"sample":"{sample}", "library":"{lib}", "read-pair":"R2"}) if config["run_fastqc"]==1 else [],
		post_fwd_paired=report(config["report_dir"]+"/trimmomatic/fastQC_{sample}_{lib}/post_fwd_paired.html", category="trimmomatic", subcategory="fastQC after trimming", labels={"sample":"{sample}", "library":"{lib}", "read-pair":"R1", "paired":"Yes"}) if config["run_fastqc"]==True else [],
		post_rev_paired=report(config["report_dir"]+"/trimmomatic/fastQC_{sample}_{lib}/post_rev_paired.html", category="trimmomatic", subcategory="fastQC after trimming", labels={"sample":"{sample}", "library":"{lib}", "read-pair":"R2", "paired":"Yes"}) if config["run_fastqc"]==True else [],
		post_fwd_unpaired=report(config["report_dir"]+"/trimmomatic/fastQC_{sample}_{lib}/post_fwd_unpaired.html", category="trimmomatic", subcategory="fastQC after trimming", labels={"sample":"{sample}", "library":"{lib}", "read-pair":"R1", "paired":"Singleton"}) if config["run_fastqc"]==True else [],
		post_rev_unpaired=report(config["report_dir"]+"/trimmomatic/fastQC_{sample}_{lib}/post_rev_unpaired.html", category="trimmomatic", subcategory="fastQC after trimming", labels={"sample":"{sample}", "library":"{lib}", "read-pair":"R2", "paired":"Singleton"}) if config["run_fastqc"]==True else []

	threads: 4
	resources:
		mem_mb=trimmomatic_mem_mb,
		disk_mb=trimmomatic_disk_mb,
		runtime=trimmomatic_runtime
	log:
		config["log_dir"]+"/trimmomatic_{sample}_{lib}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/trimmomatic_{wildcards.sample}_{wildcards.lib}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/1_mapreads.sh
		fi
		cp {input} $temp_folder
		cd $temp_folder
		fwd_reads=$(awk -F/ '{{print $NF}}' <<< {input.fwd_reads})
		rev_reads=$(awk -F/ '{{print $NF}}' <<< {input.rev_reads})
		if [ {config[run_fastqc]} -eq 1 ]
		then
			mkdir fastqc_pre_out
			fastqc -o fastqc_pre_out -t {threads} -f fastq  *.fastq.gz &>>{log}
			cp fastqc_pre_out/*R1*.html {output.pre_fwd}
			cp fastqc_pre_out/*R2*.html {output.pre_rev}
		fi

		adapter_file=$(awk -F/ '{{print $NF}}' <<< {input.adapter_file} )
		trimmomatic PE -threads {threads} -trimlog trimmomatic.log $fwd_reads $rev_reads output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:$adapter_file:2:23:10 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:50 &>> {log}
		if [ {config[run_fastqc]} -eq 1 ]
		then
			mkdir fastqc_post_out
			fastqc -o fastqc_post_out -t {threads} -f fastq  output_{{forward,reverse}}_{{paired,unpaired}}.fq.gz
			cp fastqc_post_out/*forward_paired*.html {output.post_fwd_paired}
			cp fastqc_post_out/*reverse_paired*.html {output.post_rev_paired}
			cp fastqc_post_out/*forward_unpaired*.html {output.post_fwd_unpaired}
			cp fastqc_post_out/*reverse_unpaired*.html {output.post_rev_unpaired}
		fi
		#cat trimmomatic.log >> {log}
		cp output_forward_paired.fq.gz {output.fwd_reads_trimmed}	
		cp output_forward_unpaired.fq.gz {output.fwd_reads_unpaired}	
		cp output_reverse_paired.fq.gz {output.rev_reads_trimmed}	
		cp output_reverse_unpaired.fq.gz {output.rev_reads_unpaired}	
		"""

def map_reads_mem_mb(wildcards, attempt):
	return int(config["map_reads_mem_mb"]+(config["map_reads_mem_mb"]*(attempt-1)*config["repeat_mem_mb_factor"]))
def map_reads_disk_mb(wildcards, attempt):
	return int(config["map_reads_disk_mb"]+(config["map_reads_disk_mb"]*(attempt-1)*config["repeat_disk_mb_factor"]))
def map_reads_runtime(wildcards, attempt):
	map_reads_runtime_seconds=parse_timespan(config["map_reads_runtime"])
	return str(map_reads_runtime_seconds+int((map_reads_runtime_seconds*(attempt-1))*config["repeat_runtime_factor"]))+"s"
	#	map_reads_runtime_cats=config["map_reads_runtime"].split(":")
	#return str(int(map_reads_runtime_cats[0])+int(int(map_reads_runtime_cats[0])*(attempt-1)*config["repeat_runtime_factor"]))+":"+map_reads_runtime_cats[1]+":"+map_reads_runtime_cats[2]
rule map_reads:
	input:
		reference=config["fasta_dir"]+"/{species}.fasta",
		reference_amb=config["fasta_dir"]+"/{species}.fasta.amb",
		reference_ann=config["fasta_dir"]+"/{species}.fasta.ann",
		reference_bwt=config["fasta_dir"]+"/{species}.fasta.bwt",
		reference_fai=config["fasta_dir"]+"/{species}.fasta.fai",
		reference_pac=config["fasta_dir"]+"/{species}.fasta.pac",
		reference_sa=config["fasta_dir"]+"/{species}.fasta.sa",
		fwd_reads_trimmed=config["fastq_trimmed_dir"]+"/{sample}_{lib}_trimmed_R1.fastq.gz" if config["trim_reads"]==True else get_fwd_reads,
		fwd_reads_unpaired=config["fastq_trimmed_dir"]+"/{sample}_{lib}_trimmed_U1.fastq.gz" if config["trim_reads"]==True else [],
		rev_reads_trimmed=config["fastq_trimmed_dir"]+"/{sample}_{lib}_trimmed_R2.fastq.gz" if config["trim_reads"]==True else get_rev_reads,
		rev_reads_unpaired=config["fastq_trimmed_dir"]+"/{sample}_{lib}_trimmed_U2.fastq.gz" if config["trim_reads"]==True else []

	output:
		bam_out=config["bam_dir"]+"/{species}_{sample}_{lib}.bam",
		bam_index=config["bam_dir"]+"/{species}_{sample}_{lib}.bam.bai"
	threads: 4
	resources:
		mem_mb=map_reads_mem_mb,
		disk_mb=map_reads_disk_mb,
		runtime=map_reads_runtime
	log:
		config["log_dir"]+"/map_reads_{species}_{sample}_{lib}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/map_reads_{wildcards.species}_{wildcards.sample}_{wildcards.lib}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/1_mapreads.sh
		fi
		cp {input} $temp_folder
		cd $temp_folder
		fwd_reads=$(awk -F/ '{{print $NF}}' <<< {input.fwd_reads_trimmed})
		rev_reads=$(awk -F/ '{{print $NF}}' <<< {input.rev_reads_trimmed})
		fwd_reads_unp=$(awk -F/ '{{print $NF}}' <<< {input.fwd_reads_unpaired})
		rev_reads_unp=$(awk -F/ '{{print $NF}}' <<< {input.rev_reads_unpaired})
		bwa mem -t {threads} -R '@RG\\tID:{wildcards.species}_{wildcards.sample}_{wildcards.lib}\\tLB:{wildcards.species}_{wildcards.sample}_{wildcards.lib}\\tSM:{wildcards.species}_{wildcards.sample}\\tPL:illumina' {wildcards.species}.fasta  $fwd_reads $rev_reads | samtools sort - | samtools view -bh -o aligned_pe.bam &>> {log}
		if [ {config[trim_reads]} -eq 1 ]
		then
			bwa mem -t {threads} -R '@RG\\tID:{wildcards.species}_{wildcards.sample}_{wildcards.lib}\\tLB:{wildcards.species}_{wildcards.sample}_{wildcards.lib}\\tSM:{wildcards.species}_{wildcards.sample}\\tPL:illumina' {wildcards.species}.fasta $fwd_reads_unp | samtools sort - | samtools view -bh -o aligned_fwd_unp.bam &>> {log}
			bwa mem -t {threads} -R '@RG\\tID:{wildcards.species}_{wildcards.sample}_{wildcards.lib}\\tLB:{wildcards.species}_{wildcards.sample}_{wildcards.lib}\\tSM:{wildcards.species}_{wildcards.sample}\\tPL:illumina' {wildcards.species}.fasta $rev_reads_unp | samtools sort - | samtools view -bh -o aligned_rev_unp.bam &>> {log}
			samtools merge aligned_merged.bam aligned_pe.bam aligned_fwd_unp.bam aligned_rev_unp.bam
		else
			mv aligned_pe.bam aligned_merged.bam
		fi
		samtools index aligned_merged.bam
		mkdir -p {config[bam_dir]}
		cp aligned_merged.bam {output.bam_out}
		cp aligned_merged.bam.bai {output.bam_index}
		"""

def get_bam_libs(wildcards):
	return [config["bam_dir"]+"/"+wildcards.species+"_"+wildcards.sample+"_"+i+".bam" for i in list(sample_dict[wildcards.sample][0].keys())]

def merge_bams_mem_mb(wildcards, attempt):
	return int(config["merge_bams_mem_mb"]+(config["merge_bams_mem_mb"]*(attempt-1)*config["repeat_mem_mb_factor"]))
def merge_bams_disk_mb(wildcards, attempt):
	return int(config["merge_bams_disk_mb"]+(config["merge_bams_disk_mb"]*(attempt-1)*config["repeat_disk_mb_factor"]))
def merge_bams_runtime(wildcards, attempt):
	merge_bams_runtime_seconds=parse_timespan(config["merge_bams_runtime"])
	return str(merge_bams_runtime_seconds+int((merge_bams_runtime_seconds*(attempt-1))*config["repeat_runtime_factor"]))+"s"
#merge_bams_runtime_cats=config["merge_bams_runtime"].split(":")
#	return str(int(merge_bams_runtime_cats[0])+int(int(merge_bams_runtime_cats[0])*(attempt-1)*config["repeat_runtime_factor"]))+":"+merge_bams_runtime_cats[1]+":"+merge_bams_runtime_cats[2]
rule merge_bams_deduplicate:
	input:
		get_bam_libs
	output:
		merged=config["bam_dir"]+"/{species}_{sample}.merged.bam",
		merged_index=config["bam_dir"]+"/{species}_{sample}.merged.bam.bai",
		merged_dedup=config["bam_dir"]+"/{species}_{sample}.merged.dedup.bam",
		merged_dedup_index=config["bam_dir"]+"/{species}_{sample}.merged.dedup.bam.bai",
		stats_flagstat=config["report_dir"]+"/merge_bams_dedup/{species}_{sample}.flagstats.tsv",
		#stats_flagstat_table=report(config["report_dir"]+"/merge_bams_dedup/table_flagstats/{species}_{sample}.rst", category="Merge bams and deduplicate", subcategory="Flagstat table", labels={"Sample":"{sample}"}),
		stats_flagstat_plot=report(config["report_dir"]+"/merge_bams_dedup/plot_flagstats/{species}_{sample}.flagstats.pdf", category="Merge bams and deduplicate", subcategory="Flagstat plot", labels={"Sample":"{sample}"}),
		stats_stat=config["report_dir"]+"/merge_bams_dedup/{species}_{sample}.stats",
		stats_plot_acgt_cycles=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_agct_cycles.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"ACGT content per cycle"}),
		stats_plot_coverage=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_coverage.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Coverage plot"}),
		stats_plot_gc_content=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_gc_content.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"GC content"}),
		stats_plot_gc_depth=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_gc_depth.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Mapped depth vs. GC"}),
		stats_plot_indel_cycles=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_indel_cycles.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Indels per cycle"}),
		stats_plot_indel_dist=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_indel_dist.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Indel lengths"}),
		stats_plot_insert_size=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_insert_size.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Insert size"}),
		stats_plot_quals=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_quals.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Quality per cycle"}),
		stats_plot_quals_hm=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_quals_hm.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Quality per cycle 2D plot"}),
		stats_plot_quals2=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_quals2.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Quality per cycle Plot 2"}),
		stats_plot_quals3=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_quals3.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Quality per cycle Plot 3"}),
		
	threads: 2
	resources:
		mem_mb=merge_bams_mem_mb,
		disk_mb=merge_bams_disk_mb,
		runtime=merge_bams_runtime
	log:
		config["log_dir"]+"/deduplicate_{species}_{sample}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/merge_bams_{wildcards.species}_{wildcards.sample}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/1_mapreads.sh
		fi
		cp {input} $temp_folder
		cp scripts/plot_flagstats.R $temp_folder
		cd $temp_folder
		samtools merge all_merged.bam *.bam &>> {log}
		samtools index all_merged.bam &>> {log}
		limit=$(echo `ulimit -n` - 50 | bc)
		java -jar -XX:ParallelGCThreads=2 -Xmx12g $PICARD MarkDuplicates I=all_merged.bam O=all_merged.dedup.bam MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=$limit M=dup_metrics.log ASSUME_SORTED=true TAGGING_POLICY=All &>>{log}
		samtools index all_merged.dedup.bam &>>{log}
		cat dup_metrics.log >> {log}
		samtools flagstat all_merged.dedup.bam -O tsv > all_merged.flagstat.tsv
		mkdir -p {config[report_dir]}/merge_bams_dedup 
		cp all_merged.flagstat.tsv {output.stats_flagstat}
		Rscript plot_flagstats.R all_merged.flagstat.tsv {wildcards.sample} {wildcards.species}_{wildcards.sample}.flagstats.pdf
		cp {wildcards.species}_{wildcards.sample}.flagstats.pdf {output.stats_flagstat_plot}
		samtools stats all_merged.dedup.bam > all_merged.stats
		mkdir plot_stats
		plot-bamstats -p plot_stats/{wildcards.species}_{wildcards.sample} all_merged.stats
		cp all_merged.stats {output.stats_stat}
		mkdir -p {config[report_dir]}/merge_bams_dedup/plot_bamstats 
		cp plot_stats/{wildcards.species}_{wildcards.sample}-acgt-cycles.png {output.stats_plot_acgt_cycles}
		cp plot_stats/{wildcards.species}_{wildcards.sample}-coverage.png {output.stats_plot_coverage}
		cp plot_stats/{wildcards.species}_{wildcards.sample}-gc-content.png {output.stats_plot_gc_content}
		cp plot_stats/{wildcards.species}_{wildcards.sample}-gc-depth.png {output.stats_plot_gc_depth}
		cp plot_stats/{wildcards.species}_{wildcards.sample}-indel-cycles.png {output.stats_plot_indel_cycles}
		cp plot_stats/{wildcards.species}_{wildcards.sample}-indel-dist.png {output.stats_plot_indel_dist}
		cp plot_stats/{wildcards.species}_{wildcards.sample}-insert-size.png {output.stats_plot_insert_size}
		cp plot_stats/{wildcards.species}_{wildcards.sample}-quals.png {output.stats_plot_quals}
		cp plot_stats/{wildcards.species}_{wildcards.sample}-quals-hm.png {output.stats_plot_quals_hm}
		cp plot_stats/{wildcards.species}_{wildcards.sample}-quals2.png {output.stats_plot_quals2}
		cp plot_stats/{wildcards.species}_{wildcards.sample}-quals3.png {output.stats_plot_quals3}
		cp all_merged.bam {output.merged}
		cp all_merged.bam.bai {output.merged_index}
		cp all_merged.dedup.bam {output.merged_dedup}
		cp all_merged.dedup.bam.bai {output.merged_dedup_index}
		"""
