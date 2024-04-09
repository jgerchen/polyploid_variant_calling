#import glob
#import re
configfile: workflow.source_path("../../config/1_mapreads.yaml")
import os.path
from humanfriendly import parse_timespan

sample_dict={}
if os.path.isfile(config["sample_read_file"]):
	with open(config["sample_read_file"]) as s_read_file:
		for s_read_line in s_read_file:
			s_read_line_cats=s_read_line.strip().split()
			if s_read_line_cats[0] not in sample_dict:
				sample_dict.update({s_read_line_cats[0]:({s_read_line_cats[1]:(s_read_line_cats[2], s_read_line_cats[3])},s_read_line_cats[4],s_read_line_cats[5])})
			else:
				if s_read_line_cats[1] not in sample_dict[s_read_line_cats[0]][0]:
					sample_dict[s_read_line_cats[0]][0].update({s_read_line_cats[1]:(s_read_line_cats[2], s_read_line_cats[3])})

checkpoint get_sample_reads:
	input:
		sample_list=config["sample_list"]
	output:
		sample_read_file=config["sample_read_file"]
	threads: 1
	resources:
		mem_mb=config["get_sample_reads_mem_mb"],
		disk_mb=config["get_sample_reads_disk_mb"],
		runtime=config["get_sample_reads_runtime"]
	log:
		config["log_dir"]+"/get_sample_reads.log"
	run:
		import glob
		import re
		pe_strings=[read_i.split(":") for read_i in config["paired_end_res"].split(",")]

		with open(input.sample_list) as sample_list:
			with open(output.sample_read_file, "w") as s_read_out:
				for line in sample_list:
					if len(line.strip())>0:
						sample_cats=line.strip().split()
						assert len(sample_cats)==4, "The sample list must have 4 columns (sample name,file name(s), adapter file and ploidy), but here it's %s" % len(sample_cats)
						assert "_" not in sample_cats[0], "Sample names must not contain underscores! Remove the underscore in sample name %s and try again" % sample_cats[0]
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
						#get list of wildcard pairs
						for pe_re in pe_strings:
							#1. match re1
							r1_re=re.compile(pe_re[0])
							samples_r1=sorted(list(filter(r1_re.match, sample_libs)))
							#2. match re2
							r2_re=re.compile(pe_re[1])
							samples_r2=sorted(list(filter(r2_re.match, sample_libs)))
							#ensure the length of both matches is the same
							assert len(samples_r1)==len(samples_r2), "Regular expressions for forward and reverse reads have to match the same number of times, but for sample %s RE %s matches %s times but RE %s matches %s times." % (sample_cats[0], pe_re[0], len(samples_r1), pe_re[1], len(samples_r2))
							for re_pair_i in range(len(samples_r1)):
								s_read_out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (sample_cats[0], re_pair_i, samples_r1[re_pair_i], samples_r2[re_pair_i], config["adapter_dir"]+"/"+sample_cats[2], sample_cats[3]))


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
		pre_fwd_zip=config["report_dir"]+"/trimmomatic/fastQC_{sample}_{lib}/pre_fwd.zip" if config["run_fastqc"]==1 else [],
		pre_rev=report(config["report_dir"]+"/trimmomatic/fastQC_{sample}_{lib}/pre_rev.html", category="trimmomatic", subcategory="fastQC before trimming", labels={"sample":"{sample}", "library":"{lib}", "read-pair":"R2"}) if config["run_fastqc"]==1 else [],
		pre_rev_zip=config["report_dir"]+"/trimmomatic/fastQC_{sample}_{lib}/pre_rev.zip" if config["run_fastqc"]==1 else [],
		post_fwd_paired=report(config["report_dir"]+"/trimmomatic/fastQC_{sample}_{lib}/post_fwd_paired.html", category="trimmomatic", subcategory="fastQC after trimming", labels={"sample":"{sample}", "library":"{lib}", "read-pair":"R1", "paired":"Yes"}) if config["run_fastqc"]==True else [],
		post_fwd_paired_zip=config["report_dir"]+"/trimmomatic/fastQC_{sample}_{lib}/post_fwd_paired.zip" if config["run_fastqc"]==True else [],
		post_rev_paired=report(config["report_dir"]+"/trimmomatic/fastQC_{sample}_{lib}/post_rev_paired.html", category="trimmomatic", subcategory="fastQC after trimming", labels={"sample":"{sample}", "library":"{lib}", "read-pair":"R2", "paired":"Yes"}) if config["run_fastqc"]==True else [],
		post_rev_paired_zip=config["report_dir"]+"/trimmomatic/fastQC_{sample}_{lib}/post_rev_paired.zip" if config["run_fastqc"]==True else [],
		post_fwd_unpaired=report(config["report_dir"]+"/trimmomatic/fastQC_{sample}_{lib}/post_fwd_unpaired.html", category="trimmomatic", subcategory="fastQC after trimming", labels={"sample":"{sample}", "library":"{lib}", "read-pair":"R1", "paired":"Singleton"}) if config["run_fastqc"]==True else [],
		post_fwd_unpaired_zip=config["report_dir"]+"/trimmomatic/fastQC_{sample}_{lib}/post_fwd_unpaired.zip" if config["run_fastqc"]==True else [],
		post_rev_unpaired=report(config["report_dir"]+"/trimmomatic/fastQC_{sample}_{lib}/post_rev_unpaired.html", category="trimmomatic", subcategory="fastQC after trimming", labels={"sample":"{sample}", "library":"{lib}", "read-pair":"R2", "paired":"Singleton"}) if config["run_fastqc"]==True else [],
		post_rev_unpaired_zip=config["report_dir"]+"/trimmomatic/fastQC_{sample}_{lib}/post_rev_unpaired.zip" if config["run_fastqc"]==True else []

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
		mv $fwd_reads input_fwd_reads.fastq.gz
		mv $rev_reads input_rev_reads.fastq.gz

		if [ {config[run_fastqc]} -eq 1 ]
		then
			mkdir fastqc_pre_out
			fastqc -o fastqc_pre_out -t {threads} -f fastq  *.fastq.gz &>>{log}
			cp fastqc_pre_out/*input_fwd_reads*.html {output.pre_fwd}
			cp fastqc_pre_out/*input_fwd_reads*.zip {output.pre_fwd_zip}
			cp fastqc_pre_out/*input_rev_reads*.html {output.pre_rev}
			cp fastqc_pre_out/*input_rev_reads*.zip {output.pre_rev_zip}
		fi

		adapter_file=$(awk -F/ '{{print $NF}}' <<< {input.adapter_file} )
		trimmomatic PE -threads {threads} -trimlog trimmomatic.log input_fwd_reads.fastq.gz input_rev_reads.fastq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz HEADCROP:{config[trimmomatic_headcrop]} ILLUMINACLIP:$adapter_file:{config[trimmomatic_clip_seedmismatches]}:{config[trimmomatic_clip_palindromethresh]}:{config[trimmomatic_clip_simplethresh]} TRAILING:{config[trimmomatic_trailing]} LEADING:{config[trimmomatic_leading]}  SLIDINGWINDOW:{config[trimmomatic_slidingw_size]}:{config[trimmomatic_slidingw_qual]} MINLEN:{config[trimmomatic_minlen]} &>> {log}
		if [ {config[run_fastqc]} -eq 1 ]
		then
			mkdir fastqc_post_out
			fastqc -o fastqc_post_out -t {threads} -f fastq  output_{{forward,reverse}}_{{paired,unpaired}}.fq.gz
			cp fastqc_post_out/*forward_paired*.html {output.post_fwd_paired}
			cp fastqc_post_out/*forward_paired*.zip {output.post_fwd_paired_zip}
			cp fastqc_post_out/*reverse_paired*.html {output.post_rev_paired}
			cp fastqc_post_out/*reverse_paired*.zip {output.post_rev_paired_zip}
			cp fastqc_post_out/*forward_unpaired*.html {output.post_fwd_unpaired}
			cp fastqc_post_out/*forward_unpaired*.zip {output.post_fwd_unpaired_zip}
			cp fastqc_post_out/*reverse_unpaired*.html {output.post_rev_unpaired}
			cp fastqc_post_out/*reverse_unpaired*.zip {output.post_rev_unpaired_zip}
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
		bwa mem -t {threads} -R '@RG\\tID:{wildcards.species}_{wildcards.sample}_{wildcards.lib}\\tLB:{wildcards.species}_{wildcards.sample}_{wildcards.lib}\\tSM:{wildcards.species}_{wildcards.sample}\\tPL:illumina' {wildcards.species}.fasta  $fwd_reads $rev_reads | samtools sort - | samtools view -bh -o aligned_pe.bam &>> {log}
		if [ {config[trim_reads]} -eq 1 ]
		then
			fwd_reads_unp=$(awk -F/ '{{print $NF}}' <<< {input.fwd_reads_unpaired})
			rev_reads_unp=$(awk -F/ '{{print $NF}}' <<< {input.rev_reads_unpaired})
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
		stats_plot_acgt_cycles=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_agct_cycles.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"ACGT content per cycle"}) if config["plot_bamstats"]==True else [],
		stats_plot_coverage=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_coverage.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Coverage plot"}) if config["plot_bamstats"]==True else [],
		stats_plot_gc_content=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_gc_content.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"GC content"}) if config["plot_bamstats"]==True else [],
		stats_plot_gc_depth=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_gc_depth.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Mapped depth vs. GC"}) if config["plot_bamstats"]==True else [],
		stats_plot_indel_cycles=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_indel_cycles.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Indels per cycle"}) if config["plot_bamstats"]==True else [],
		stats_plot_indel_dist=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_indel_dist.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Indel lengths"})  if config["plot_bamstats"]==True else [],
		stats_plot_insert_size=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_insert_size.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Insert size"})  if config["plot_bamstats"]==True else [],
		stats_plot_quals=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_quals.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Quality per cycle"})  if config["plot_bamstats"]==True else [],
		stats_plot_quals_hm=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_quals_hm.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Quality per cycle 2D plot"})  if config["plot_bamstats"]==True else [],
		stats_plot_quals2=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_quals2.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Quality per cycle Plot 2"}) if config["plot_bamstats"]==True else [],
		stats_plot_quals3=report(config["report_dir"]+"/merge_bams_dedup/plot_bamstats/{species}_{sample}_quals3.png", category="Merge bams and deduplicate", subcategory="Plot bamstats", labels={"Sample":"{sample}", "Plot":"Quality per cycle Plot 3"})  if config["plot_bamstats"]==True else [],
		
	threads: 1
	resources:
		mem_mb=merge_bams_mem_mb,
		disk_mb=merge_bams_disk_mb,
		runtime=merge_bams_runtime
	params:
		plot_flagstats_script=workflow.source_path("../scripts/plot_flagstats.R")
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
		cp {params.plot_flagstats_script} $temp_folder
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
		cp all_merged.stats {output.stats_stat}
		if [ {config[plot_bamstats]} -eq 1 ]
		then

			mkdir plot_stats
			plot-bamstats -p plot_stats/{wildcards.species}_{wildcards.sample} all_merged.stats
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
		fi	
		cp all_merged.bam {output.merged}
		cp all_merged.bam.bai {output.merged_index}
		cp all_merged.dedup.bam {output.merged_dedup}
		cp all_merged.dedup.bam.bai {output.merged_dedup_index}
		"""

def bam_depth_mem_mb(wildcards, attempt):
	return int(config["bam_depth_mem_mb"]+(config["bam_depth_mem_mb"]*(attempt-1)*config["repeat_mem_mb_factor"]))
def bam_depth_disk_mb(wildcards, attempt):
	return int(config["bam_depth_disk_mb"]+(config["bam_depth_disk_mb"]*(attempt-1)*config["repeat_disk_mb_factor"]))
def bam_depth_runtime(wildcards, attempt):
	bam_depth_runtime_seconds=parse_timespan(config["bam_depth_runtime"])
	return str(bam_depth_runtime_seconds+int((bam_depth_runtime_seconds*(attempt-1))*config["repeat_runtime_factor"]))+"s"


rule bam_depth:
	input:
		merged_dedup=config["bam_dir"]+"/{species}_{sample}.merged.dedup.bam",
		merged_dedup_index=config["bam_dir"]+"/{species}_{sample}.merged.dedup.bam.bai"
	output:
		stats_depth_contig=config["report_dir"]+"/bam_depth/{species}_{sample}.chr.stat.gz",
		stats_depth_contig_clean=config["report_dir"]+"/bam_depth/{species}_{sample}.chr.clean.stat.gz",
		stats_depth_window=config["report_dir"]+"/bam_depth/{species}_{sample}.win.stat.gz"
	threads: 6
	resources:
		mem_mb=bam_depth_mem_mb,
		disk_mb=bam_depth_disk_mb,
		runtime=bam_depth_runtime
	log:
		config["log_dir"]+"/bam_depth_{species}_{sample}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/bam_depth_{wildcards.species}_{wildcards.sample}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/1_mapreads.sh
		fi
		cp {input} $temp_folder
		cd $temp_folder
		$PANDEPTH -i {wildcards.species}_{wildcards.sample}.merged.dedup.bam -o {wildcards.species}_{wildcards.sample} -t {threads} -x 4
		cp {wildcards.species}_{wildcards.sample}.chr.stat.gz {output.stats_depth_contig}
		$PANDEPTH -i {wildcards.species}_{wildcards.sample}.merged.dedup.bam -o {wildcards.species}_{wildcards.sample}.clean -q {config[pandepth_min_q]} -t {threads}
		cp {wildcards.species}_{wildcards.sample}.clean.chr.stat.gz {output.stats_depth_contig_clean}
		$PANDEPTH -i {wildcards.species}_{wildcards.sample}.merged.dedup.bam -w {config[pandepth_window_size]} -o {wildcards.species}_{wildcards.sample} -q {config[pandepth_min_q]} -t {threads}
		cp {wildcards.species}_{wildcards.sample}.win.stat.gz {output.stats_depth_window}
		"""


def get_samples_bam_stats(wildcards):
	#bam_stats_output_dict={"stats_depth_contig":[], "stats_depth_contig_clean":[], "stats_flagstat":[]}
	bam_stats_output_dict={"stats_depth_contig":[], "stats_depth_contig_clean":[], "stats_flagstat":[]}
	with checkpoints.get_sample_reads.get(species=wildcards.species).output[0].open() as f:
		#TODO: check if dictionary is already populated, then there's no need to read the file again!
		if len(sample_dict)==0:
			for sample_line in f:
				samp_cats=sample_line.strip().split()
				if samp_cats[0] not in sample_dict:
					sample_dict.update({samp_cats[0]:({samp_cats[1]:(samp_cats[2], samp_cats[3])},samp_cats[4],samp_cats[5])})
				else:
					if samp_cats[1] not in sample_dict[samp_cats[0]][0]:
						sample_dict[samp_cats[0]][0].update({samp_cats[1]:(samp_cats[2], samp_cats[3])})
		for dict_sample in sample_dict:
			bam_stats_output_dict["stats_depth_contig"].append(config["report_dir"]+"/bam_depth/{species}_"+dict_sample+".chr.stat.gz")
			bam_stats_output_dict["stats_depth_contig_clean"].append(config["report_dir"]+"/bam_depth/{species}_"+dict_sample+".chr.clean.stat.gz")
			bam_stats_output_dict["stats_flagstat"].append(config["report_dir"]+"/merge_bams_dedup/{species}_"+dict_sample+".flagstats.tsv")
	return bam_stats_output_dict


rule merge_bam_stats:
	input:
		unpack(get_samples_bam_stats)
		#stats_depth_contig=expand(config["report_dir"]+"/bam_depth/{{species}}_{sample}.chr.stat.gz", sample=list(sample_dict.keys())),
		#stats_depth_contig_clean=expand(config["report_dir"]+"/bam_depth/{{species}}_{sample}.chr.clean.stat.gz", sample=list(sample_dict.keys())),
		#stats_flagstat=expand(config["report_dir"]+"/merge_bams_dedup/{{species}}_{sample}.flagstats.tsv", sample=list(sample_dict.keys()))
	output:
		stats_table=config["report_dir"]+"/merge_bam_stats/{species}.bam.table.tsv",
		plot_depth=report(config["report_dir"]+"/merge_bam_stats/{species}.depths.boxplot.pdf", category="Merge bam stats", subcategory="Depth plot", labels={"Species":"{species}"}),
		plot_flagstats=report(config["report_dir"]+"/merge_bam_stats/{species}.flagstats.boxplot.pdf", category="Merge bam stats", subcategory="Flagstats plot", labels={"Species":"{species}"})
	threads: 1
	resources:
		mem_mb=16000,
		disk_mb=4000,
		runtime="3600s"
	params:
		collect_bamstats_script=workflow.source_path("../scripts/collect_bam_stats.py"),
		plot_merged_bamstats_script=workflow.source_path("../scripts/plot_merge_bam_stats.R")
	log:
		config["log_dir"]+"/merge_bam_stats_{species}.log"
	shell:
		"""
		temp_folder={config[temp_dir]}/merge_bam_stats_{wildcards.species}
		mkdir -p $temp_folder
		trap 'rm -rf $temp_folder' TERM EXIT
		if [ {config[load_cluster_code]} -eq 1 ]
		then
			source {config[cluster_code_dir]}/1_mapreads.sh
		fi
		cp {input} $temp_folder
		cp {params.collect_bamstats_script} $temp_folder
		cp {params.plot_merged_bamstats_script} $temp_folder
		cd $temp_folder
		python3 collect_bam_stats.py --depth *.chr.stat.gz --clean_depth *.chr.clean.stat.gz --flagstats *.flagstats.tsv -o {wildcards.species}.bam.table.tsv
		cp {wildcards.species}.bam.table.tsv {output.stats_table}
		Rscript plot_merge_bam_stats.R {wildcards.species}.bam.table.tsv {wildcards.species}.depths.boxplot.pdf {wildcards.species}.flagstats.boxplot.pdf 
		cp {wildcards.species}.depths.boxplot.pdf {output.plot_depth}
		cp {wildcards.species}.flagstats.boxplot.pdf {output.plot_flagstats}
		"""
