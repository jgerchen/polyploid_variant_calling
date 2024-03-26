import gzip
import subprocess
import argparse
import re
parser = argparse.ArgumentParser()

parser.add_argument('--depth', nargs='+')
parser.add_argument('--clean_depth', nargs='+')
parser.add_argument('--flagstats', nargs='+')
parser.add_argument('-o', '--output_file')
args=parser.parse_args()
chr_file_list=args.depth
clean_chr_file_list=args.clean_depth
flagstat_list=args.flagstats
output_file=args.output_file

sample_stats_dict={}

for sample_flagstat in flagstat_list:
	sample_name=sample_flagstat.split("/")[-1].split("_")[-1].replace(".flagstats.tsv", "")
	s_dict={}
	with open(sample_flagstat) as s_flagstat_input:
		for line in s_flagstat_input:
			s_fs_value=line.strip().split("\t")[0]
			s_fs_cat=line.strip().split("\t")[2]
			if "total" in s_fs_cat:
				s_dict.update({"totalreads":s_fs_value})
			elif s_fs_cat=="mapped":
				s_dict.update({"mapped":s_fs_value})
			elif s_fs_cat=="mapped %":
				s_dict.update({"perc_mapped":s_fs_value})
			elif s_fs_cat=="duplicates":
				s_dict.update({"duplicates":s_fs_value})
			elif s_fs_cat=="paired in sequencing":
				s_dict.update({"paired":s_fs_value})
			elif s_fs_cat=="properly paired":
				s_dict.update({"properly_paired":s_fs_value})
			elif s_fs_cat=="properly paired %":
				s_dict.update({"perc_properly_paired":s_fs_value})
			elif s_fs_cat=="singletons":
				s_dict.update({"singletons":s_fs_value})
			elif s_fs_cat=="singletons %":
				s_dict.update({"perc_singletons":s_fs_value})
			elif s_fs_cat=="with mate mapped to a different chr":
				s_dict.update({"mate_diff_chr":s_fs_value})
			elif s_fs_cat=="with mate mapped to a different chr (mapQ>=5)":
				s_dict.update({"mate_diff_chr_mq5":s_fs_value})
		s_dict.update({"perc_duplicates":str(round((int(s_dict["duplicates"])/int(s_dict["mapped"]))*100, 2))+"%"})

		sample_re=re.compile(".*"+sample_name)
		sample_chr_matches=list(filter(sample_re.match, chr_file_list))
		assert len(sample_chr_matches)==1, "Error: there should be exactly one raw depth file matching for sample %s, but it's %s" % (sample_name, len(sample_chr_matches))
		with gzip.open(sample_chr_matches[0], "rt") as depth_raw_file:
			for dr_line in depth_raw_file:
				if dr_line[0:2]=="##":
					dr_dict={dr_i.split()[0]:dr_i.split()[1] for dr_i in dr_line.strip().split("\t")}
					s_dict.update({"raw_mean_depth":dr_dict["MeanDepth:"]})
					s_dict.update({"raw_perc_cov":dr_dict["Coverage(%):"]})
		sample_chr_clean_matches=list(filter(sample_re.match, clean_chr_file_list))
		assert len(sample_chr_clean_matches)==1, "Error: there should be exactly one clean depth file matching for sample %s, but it's %s" % (sample_name, len(sample_chr_clean_matches))
		with gzip.open(sample_chr_clean_matches[0], "rt") as depth_clean_file:
			for dc_line in depth_clean_file:
				if dc_line[0:2]=="##":
					dc_dict={dc_i.split()[0]:dc_i.split()[1] for dc_i in dc_line.strip().split("\t")}
					s_dict.update({"clean_mean_depth":dc_dict["MeanDepth:"]})
					s_dict.update({"clean_perc_cov":dc_dict["Coverage(%):"]})
	sample_stats_dict.update({sample_name:s_dict})
with open(output_file,"w") as output_table:
	output_table.write("sample\traw_mean_depth\traw_perc_cov\tclean_mean_depth\tclean_perc_cov\ttotal\tmapped\tperc_mapped\tduplicates\tperc_duplicates\tpaired\tproperly_paired\tperc_properly_paired\tsingletons\tperc_singletons\tmate_diff_chr\tmate_diff_chr_mq5\n")
	for out_sample in list(sample_stats_dict.keys()):
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (out_sample, sample_stats_dict[out_sample]["raw_mean_depth"], sample_stats_dict[out_sample]["raw_perc_cov"],sample_stats_dict[out_sample]["clean_mean_depth"], sample_stats_dict[out_sample]["clean_perc_cov"], sample_stats_dict[out_sample]["totalreads"], sample_stats_dict[out_sample]["mapped"], sample_stats_dict[out_sample]["perc_mapped"], sample_stats_dict[out_sample]["duplicates"], sample_stats_dict[out_sample]["perc_duplicates"], sample_stats_dict[out_sample]["paired"], sample_stats_dict[out_sample]["properly_paired"], sample_stats_dict[out_sample]["perc_properly_paired"],sample_stats_dict[out_sample]["singletons"],sample_stats_dict[out_sample]["perc_singletons"],sample_stats_dict[out_sample]["mate_diff_chr"],sample_stats_dict[out_sample]["mate_diff_chr_mq5"] ))
#rplot = subprocess.call("Rscript scripts/plot_merge_bam_stats.R %s %s %s" % (output.stats_table, output.plot_depth, output.plot_flagstats), shell=True)
