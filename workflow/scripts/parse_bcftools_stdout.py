import argparse
import sys
import numpy
import matplotlib.pyplot as plt

parser=argparse.ArgumentParser()
parser.add_argument('--n_sites')
parser.add_argument('--output')
parser.add_argument('--histogram_bins')
parser.add_argument('--invariants', action='store_true')
parser.add_argument('--inv_INFO_fields', nargs='*')
parser.add_argument('--inv_GT_fields', nargs='*')
parser.add_argument('--biallelic', action='store_true')
parser.add_argument('--bi_INFO_fields', nargs='*')
parser.add_argument('--bi_GT_fields', nargs='*')
parser.add_argument('--multiallelic', action='store_true')
parser.add_argument('--ma_INFO_fields', nargs='*')
parser.add_argument('--ma_GT_fields', nargs='*')
args=parser.parse_args()

n_sites=int(args.n_sites)
hist_bins=int(args.histogram_bins)

output=args.output
check_invariants=args.invariants
if args.inv_INFO_fields and len(args.inv_INFO_fields)>0:
	get_info_stats_inv=True
	inv_info_stats=args.inv_INFO_fields
else:
	get_info_stats_inv=False
if args.inv_GT_fields and len(args.inv_GT_fields)>0:
	get_gt_stats_inv=True
	inv_gt_stats=args.inv_GT_fields
else:
	get_gt_stats_inv=False

check_biallelic=args.biallelic
if args.bi_INFO_fields and len(args.bi_INFO_fields)>0:
	get_info_stats_bi=True
	bi_info_stats=args.bi_INFO_fields
else:
	get_info_stats_bi=False
if args.bi_GT_fields and len(args.bi_GT_fields)>0:
	get_gt_stats_bi=True
	bi_gt_stats=args.bi_GT_fields
else:
	get_gt_stats_bi=False

check_multiallelic=args.multiallelic
if args.bi_INFO_fields and len(args.bi_INFO_fields)>0:
	get_info_stats_ma=True
	ma_info_stats=args.ma_INFO_fields
else:
	get_info_stats_ma=False
if args.ma_GT_fields and len(args.ma_GT_fields)>0:
	get_gt_stats_ma=True
	ma_gt_stats=args.bi_GT_fields
else:
	get_gt_stats_ma=False

def get_qual(s_qual, qual_count_dict, qual_list):
	if s_qual=="Infinity" or s_qual=="inf":
		if "Infinity" not in qual_count_dict:
			qual_count_dict.update({"Infinity":1})
		else:
			qual_count_dict["Infinity"]+=1
	else:
		try:
			qual_fl=float(s_qual)
			if numpy.isinf(qual_fl):
				print("%s equals to inf" % s_qual)
			qual_list[qual_count_dict["float"]]=qual_fl
			qual_count_dict["float"]+=1
		except:			
			if s_qual not in qual_count_dict:
				qual_count_dict.update({s_qual:1})
			else:
				qual_count_dict[s_qual]+=1

def get_info_stats(s_info, info_dict, info_counts):
	info_cats=s_info.split(";")
	for info_cat in info_cats:
		info_cat_subs=info_cat.split("=")
		if info_cat_subs[0] in info_dict:
			info_dict[info_cat_subs[0]][info_counts[info_cat_subs[0]]]=info_cat_subs[1]
			info_counts[info_cat_subs[0]]+=1

def get_gt_stats(s_gt_info, s_gts, gt_gen_dict, gt_dict, gt_counts):
	gt_cats=s_gt_info.split(":")
	for gt_i in range(len(s_gts)):
		gt_subs=s_gts[gt_i].split(":")
		gt_gt=gt_subs[0]
		#counts genotypes
		if gt_gt not in gt_gen_dict:
			gt_gen_dict.update({gt_gt:[0]*len(s_gts)})
		gt_gen_dict[gt_gt][gt_i]+=1
		#should collect genotype stats
		#{STAT1:[[sample1,sample1,sample1],[sample2,sample2,sample2],[sample3,sample3,sample3]], STAT2:[[sample1,sample1,sample1],[sample2,sample2,sample2],[sample3,sample3,sample3]]  }

		for gt_cat_i in range(len(gt_cats)):
			if gt_cats[gt_cat_i] in gt_dict:
				try:
					gt_dict[gt_cats[gt_cat_i]][gt_i][gt_counts[gt_cats[gt_cat_i]][gt_i]]=gt_subs[gt_cat_i]
					#{STAT: ok                 sample   site                           }	
					#gt_counts={STAT1:[sample1, sample2, sample3], STAT2:[sample1, sample2, sample3] }

					gt_counts[gt_cats[gt_cat_i]][gt_i]+=1
				except ValueError:
						
					pass
if check_invariants==True:
	invariant_counter=0
	invariant_long_ref_dict={}
	invariant_filter_dict={}
	invariant_qual_list=numpy.empty(n_sites, dtype=float)
	invariant_qual_count_dict={"float":0}

	if get_info_stats_inv==True:
		inv_info_dict={inv_info_stat.split("_")[0]:numpy.empty(n_sites, dtype=inv_info_stat.split("_")[1]) for inv_info_stat in inv_info_stats}
		inv_info_counts={inv_info_count.split("_")[0]:0 for inv_info_count in inv_info_stats}
	
#biallelic
if check_biallelic==True:
	biallelic_counter=0
	biallelic_long_ref_dict={}
	biallelic_long_alt_dict={}
	biallelic_filter_dict={}
	biallelic_qual_list=numpy.empty(n_sites, dtype=float)
	biallelic_qual_count_dict={"float":0}
	if get_info_stats_bi==True:
		bi_info_dict={bi_info_stat.split("_")[0]:numpy.empty(n_sites, dtype=bi_info_stat.split("_")[1]) for bi_info_stat in bi_info_stats}
		bi_info_counts={bi_info_count.split("_")[0]:0 for bi_info_count in bi_info_stats}
	
#multiallellic
if check_multiallelic==True:
	multiallelic_counter=0
	multiallelic_long_ref_dict={}
	multiallelic_long_alt_dict={}
	multiallelic_filter_dict={}
	multiallelic_n_alt_alleles={2:0, 3:0, 4:0,5:0,6:0,7:0,8:0,9:0,10:0}
	multiallelic_qual_list=numpy.empty(n_sites, dtype=float)
	multiallelic_qual_count_dict={"float":0}
	if get_info_stats_ma==True:
		ma_info_dict={ma_info_stat.split("_")[0]:numpy.empty(n_sites, dtype=ma_info_stat.split("_")[1]) for ma_info_stat in ma_info_stats}
		ma_info_counts={ma_info_count.split("_")[0]:0 for ma_info_count in ma_info_stats}
	

line = sys.stdin.readline()
line_counter=0
while line:
	if line[0]=="#":
		if line[0:6]=="#CHROM":
			sample_list=line.strip().split("\t")[9:]
			n_samples=len(sample_list)		
			if get_gt_stats_inv==True:
				inv_gt_gen_dict={}
				inv_gt_dict={inv_gt_stat.split("_")[0]:[numpy.empty(n_sites, dtype=inv_gt_stat.split("_")[1]) for sample in sample_list] for inv_gt_stat in inv_gt_stats}
				inv_gt_counts={inv_gt_count.split("_")[0]:[0 for u in sample_list] for inv_gt_count in inv_gt_stats}
			if get_gt_stats_bi==True:
				bi_gt_gen_dict={}
				bi_gt_dict={bi_gt_stat.split("_")[0]:[numpy.empty(n_sites, dtype=bi_gt_stat.split("_")[1]) for sample in sample_list] for bi_gt_stat in bi_gt_stats}
				bi_gt_counts={bi_gt_count.split("_")[0]:[0 for u in sample_list] for bi_gt_count in bi_gt_stats}
			if get_gt_stats_ma==True:
				ma_gt_gen_dict={}
				ma_gt_dict={ma_gt_stat.split("_")[0]:[numpy.empty(n_sites, dtype=ma_gt_stat.split("_")[1]) for sample in sample_list] for ma_gt_stat in ma_gt_stats}
				ma_gt_counts={ma_gt_count.split("_")[0]:[0 for u in sample_list] for ma_gt_count in ma_gt_stats}

	else:
		line_cats=line.strip().split("\t")
		assert(len(line_cats)>8),"Variant line sould be at least 8 categories, but line %s has only %s" % (line, len(line_cats))
		line_counter+=1
		ref=line_cats[3]
		alt=line_cats[4]
		qual=line_cats[5]
		filter_cat=line_cats[6]
		info=line_cats[7]
		gt_info=line_cats[8]
		gts=line_cats[9:]
		#invariant stuff here
		if check_invariants==True:
			if alt==".":
				#invariant stuff here
				invariant_counter+=1
				l_ref=len(ref)
				if l_ref not in invariant_long_ref_dict:
					invariant_long_ref_dict.update({l_ref:1})
				else:
					invariant_long_ref_dict[l_ref]+=1
				if filter_cat not in invariant_filter_dict:
					invariant_filter_dict.update({filter_cat:1})
				else:
					invariant_filter_dict[filter_cat]+=1
				get_qual(qual, invariant_qual_count_dict, invariant_qual_list)
				if get_info_stats_inv==True:
					get_info_stats(info, inv_info_dict, inv_info_counts)	
				if get_gt_stats_inv==True:
					get_gt_stats(gt_info, gts, inv_gt_gen_dict, inv_gt_dict, inv_gt_counts)
#		#multivariant stuff here
		if check_multiallelic==True:
			if  alt!="." and len(alt.split(","))>1:
				multiallelic_counter+=1
				n_alt_alleles=len(alt.split(","))
				if len(ref) not in multiallelic_long_ref_dict:
					multiallelic_long_ref_dict.update({len(ref):1})
				else:
					multiallelic_long_ref_dict[len(ref)]+=1

				for alt_allele in alt.split(","):
					if len(alt_allele) not in multiallelic_long_alt_dict:
						multiallelic_long_alt_dict.update({len(alt_allele):1})
					else:
						multiallelic_long_alt_dict[len(alt_allele)]+=1
				if n_alt_alleles<=10:
					multiallelic_n_alt_alleles[n_alt_alleles]+=1
				else:
					multiallelic_n_alt_alleles[10]+=1
				if filter_cat not in multiallelic_filter_dict:
					multiallelic_filter_dict.update({filter_cat:1})
				else:
					multiallelic_filter_dict[filter_cat]+=1
				get_qual(qual, multiallelic_qual_count_dict, multiallelic_qual_list)
				if get_info_stats_ma==True:
					get_info_stats(info, ma_info_dict, ma_info_counts)	
				if get_gt_stats_ma==True:
					get_gt_stats(gt_info, gts, ma_gt_gen_dict, ma_gt_dict, ma_gt_counts)
		if check_biallelic==True:
			if alt!="." and len(alt.split(","))==1:
				biallelic_counter+=1
				if len(ref) not in biallelic_long_ref_dict:
					biallelic_long_ref_dict.update({len(ref):1})
				else:
					biallelic_long_ref_dict[len(ref)]+=1
				if len(alt) not in biallelic_long_alt_dict:
					biallelic_long_alt_dict.update({len(alt):1})
				else:
					biallelic_long_alt_dict[len(alt)]+=1

				if filter_cat not in biallelic_filter_dict:
					biallelic_filter_dict.update({filter_cat:1})
				else:
					biallelic_filter_dict[filter_cat]+=1

				get_qual(qual, biallelic_qual_count_dict, biallelic_qual_list)
				if get_info_stats_bi==True:
					get_info_stats(info, bi_info_dict, bi_info_counts)	
				if get_gt_stats_bi==True:
					get_gt_stats(gt_info, gts, bi_gt_gen_dict, bi_gt_dict, bi_gt_counts)
	
	line = sys.stdin.readline()

#Output functions

def write_histogram(hist_array, var_type, statistic_category, statistic_type, sample, dict_output_table):
	hist_mean=numpy.mean(hist_array)
	output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(var_type, statistic_category,statistic_type+"_mean", "NA", "NA", sample,hist_mean))
	hist_median=numpy.median(hist_array)
	output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(var_type, statistic_category,statistic_type+"_median", "NA", "NA", sample,hist_median))
	hist_std=numpy.std(hist_array)
	output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(var_type, statistic_category,statistic_type+"_std", "NA", "NA", sample,hist_std))

	hist=numpy.histogram(hist_array, hist_bins)
	for hist_bin in range(len(hist[0])):
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(var_type, statistic_category, statistic_type+"_histogram", hist[1][hist_bin], hist[1][hist_bin+1],sample, hist[0][hist_bin]))

def plot_hist(plot_array, plot_title ,plot_file, number_bins):
	fig,ax=plt.subplots()
	plt.hist(plot_array, number_bins, log=True)
	ax.set_title(plot_title)
	fig.tight_layout()
	plt.savefig(plot_file)

def plot_gt_violin(violin_arrays, violin_counts, violin_title, violin_samples, violin_file):
	plot_arrays=numpy.array([violin_arrays[vc][0:violin_counts[vc]] for vc in range(len(violin_counts))], dtype=object, copy=False)
	fig,ax=plt.subplots()
	ax.violinplot(plot_arrays, showmeans=True, showmedians=True, showextrema=True, widths=0.7)
	ax.set_title(violin_title)
	ax.set_xticks(numpy.arange(1, len(violin_samples) + 1), labels=violin_samples)
	ax.set_xlim(0.25, len(violin_samples) + 0.75)
	ax.set_xlabel('Sample')
	plt.savefig(violin_file)

def plot_gt_types(gt_type_gen_dict, gt_type_samples, gt_type_title, gt_type_file):
	fig, ax = plt.subplots(layout="constrained")
	x = numpy.arange(len(gt_type_samples))
	width = 1/(len(gt_type_gen_dict)+1)  # the width of the bars
	multiplier = 0
	for gt_type in sorted(list(gt_type_gen_dict.keys()), key=len):
		offset=width*multiplier
		rects =ax.bar(x + offset, gt_type_gen_dict[gt_type], width, label=gt_type, log=True)
		ax.bar_label(rects, labels=[gt_type for u in gt_type_samples], padding=3, fontsize=6, rotation='vertical')
		multiplier += 1
	ax.set_ylabel('GT count')
	ax.set_title(gt_type_title)
	ax.set_xticks(x + (width*(len(gt_type_gen_dict)/2)), gt_type_samples)
	plt.gcf().set_size_inches(len(gt_type_samples)*3, 5)
	plt.savefig(gt_type_file,dpi=200)

def plot_general_data(general_title, general_sub_titles, general_out_file, *general_input_dicts):
	inp_dicts=[inp_dict for inp_dict in general_input_dicts]
	fig,ax=plt.subplots(len(inp_dicts), layout="constrained")
	fig.suptitle(general_title)
	for inp_dict_i in range(len(inp_dicts)):
		general_plot_cats=sorted(inp_dicts[inp_dict_i])
		general_plot_values=[inp_dicts[inp_dict_i][gp_cat] for gp_cat in general_plot_cats]
		ax[inp_dict_i].bar(general_plot_cats, general_plot_values, log=True)
		ax[inp_dict_i].set_title(general_sub_titles[inp_dict_i])
	plt.gcf().set_size_inches(8,len(inp_dicts)*2.5)
	plt.savefig(general_out_file)

def plot_info_stats(info_title, info_dict, info_counts, info_out_file):
	n_info_stats=len(info_counts)
	info_stats=[i for i in info_counts]
	if n_info_stats>1:
		fig,axs=plt.subplots(n_info_stats, layout="constrained")
		fig.suptitle(info_title)
		for info_stat_i in range(n_info_stats):	
			axs[info_stat_i].hist(info_dict[info_stats[info_stat_i]][0:info_counts[info_stats[info_stat_i]]], bins=hist_bins, log=True)
			axs[info_stat_i].set_title(info_stats[info_stat_i]+": "+str(info_counts[info_stats[info_stat_i]])+" observations" )
		plt.gcf().set_size_inches(8,n_info_stats*2.5)
	else:
		fig,ax=plt.subplots()
		ax.hist(info_dict[info_stats[0]][0:info_counts[info_stats[0]]], bins=hist_bins, log=True)
		fig.suptitle(info_title)
		ax.set_title(info_stats[0]+": "+str(info_counts[info_stats[0]])+" observations" )
	plt.savefig(info_out_file)

def write_output_info(var_type, info_dict, info_counts ,info_out_file):
	for info_stat in info_counts:
		info_out_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(var_type, "INFO",info_stat+"_count", "NA", "NA","all" ,info_counts[info_stat]))
		write_histogram(info_dict[info_stat][0:info_counts[info_stat]],var_type, "INFO", info_stat, "all", info_out_file)

def write_dict(w_dict, var_type, statistic_category, statistic_type, sample, dict_output_table):
	for w_dict_entry in sorted(w_dict):	
		dict_output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(var_type, statistic_category, statistic_type, w_dict_entry, "NA", sample ,w_dict[w_dict_entry]))

def write_output_gt(var_type, gt_dict, gt_gen_dict, gt_counts, gt_samples, gt_out_file):
	for gt_sample_i in range(len(gt_samples)):
		for gt_type in gt_gen_dict:
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(var_type, "GT","GT_count", gt_type, "NA", gt_samples[gt_sample_i], gt_gen_dict[gt_type][gt_sample_i]))
	for gt_stat in gt_counts:
		for gt_sample_i in range(len(gt_counts[gt_stat])):
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(var_type, "GT",gt_stat+"_count", "NA", "NA",sample_list[gt_sample_i],gt_counts[gt_stat][gt_sample_i]))
			write_histogram(gt_dict[gt_stat][gt_sample_i][0:gt_counts[gt_stat][gt_sample_i]], var_type ,"GT", gt_stat, gt_samples[gt_sample_i], gt_out_file)



#write output file
with open(output+"_table.tsv","w") as output_table:
	output_table.write("var_type\tstatistic_category\tstatistic_type\tvar_or_bin_start\tbin_end\tsample\tstatistic\n")
	#biallelic sites
	if check_biallelic==True:
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "general", "site_count" ,"NA", "NA","all",biallelic_counter))
		write_dict(biallelic_long_ref_dict, "biallelic", "general", "ref_length_count", "all", output_table)
		write_dict(biallelic_long_alt_dict, "biallelic", "general", "alt_length_count", "all", output_table)
		write_dict(biallelic_filter_dict, "biallelic", "general", "filter_category_count", "all", output_table)
		write_dict(biallelic_qual_count_dict, "biallelic", "general", "qual_category_count", "all", output_table)
		#plot general data
		plot_general_data("Biallelic: %s sites" % biallelic_counter, ["Ref length", "Alt length", "filter categories", "QUAL categories"], "biallelic_general.pdf", biallelic_long_ref_dict, biallelic_long_alt_dict, biallelic_filter_dict, biallelic_qual_count_dict)
		#plot qual histogram
		write_histogram(biallelic_qual_list[0:biallelic_qual_count_dict["float"]], "biallelic", "general", "qual_float_histogram", "all", output_table)
		plot_hist(biallelic_qual_list[0:biallelic_qual_count_dict["float"]], "QUAL" ,output+"_QUAL_biallelic.pdf", hist_bins)

	#info stats
		if get_info_stats_bi==True:
			write_output_info("biallelic", bi_info_dict, bi_info_counts, output_table)
			plot_info_stats("biallelic", bi_info_dict, bi_info_counts, output+"_INFO_biallelic.pdf")
		#GT stats
		if get_gt_stats_bi==True:
			write_output_gt("biallelic", bi_gt_dict, bi_gt_gen_dict, bi_gt_counts, sample_list, output_table)
			plot_gt_types(bi_gt_gen_dict, sample_list, "Genotype counts", output+"_GT_counts_biallelic.pdf")
			for bi_gt_stat in bi_gt_counts:
				plot_gt_violin(bi_gt_dict[bi_gt_stat], bi_gt_counts[bi_gt_stat], bi_gt_stat, sample_list, output+"_"+bi_gt_stat+"_comp_biallelic.pdf")

	if check_invariants==True:
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("invariant", "general", "site_count" ,"NA", "NA","all",invariant_counter))
		write_dict(invariant_long_ref_dict, "invariant", "general", "ref_length_count", "all", output_table)
		write_dict(invariant_filter_dict, "invariant", "general", "filter_category_count", "all", output_table)
		write_dict(invariant_qual_count_dict, "invariant", "general", "qual_category_count", "all", output_table)
		plot_general_data("Invariant: %s sites" % invariant_counter, ["Ref length", "filter categories", "QUAL categories"], "invariant_general.pdf", invariant_long_ref_dict,invariant_filter_dict, invariant_qual_count_dict)
		#plot qual histogram
		write_histogram(invariant_qual_list[0:invariant_qual_count_dict["float"]], "invariant", "general", "qual_float_histogram", "all", output_table)
		plot_hist(invariant_qual_list[0:invariant_qual_count_dict["float"]], "QUAL" ,output+"_QUAL_invariant.pdf", hist_bins)
		if get_info_stats_inv==True:
			write_output_info("invariant", inv_info_dict, inv_info_counts, output_table)
			plot_info_stats("invariant", inv_info_dict, inv_info_counts, output+"_INFO_invariant.pdf")
		if get_gt_stats_inv==True:
			write_output_gt("invariant", inv_gt_dict, inv_gt_gen_dict, inv_gt_counts, sample_list, output_table)
			#plot_gt_types(inv_gt_gen_dict, sample_list, "Genotype counts", output+"_GT_counts_invariant.pdf")
			for inv_gt_stat in inv_gt_counts:
				plot_gt_violin(inv_gt_dict[inv_gt_stat], inv_gt_counts[inv_gt_stat], inv_gt_stat, sample_list, output+"_"+inv_gt_stat+"_comp_invariant.pdf")

	if check_multiallelic==True:
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "general", "site_count" ,"NA", "NA","all",multiallelic_counter))
		write_dict(multiallelic_long_ref_dict, "multiallelic", "general", "ref_length_count", "all", output_table)
		write_dict(multiallelic_long_alt_dict, "multiallelic", "general", "alt_length_count", "all", output_table)
		write_dict(multiallelic_n_alt_alleles, "multiallelic", "general", "number_alt_alleles", "all", output_table)
		write_dict(multiallelic_filter_dict, "multiallelic", "general", "filter_category_count", "all", output_table)
		write_dict(multiallelic_qual_count_dict, "multiallelic", "general", "qual_category_count", "all", output_table)
		plot_general_data("Multiallelic: %s sites" % multiallelic_counter, ["Ref length", "Alt length","Number of alt alleles", "filter categories", "QUAL categories"], "multiallelic_general.pdf", multiallelic_long_ref_dict, multiallelic_long_alt_dict, multiallelic_n_alt_alleles, multiallelic_filter_dict, multiallelic_qual_count_dict)
		write_histogram(multiallelic_qual_list[0:multiallelic_qual_count_dict["float"]], "multiallelic", "general", "qual_float_histogram", "all", output_table)
		plot_hist(multiallelic_qual_list[0:multiallelic_qual_count_dict["float"]], "QUAL" ,output+"_QUAL_multiallelic.pdf", hist_bins)
		if get_info_stats_ma==True:
			write_output_info("multiallelic", ma_info_dict, ma_info_counts, output_table)
			plot_info_stats("multiallelic", ma_info_dict, ma_info_counts, output+"_INFO_multiallelic.pdf")
		if get_gt_stats_ma==True:
			write_output_gt("multiallelic", ma_gt_dict, ma_gt_gen_dict, ma_gt_counts, sample_list, output_table)
			for ma_gt_stat in ma_gt_counts:
				plot_gt_violin(ma_gt_dict[ma_gt_stat], ma_gt_counts[ma_gt_stat], ma_gt_stat, sample_list, output+"_"+ma_gt_stat+"_comp_multiallelic.pdf")


