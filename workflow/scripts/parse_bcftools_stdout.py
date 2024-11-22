import argparse
import sys
import numpy
import numpy.ma as ma
import matplotlib.pyplot as plt
import time
import itertools

parser=argparse.ArgumentParser()
parser.add_argument('--n_sites')
parser.add_argument('--output')
parser.add_argument('--histogram_bins')
parser.add_argument('--invariants', action='store_true')
parser.add_argument('--biallelic', action='store_true')
parser.add_argument('--multiallelic', action='store_true')
parser.add_argument('--plot_filter', action='store_true')
args=parser.parse_args()

n_sites=int(args.n_sites)
hist_bins=int(args.histogram_bins)
output=args.output
check_invariants=args.invariants
check_biallelic=args.biallelic
check_multiallelic=args.multiallelic
plot_filter=args.plot_filter
line_counter=0

#print("Allocating main array...")
main_array_time=time.time()
#create arrays
#number of alt alleles
n_alt_array=numpy.zeros(n_sites, dtype=numpy.uint16)
#QUAL
QUAL_array=numpy.empty(n_sites)
QUAL_array.fill(numpy.nan)
if check_biallelic==True:
	QUAL_count_dict_bi={}
if check_invariants==True:
	QUAL_count_dict_inv={}
if check_multiallelic==True:
	QUAL_count_dict_ma={}
if plot_filter==True:
	filter_dict={}
#Not included: FIlters -> we typically only use hard filtering


#INFO do we need more stats???
#Invariants only DP, uint32 with 0 for missing?
info_DP_array=numpy.zeros(n_sites, dtype=numpy.uint32)

#others DP+GATK best practices: DP, FS, QD, MQ,
if check_biallelic==True or check_multiallelic==True:
	info_array_GATK_bp=numpy.empty((n_sites, 6))
	info_array_GATK_bp.fill(numpy.nan)
	GATK_index_array={"FS":0, "QD":1, "MQ":2, "MQRankSum":3, "ReadPosRankSum":4, "SOR":5}

#print("Done allocating main arrays in %s s" % str(time.time()-main_array_time))

#print("Reading vcf file...")
vcf_read_time=time.time()

line = sys.stdin.readline()
while line:
	if line[0]=="#":
		if line[0:6]=="#CHROM":
			sample_list=line.strip().split("\t")[9:]
			n_samples=len(sample_list)		
			#DEPTH array					
			#GT_depth_array=numpy.zeros((n_sites, n_samples), dtype=numpy.uint32)
			#DEPTH count dicts instead! But for each type of variant
			#Try to put GQ/RGQ in arrays! No! use count dicts!




				#GQ_array=numpy.empty((n_sites, n_samples))
				#GQ_array.fill(numpy.nan)
				#RGQ_array=numpy.empty((n_sites, n_samples))
				#RGQ_array.fill(numpy.nan)

			if check_biallelic==True:
			#	gt_bi_GQ_dict=[{i_gq:0 for i_gq in range(100)} for i_s in sample_list]
				gt_bi_GT_dict=[{} for i_s in sample_list]
				gt_bi_depth_count_dict=[{d:0 for d in range(100)} for u in range(n_samples)]
				gt_bi_max_depth=0
				gt_bi_GQ_count_dict=[{d:0 for d in range(100)} for u in range(n_samples)]
			if check_multiallelic==True:
			#	gt_ma_GQ_dict=[{i_gq:0 for i_gq in range(100)} for i_s in sample_list]
				gt_ma_GT_dict=[{} for i_s in sample_list]
				gt_ma_depth_count_dict=[{d:0 for d in range(100)} for u in range(n_samples)]
				gt_ma_max_depth=0
				gt_ma_GQ_count_dict=[{d:0 for d in range(100)} for u in range(n_samples)]
			if check_invariants==True:
			#	gt_inv_RGQ_dict=gt_bi_GQ_dict=[{i_gq:0 for i_gq in range(100)} for i_s in sample_list]
				gt_inv_GT_dict=[{} for i_s in sample_list]
				gt_inv_depth_count_dict=[{d:0 for d in range(100)} for u in range(n_samples)]
				gt_inv_max_depth=0
				gt_inv_RGQ_count_dict=[{d:0 for d in range(100)} for u in range(n_samples)]

	else:
		line_cats=line.strip().split("\t")
		assert(len(line_cats)>8),"Variant line sould be at least 8 categories, but line %s has only %s" % (line, len(line_cats))
		ref=line_cats[3]
		alt=line_cats[4]
		if alt==".":
			n_alt_array[line_counter]=0
			n_alt_alleles=0
		else:
			n_alt_alleles=len(alt.split(","))
			n_alt_array[line_counter]=n_alt_alleles			
		#If no alt alleles site must be invariant

		qual=line_cats[5]
		#parse QUAL...leave inf QUALS in main array
		try:
			QUAL_array[line_counter]=float(qual)
		except:
			if n_alt_alleles==0:
				if check_invariants==True:
					if qual in QUAL_count_dict_inv:
						QUAL_count_dict_inv[qual]+=1
					else:	
						QUAL_count_dict_inv.update({qual:1})
			elif n_alt_alleles==1:
				if check_biallelic==True:
					if qual in QUAL_count_dict_bi:
						QUAL_count_dict_bi[qual]+=1
					else:	
						QUAL_count_dict_bi.update({qual:1})
			else:
				if check_multiallelic==True:
					if qual in QUAL_count_dict_ma:
						QUAL_count_dict_ma[qual]+=1
					else:	
						QUAL_count_dict_ma.update({qual:1})

		filter_cat=line_cats[6]
		if plot_filter:
			if filter_cat not in filter_dict:
				filter_dict.update({filter_cat:1})
			else:
				filter_dict[filter_cat]+=1
		#parse INFO
		info=line_cats[7]
		info_cats=info.split(";")
		for info_cat in info_cats:
			info_cat_subs=info_cat.split("=")
			if info_cat_subs[0]=="DP":
				info_DP_array[line_counter]=int(info_cat_subs[1])
			elif (check_biallelic==True or check_multiallelic==True) and info_cat_subs[0] in GATK_index_array:
				info_array_GATK_bp[line_counter][GATK_index_array[info_cat_subs[0]]]=float(info_cat_subs[1])
		gt_info=line_cats[8].split(":")
		try:
			gt_DP_i=gt_info.index("DP")
		except:
			gt_DP_i=None
		try:
			gt_GQ_i=gt_info.index("GQ")
		except:
			gt_GQ_i=None
		try:
			gt_RGQ_i=gt_info.index("RGQ")
		except:
			gt_RGQ_i=None

		gts=line_cats[9:]
		#Iterate over GTs once
		for gt_i in range(len(gts)):
			gt_all=gts[gt_i].split(":")
			#get DP
			try:
				gt_depth=int(gt_all[gt_DP_i])
			except:
				gt_depth=0	
			#GT_depth_array[line_counter][gt_i]=gt_depth
			#if gt_depth in GT_depth_count_dicts[gt_i]:
			#	GT_depth_count_dicts[gt_i][gt_depth]+=1
			#else:	
			#	GT_depth_count_dicts[gt_i].update({gt_depth:1})
			gt_GT=gt_all[0]
			#test if invariant, bi or ma
			if n_alt_alleles==0:
				if check_invariants==True:
					if gt_GT not in gt_inv_GT_dict[gt_i]:
						gt_inv_GT_dict[gt_i].update({gt_GT:1})
					else:
						gt_inv_GT_dict[gt_i][gt_GT]+=1
					if gt_depth in gt_inv_depth_count_dict[gt_i]:
						gt_inv_depth_count_dict[gt_i][gt_depth]+=1
					else:	
						gt_inv_depth_count_dict[gt_i].update({gt_depth:1})
					#if gt_depth>gt_inv_max_depth:
					#	gt_inv_max_depth=gt_depth
					#	try:
					#		gt_inv_RGQ_dict[gt_i]=int(gt_all[gt_RGQ_i])
							
					#	except:
					#		pass
					if gt_RGQ_i:
						try:
							#RGQ_array[line_counter][gt_i]=float(gt_all[gt_RGQ_i])
							gt_inv_RGQ_count_dict[gt_i][int(gt_all[gt_RGQ_i])]+=1
						except:
							pass
			elif n_alt_alleles==1:
				if check_biallelic==True:
					if gt_GT not in gt_bi_GT_dict[gt_i]:
						gt_bi_GT_dict[gt_i].update({gt_GT:1})
					else:
						gt_bi_GT_dict[gt_i][gt_GT]+=1
					if gt_depth in gt_bi_depth_count_dict[gt_i]:
						gt_bi_depth_count_dict[gt_i][gt_depth]+=1
					else:	
						gt_bi_depth_count_dict[gt_i].update({gt_depth:1})
					#if gt_depth>gt_bi_max_depth:
					#	gt_bi_max_depth=gt_depth
					#if gt_GQ_i:
						#try:
						#	gt_bi_GQ_dict[gt_i]=int(gt_all[gt_GQ_i])
						#except:
						#	pass
					if gt_GQ_i:
						try:
							#GQ_array[line_counter][gt_i]=float(gt_all[gt_GQ_i])
							gt_bi_GQ_count_dict[gt_i][int(gt_all[gt_GQ_i])]+=1
						except:
							pass
			else:
				if check_multiallelic==True:
					if gt_GT not in gt_ma_GT_dict[gt_i]:
						gt_ma_GT_dict[gt_i].update({gt_GT:1})
					else:
						gt_ma_GT_dict[gt_i][gt_GT]+=1
					if gt_depth in gt_ma_depth_count_dict[gt_i]:
						gt_ma_depth_count_dict[gt_i][gt_depth]+=1
					else:	
						gt_ma_depth_count_dict[gt_i].update({gt_depth:1})
					#if gt_depth>gt_ma_max_depth:
					#	gt_ma_max_depth=gt_depth
					if gt_GQ_i:
						try:
							#GQ_array[line_counter][gt_i]=float(gt_all[gt_GQ_i])
							gt_ma_GQ_count_dict[gt_i][int(gt_all[gt_GQ_i])]+=1
						except:
							pass
					#if gt_GQ_i:
					#	try:
					#		gt_ma_GQ_dict[gt_i]=int(gt_all[gt_GQ_i])
					#	except:
					#		pass
				
		line_counter+=1
	line = sys.stdin.readline()

#print("Done reading VCF in %s s" % str(time.time()-vcf_read_time))




##########function to estimate mean, SD, median and histogram from count dictionaries

def analyze_count_dict(input_dict, n_bins, max_bin):
	i_dict_n_items=0
	i_dict_total=0
	#mean
	for i_dict_i in input_dict:
		i_dict_n_items+=input_dict[i_dict_i]
		i_dict_total+=i_dict_i*input_dict[i_dict_i]
	i_dict_mean=i_dict_total/i_dict_n_items
	#median?
	if i_dict_n_items%2==1:
		median_count=int(i_dict_n_items/2)+1
	else:
		median_count=i_dict_n_items/2
	i_median_count=0
	i_median_index=0
	sorted_inp_keys=sorted(list(input_dict.keys()))
	while i_median_count<median_count:
		i_median_count+=input_dict[sorted_inp_keys[i_median_index]]
		i_median_index+=1
	if i_median_count==median_count:
		i_dict_median=sorted_inp_keys[i_median_index]
	elif i_median_index==0:
		i_dict_median=sorted_inp_keys[i_median_index]
	else:
		#check which index is closer to the median
		if abs(median_count-i_median_count)>abs(median_count-i_median_count-input_dict[sorted_inp_keys[i_median_index-1]]):
			i_dict_median=sorted_inp_keys[i_median_index]
		elif abs(median_count-i_median_count)<abs(median_count-i_median_count-input_dict[sorted_inp_keys[i_median_index-1]]):
			i_dict_median=sorted_inp_keys[i_median_index-1]
		else:
			i_dict_median=numpy.mean(sorted_inp_keys[i_median_index-1:i_median_index])
	#SD
	std_sum_sq=0
	for i_std in sorted_inp_keys:
		std_sum_sq+=input_dict[i_std]*((i_std-i_dict_mean)**2)
	i_dict_std=numpy.sqrt((1/(i_dict_n_items))*std_sum_sq)
	#Histogram
	#i_dict_bin_width=float(sorted_inp_keys[-1])/n_bins
	#i_dict_bins=[i_dict_bin_width*i for i in range(n_bins+1)]
	if max_bin==0:
		i_dict_bins=numpy.linspace(0, sorted_inp_keys[-1], n_bins+1)
	else:
		i_dict_bins=numpy.linspace(0, max_bin, n_bins+1)
	i_dict_hist=[0 for c in range(n_bins)]
	curr_bin=0
	#print(len(sorted_inp_keys))
	for i_bin_value in sorted_inp_keys:
		try:
			while i_bin_value>i_dict_bins[curr_bin+1]:
				curr_bin+=1
			i_dict_hist[curr_bin]+=input_dict[i_bin_value]
		except:
			print("ERROR!!!!!!!!!!!!!!!!!!!!!!!!")
			print(i_dict_hist)
			print(i_dict_bins)
			print(i_bin_value)
			print(curr_bin)
			exit(1)
	return(i_dict_mean, i_dict_median, i_dict_std, i_dict_bins, i_dict_hist)

QUAL_not_na_inf=numpy.isfinite(QUAL_array)
#Output functions
with open(output+"_table.tsv", "w") as output_table:
	output_table.write("var_type\tstatistic_category\tstatistic_type\tvar_or_bin_start\tbin_end\tsample\tstatistic\n")
	#write/plot biallelic
	if check_biallelic==True:
		biallelic_plot_time=time.time()
		#get array of biallelic sites and count
		biallelic_sites=n_alt_array==1
		biallelic_counter=sum(biallelic_sites)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "general", "site_count" ,"NA", "NA","all",biallelic_counter))
		#Write/plot QUAL scores
		#1. Subset QUAL scores for biallelic sites, mask inf and count them!
		#QUAL sites
		#count inf!
		biallelic_inf_sites=sum(numpy.logical_and(biallelic_sites, numpy.isinf(QUAL_array)))
		if biallelic_inf_sites>0:
			QUAL_count_dict_bi.update({"Infinite":biallelic_inf_sites})
		biallelic_QUAL_float=numpy.logical_and(QUAL_not_na_inf, biallelic_sites)
		QUAL_count_dict_bi.update({"float":sum(biallelic_QUAL_float)})
		biallelic_QUAL=QUAL_array[biallelic_QUAL_float]
		#get mean, median, std
		bi_QUAL_mean=numpy.mean(biallelic_QUAL)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "general","QUAL_mean", "NA", "NA", "all", bi_QUAL_mean))
		bi_QUAL_median=numpy.median(biallelic_QUAL)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "general","QUAL_median", "NA", "NA", "all",bi_QUAL_median))
		bi_QUAL_std=numpy.std(biallelic_QUAL)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "general", "QUAL_std", "NA", "NA", "all",bi_QUAL_std))
		#get/plot histogram-> write output array of the function to file
		fig,ax=plt.subplots()
		bi_QUAL_hist=plt.hist(biallelic_QUAL, hist_bins, log=True)
		ax.set_title("biallelic QUAL")
		fig.tight_layout()
		plt.savefig(output+"_QUAL_biallelic.pdf")
		#get hist from plot object and write to output tsv n-> counts for each bin (as array), bins -> start and end of each bin (length: n+1)
		for bi_QUAL_hist_bin_i in range(len(bi_QUAL_hist[0])):
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "general", "QUAL_hist",bi_QUAL_hist[1][bi_QUAL_hist_bin_i] ,bi_QUAL_hist[1][bi_QUAL_hist_bin_i+1] , "all",bi_QUAL_hist[0][bi_QUAL_hist_bin_i]))
		#plot other QUAL counts
		fig,ax=plt.subplots()
		bi_qual_plot_cats=sorted(QUAL_count_dict_bi)
		bi_qual_plot_values=[QUAL_count_dict_bi[bi_qual_cat] for bi_qual_cat in bi_qual_plot_cats]
		plt.bar(bi_qual_plot_cats, bi_qual_plot_values, log=True)
		ax.set_title("Biallelic QUAL categories")
		fig.tight_layout()
		plt.savefig(output+"_QUAL_categories_biallelic.pdf")
		for bi_q_cat in bi_qual_plot_cats:
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "general", "QUAL_category", bi_q_cat, "NA", "all",QUAL_count_dict_bi[bi_q_cat] ))
			
		#Write plot INFO stats
		bi_INFO_DP=info_DP_array[biallelic_sites]
		bi_INFO_DP_mean=numpy.mean(bi_INFO_DP)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "INFO","DP_mean", "NA", "NA", "all", bi_INFO_DP_mean))
		bi_INFO_DP_median=numpy.median(bi_INFO_DP)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "INFO","DP_median", "NA", "NA", "all",bi_INFO_DP_median))
		bi_INFO_DP_std=numpy.std(bi_INFO_DP)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "INFO", "DP_std", "NA", "NA", "all",bi_INFO_DP_std))
		#Plot GATK stats
		fig,axs=plt.subplots(7, layout="constrained")
		fig.suptitle("biallelic INFO")
		bi_INFO_DP_hist=axs[0].hist(bi_INFO_DP, hist_bins, log=True)
		axs[0].set_title("biallelic INFO DP")		
		for bi_INFO_DP_hist_bin_i in range(len(bi_INFO_DP_hist[0])):
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "INFO", "DP_hist",bi_INFO_DP_hist[1][bi_INFO_DP_hist_bin_i] ,bi_INFO_DP_hist[1][bi_INFO_DP_hist_bin_i+1] , "all",bi_INFO_DP_hist[0][bi_INFO_DP_hist_bin_i]  ))
		
		bi_INFO_array_GATK=info_array_GATK_bp[biallelic_sites,:]
		for bi_GATK_i in GATK_index_array:
			#select column
			bi_INFO_array_GATK_i=bi_INFO_array_GATK[:,GATK_index_array[bi_GATK_i]]
			#count NAs
			bi_INFO_GATK_i_not_na=numpy.isfinite(bi_INFO_array_GATK_i)
			bi_INFO_GATK_i_not_na_array=bi_INFO_array_GATK_i[bi_INFO_GATK_i_not_na]
			bi_INFO_GATK_i_mean=numpy.mean(bi_INFO_GATK_i_not_na_array)
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "INFO","%s_mean" % bi_GATK_i, "NA", "NA", "all", bi_INFO_GATK_i_mean))
			bi_INFO_GATK_i_median=numpy.median(bi_INFO_GATK_i_not_na_array)
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "INFO","%s_median" % bi_GATK_i, "NA", "NA", "all",bi_INFO_GATK_i_median))
			bi_INFO_GATK_i_std=numpy.std(bi_INFO_GATK_i_not_na_array)
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "INFO", "%s_std" % bi_GATK_i, "NA", "NA", "all",bi_INFO_GATK_i_std))
			bi_GATK_i_hist=axs[GATK_index_array[bi_GATK_i]+1].hist(bi_INFO_GATK_i_not_na_array, hist_bins, log=True)
			axs[GATK_index_array[bi_GATK_i]+1].set_title("biallelic INFO %s, %s observations" % (bi_GATK_i, sum(bi_INFO_GATK_i_not_na)))
			for bi_INFO_GATK_hist_bin_i in range(len(bi_GATK_i_hist[0])):
				output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "INFO", "%s_hist" % bi_GATK_i, bi_GATK_i_hist[1][bi_INFO_GATK_hist_bin_i] ,bi_GATK_i_hist[1][bi_INFO_GATK_hist_bin_i+1] , "all",bi_GATK_i_hist[0][bi_INFO_GATK_hist_bin_i] ))
		plt.gcf().set_size_inches(8, 17.5)
		plt.savefig(output+"_INFO_biallelic.pdf")

		#Write/plot GT stats
		#plot GTs
		fig, ax = plt.subplots(layout="constrained")
		#x = numpy.arange(len(sample_list))
		#get all gts
		#gt_bi_gts=set()
		samples_gts=0
		for gt_bi_sub_dict in gt_bi_GT_dict:
		#	for gt_bi_sub_key in gt_bi_sub_dict:
		#		gt_bi_gts.add(gt_bi_sub_key)
			samples_gts+=1
		width = 1/(samples_gts+1+len(sample_list))  # the width of the bars
		#multiplier = 0
		#iterate over samples
		gt_counter=0
		sample_positions=[]
		for gt_dict in gt_bi_GT_dict:
			#for gt_type in sorted(list(gt_dict.keys()), key=len):
			#iterate over genotypes
			sample_counter=0
			for gt_type in gt_dict:
				output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "GT","GT_count", gt_type, "NA", sample_list[len(sample_positions)], gt_dict[gt_type]))
				offset=width*gt_counter
				rects =ax.bar(offset, gt_dict[gt_type], width, label=gt_type, log=True)
				ax.bar_label(rects, labels=[gt_type], padding=3, fontsize=6, rotation='vertical')
				gt_counter += 1
				sample_counter+=1
			sample_positions.append((gt_counter-(sample_counter/2))*width)
			gt_counter+=1
		ax.set_ylabel('GT count')
		ax.set_title("Biallelic Genotypes")
		ax.set_xticks(sample_positions, sample_list)
		plt.gcf().set_size_inches(len(sample_list)*3, 5)
		plt.savefig(output+"_GT_counts_biallelic.pdf",dpi=200)




#TODO this differently!
		#1. iterate over samples again


# The bin_edges are the same for all of the histograms

#bin_edges = numpy.linspace(0, hist_range[1], histogram_bins + 1)


#heights = np.diff(bin_edges)


#centers = bin_edges[:-1] + heights / 2


#GET DP
		
		bi_GT_sample_results=[analyze_count_dict(i, hist_bins, 0) for i in gt_bi_depth_count_dict]
		#bi_GT_binned_maximums=[numpy.max(m[4]) for m in bi_GT_sample_results]
		#bi_GT_x_location_sum=0
		#bi_GT_x_locations=[0]*n_samples
		#for bi_GT_l_i in range(n_samples):
		#	bi_GT_x_locations[bi_GT_l_i]=bi_GT_binned_maximums[bi_GT_l_i]+bi_GT_x_location_sum
		#	bi_GT_x_location_sum+=bi_GT_binned_maximums[bi_GT_l_i]
		fig, axs = plt.subplots(n_samples, layout="constrained")
		for bi_GT_plot_hist_i in range(n_samples):
			#bins ->  bi_GT_sample_results[bi_GT_plot_hist_i][3]
			#values ->  bi_GT_sample_results[bi_GT_plot_hist_i][4]
			for bi_GT_DP_bin_i in range(len(bi_GT_sample_results[bi_GT_plot_hist_i][4])):
				output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "GT","GT_DP",bi_GT_sample_results[bi_GT_plot_hist_i][3][bi_GT_DP_bin_i] ,bi_GT_sample_results[bi_GT_plot_hist_i][3][bi_GT_DP_bin_i+1] , sample_list[bi_GT_plot_hist_i],bi_GT_sample_results[bi_GT_plot_hist_i][4][bi_GT_DP_bin_i] ))
			#bi_GT_bin_edges = bi_GT_sample_results[bi_GT_plot_hist_i][3]
			bi_GT_DP_bin_width=bi_GT_sample_results[bi_GT_plot_hist_i][3][1]-bi_GT_sample_results[bi_GT_plot_hist_i][3][0]
			bi_GT_DP_hist_i=axs[bi_GT_plot_hist_i].bar(bi_GT_sample_results[bi_GT_plot_hist_i][3][:-1],bi_GT_sample_results[bi_GT_plot_hist_i][4],width=bi_GT_DP_bin_width,align="edge", log=True)
			axs[bi_GT_plot_hist_i].set_title(sample_list[bi_GT_plot_hist_i])		
			#bi_GT_heights = numpy.diff(bi_GT_bin_edges)
			#bi_GT_centers = bi_GT_bin_edges[:-1] + bi_GT_heights / 2 
			#bi_GT_lefts = [bi_GT_x_locations[bi_GT_plot_hist_i] - 0.5 * u for u in bi_GT_sample_results[bi_GT_plot_hist_i][4]]
			#ax.barh(bi_GT_centers, bi_GT_sample_results[bi_GT_plot_hist_i][4], height=bi_GT_heights, left=bi_GT_lefts, color='b')
		#ax.set_xticks(bi_GT_x_locations, sample_list)
		plt.gcf().set_size_inches(8, 2.5*n_samples)
		plt.savefig(output+"_GT_DP_biallelic.pdf")
		#3. Calculate histograms from dictionaries
		#write and plot histograms...
		#Get GQ

		bi_GQ_sample_results=[analyze_count_dict(i, 99, 99) for i in gt_bi_GQ_count_dict]
		fig, axs = plt.subplots(n_samples, layout="constrained")
		for bi_GQ_plot_hist_i in range(n_samples):
			for bi_GQ_DP_bin_i in range(len(bi_GQ_sample_results[bi_GQ_plot_hist_i][4])):
				output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "GT","GT_GQ",bi_GQ_sample_results[bi_GQ_plot_hist_i][3][bi_GQ_DP_bin_i] ,bi_GQ_sample_results[bi_GQ_plot_hist_i][3][bi_GQ_DP_bin_i+1] , sample_list[bi_GQ_plot_hist_i],bi_GQ_sample_results[bi_GQ_plot_hist_i][4][bi_GQ_DP_bin_i] ))
			#bi_GQ_bin_edges = bi_GQ_sample_results[bi_GQ_plot_hist_i][3]
			bi_GQ_DP_bin_width=bi_GQ_sample_results[bi_GQ_plot_hist_i][3][1]-bi_GQ_sample_results[bi_GQ_plot_hist_i][3][0]
			bi_GQ_DP_hist_i=axs[bi_GQ_plot_hist_i].bar(bi_GQ_sample_results[bi_GQ_plot_hist_i][3][:-1],bi_GQ_sample_results[bi_GQ_plot_hist_i][4],width=bi_GQ_DP_bin_width,align="edge", log=True)
			axs[bi_GQ_plot_hist_i].set_title(sample_list[bi_GQ_plot_hist_i])		
		plt.gcf().set_size_inches(8, 2.5*n_samples)
		plt.savefig(output+"_GT_GQ_biallelic.pdf")
		#make violin plot for DP
		#subset DP array
#		bi_GT_depth_array=GT_depth_array[biallelic_sites,:]
#		fig,ax=plt.subplots()
#		ax.violinplot(bi_GT_depth_array, showmeans=True, showmedians=True, showextrema=True, widths=0.7)
#		#ax.boxplot(bi_GT_depth_array)
#		ax.set_title("GT DP")
#		ax.set_xticks(numpy.arange(1, len(sample_list) + 1), labels=sample_list)
#		#ax.set_xlim(0.25, len(sample_list) + 0.75)
#		ax.set_xlabel('Sample')
#		plt.savefig(output+"_GT_DP_biallelic.pdf")


		#make violin plot for GQ
#		bi_GQ_array=GQ_array[biallelic_sites,:]
#		fig,ax=plt.subplots()
#		bi_GQ_mask = ~numpy.isnan(bi_GQ_array)
#		bi_GQ_filtered_data = [d[m] for d, m in zip(bi_GQ_array.T, bi_GQ_mask.T)]
#		ax.violinplot(bi_GQ_filtered_data, showmeans=True, showmedians=True, showextrema=True, widths=0.7)
		#ax.boxplot(bi_GT_depth_array)
#		ax.set_title("GT GQ")
#		ax.set_xticks(numpy.arange(1, len(sample_list) + 1), labels=sample_list)
#		ax.set_xlim(0.25, len(sample_list) + 0.75)
#		ax.set_xlabel('Sample')
#		plt.savefig(output+"_GT_GQ_biallelic.pdf")
		#plot 2d hist of DP and GQ


		#print("Done writing/plotting biallelic sites in %s s" % str(time.time()-biallelic_plot_time))


	#write/plot biallelic
	if check_multiallelic==True:
		multiallelic_plot_time=time.time()
		#get array of multiallelic sites and count
		multiallelic_sites=n_alt_array>1
		multiallelic_counter=sum(multiallelic_sites)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "general", "site_count" ,"NA", "NA","all",multiallelic_counter))
		#Make histogram of number of alt alleles
		ma_n_alt_alleles=n_alt_array[multiallelic_sites]
		ma_n_alt_alleles_mean=numpy.mean(ma_n_alt_alleles)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "general","n_alt_alleles_mean", "NA", "NA", "all", ma_n_alt_alleles_mean))
		ma_n_alt_alleles_median=numpy.median(ma_n_alt_alleles)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "general","n_alt_alleles_median", "NA", "NA", "all",ma_n_alt_alleles_median))
		ma_n_alt_alleles_std=numpy.std(ma_n_alt_alleles)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "general", "n_alt_alleles_std", "NA", "NA", "all",ma_n_alt_alleles_std))
		fig,ax=plt.subplots()
		#TODO: improve plot
		ma_n_alt_alleles_hist=plt.hist(ma_n_alt_alleles, hist_bins, log=True)
		ax.set_title("multiallelic n_alt_alleles")
		fig.tight_layout()
		plt.savefig(output+"_n_alt_alleles_multiallelic.pdf")
		#get hist from plot object and write to output tsv n-> counts for each bin (as array), bins -> start and end of each bin (length: n+1)
		for ma_n_alt_alleles_hist_bin_i in range(len(ma_n_alt_alleles_hist[0])):
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "general", "n_alt_alleles_hist",ma_n_alt_alleles_hist[1][ma_n_alt_alleles_hist_bin_i] ,ma_n_alt_alleles_hist[1][ma_n_alt_alleles_hist_bin_i+1] , "all",ma_n_alt_alleles_hist[0][ma_n_alt_alleles_hist_bin_i]))

		#Write/plot QUAL scores
		#1. Subset QUAL scores for multiallelic sites, mask inf and count them!
		#QUAL sites
		#count inf!
		multiallelic_inf_sites=sum(numpy.logical_and(multiallelic_sites, numpy.isinf(QUAL_array)))
		if multiallelic_inf_sites>0:
			QUAL_count_dict_ma.update({"Infinite":multiallelic_inf_sites})
		multiallelic_QUAL_float=numpy.logical_and(QUAL_not_na_inf, multiallelic_sites)
		QUAL_count_dict_ma.update({"float":sum(multiallelic_QUAL_float)})
		multiallelic_QUAL=QUAL_array[multiallelic_QUAL_float]
		#get mean, median, std
		ma_QUAL_mean=numpy.mean(multiallelic_QUAL)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "general","QUAL_mean", "NA", "NA", "all", ma_QUAL_mean))
		ma_QUAL_median=numpy.median(multiallelic_QUAL)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "general","QUAL_median", "NA", "NA", "all",ma_QUAL_median))
		ma_QUAL_std=numpy.std(multiallelic_QUAL)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "general", "QUAL_std", "NA", "NA", "all",ma_QUAL_std))
		#get/plot histogram-> write output array of the function to file
		fig,ax=plt.subplots()
		ma_QUAL_hist=plt.hist(multiallelic_QUAL, hist_bins, log=True)
		ax.set_title("multiallelic QUAL")
		fig.tight_layout()
		plt.savefig(output+"_QUAL_multiallelic.pdf")
		#get hist from plot object and write to output tsv n-> counts for each bin (as array), bins -> start and end of each bin (length: n+1)
		for ma_QUAL_hist_bin_i in range(len(ma_QUAL_hist[0])):
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "general", "QUAL_hist",ma_QUAL_hist[1][ma_QUAL_hist_bin_i] ,ma_QUAL_hist[1][ma_QUAL_hist_bin_i+1] , "all",ma_QUAL_hist[0][ma_QUAL_hist_bin_i]))
		#plot other QUAL counts
		fig,ax=plt.subplots()
		ma_qual_plot_cats=sorted(QUAL_count_dict_ma)
		ma_qual_plot_values=[QUAL_count_dict_ma[ma_qual_cat] for ma_qual_cat in ma_qual_plot_cats]
		plt.bar(ma_qual_plot_cats, ma_qual_plot_values, log=True)
		ax.set_title("Multiallelic QUAL categories")
		fig.tight_layout()
		plt.savefig(output+"_QUAL_categories_multiallelic.pdf")
		for ma_q_cat in ma_qual_plot_cats:
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "general", "QUAL_category", ma_q_cat, "NA", "all",QUAL_count_dict_ma[ma_q_cat] ))
			
		#Write plot INFO stats
		ma_INFO_DP=info_DP_array[multiallelic_sites]
		ma_INFO_DP_mean=numpy.mean(ma_INFO_DP)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "INFO","DP_mean", "NA", "NA", "all", ma_INFO_DP_mean))
		ma_INFO_DP_median=numpy.median(ma_INFO_DP)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "INFO","DP_median", "NA", "NA", "all",ma_INFO_DP_median))
		ma_INFO_DP_std=numpy.std(ma_INFO_DP)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "INFO", "DP_std", "NA", "NA", "all",ma_INFO_DP_std))
		#Plot GATK stats
		fig,axs=plt.subplots(7, layout="constrained")
		fig.suptitle("multiallelic INFO")
		ma_INFO_DP_hist=axs[0].hist(ma_INFO_DP, hist_bins, log=True)
		axs[0].set_title("multiallelic INFO DP")		
		for ma_INFO_DP_hist_bin_i in range(len(ma_INFO_DP_hist[0])):
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "INFO", "DP_hist",ma_INFO_DP_hist[1][ma_INFO_DP_hist_bin_i] ,ma_INFO_DP_hist[1][ma_INFO_DP_hist_bin_i+1] , "all",ma_INFO_DP_hist[0][ma_INFO_DP_hist_bin_i]  ))
		
		ma_INFO_array_GATK=info_array_GATK_bp[multiallelic_sites,:]
		for ma_GATK_i in GATK_index_array:
			#select column
			ma_INFO_array_GATK_i=ma_INFO_array_GATK[:,GATK_index_array[ma_GATK_i]]
			#count NAs
			ma_INFO_GATK_i_not_na=numpy.isfinite(ma_INFO_array_GATK_i)
			ma_INFO_GATK_i_not_na_array=ma_INFO_array_GATK_i[ma_INFO_GATK_i_not_na]
			ma_INFO_GATK_i_mean=numpy.mean(ma_INFO_GATK_i_not_na_array)
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "INFO","%s_mean" % ma_GATK_i, "NA", "NA", "all", ma_INFO_GATK_i_mean))
			ma_INFO_GATK_i_median=numpy.median(ma_INFO_GATK_i_not_na_array)
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "INFO","%s_median" % ma_GATK_i, "NA", "NA", "all",ma_INFO_GATK_i_median))
			ma_INFO_GATK_i_std=numpy.std(ma_INFO_GATK_i_not_na_array)
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "INFO", "%s_std" % ma_GATK_i, "NA", "NA", "all",ma_INFO_GATK_i_std))
			ma_GATK_i_hist=axs[GATK_index_array[ma_GATK_i]+1].hist(ma_INFO_GATK_i_not_na_array, hist_bins, log=True)
			axs[GATK_index_array[ma_GATK_i]+1].set_title("multiallelic INFO %s, %s observations" % (ma_GATK_i, sum(ma_INFO_GATK_i_not_na)))
			for ma_INFO_GATK_hist_bin_i in range(len(ma_GATK_i_hist[0])):
				output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "INFO", "%s_hist" % ma_GATK_i, ma_GATK_i_hist[1][ma_INFO_GATK_hist_bin_i] ,ma_GATK_i_hist[1][ma_INFO_GATK_hist_bin_i+1] , "all",ma_GATK_i_hist[0][ma_INFO_GATK_hist_bin_i] ))
		plt.gcf().set_size_inches(8, 17.5)
		plt.savefig(output+"_INFO_multiallelic.pdf")

		#Write/plot GT stats
		#plot GTs
		fig, ax = plt.subplots(layout="constrained")
		#x = numpy.arange(len(sample_list))
		#get all gts
		#gt_ma_gts=set()
		samples_gts=0
		for gt_ma_sub_dict in gt_ma_GT_dict:
		#	for gt_ma_sub_key in gt_ma_sub_dict:
		#		gt_ma_gts.add(gt_ma_sub_key)
			samples_gts+=1
		width = 1/(samples_gts+1+len(sample_list))  # the width of the bars
		#multiplier = 0
		#iterate over samples
		gt_counter=0
		sample_positions=[]
		for gt_dict in gt_ma_GT_dict:
			#for gt_type in sorted(list(gt_dict.keys()), key=len):
			#iterate over genotypes
			sample_counter=0
			for gt_type in gt_dict:
				output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "GT","GT_count", gt_type, "NA", sample_list[len(sample_positions)], gt_dict[gt_type]))
				offset=width*gt_counter
				rects =ax.bar(offset, gt_dict[gt_type], width, label=gt_type, log=True)
				ax.bar_label(rects, labels=[gt_type], padding=3, fontsize=6, rotation='vertical')
				gt_counter += 1
				sample_counter+=1
			sample_positions.append((gt_counter-(sample_counter/2))*width)
			gt_counter+=1
		ax.set_ylabel('GT count')
		ax.set_title("Multiallelic Genotypes")
		ax.set_xticks(sample_positions, sample_list)
		plt.gcf().set_size_inches(len(sample_list)*3, 5)
		plt.savefig(output+"_GT_counts_multiallelic.pdf",dpi=200)



		ma_GT_sample_results=[analyze_count_dict(i, hist_bins, 0) for i in gt_ma_depth_count_dict]

		fig, axs = plt.subplots(n_samples, layout="constrained")
		for ma_GT_plot_hist_i in range(n_samples):
			for ma_GT_DP_bin_i in range(len(ma_GT_sample_results[ma_GT_plot_hist_i][4])):
				output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "GT","GT_DP",ma_GT_sample_results[ma_GT_plot_hist_i][3][ma_GT_DP_bin_i] ,ma_GT_sample_results[ma_GT_plot_hist_i][3][ma_GT_DP_bin_i+1] , sample_list[ma_GT_plot_hist_i],ma_GT_sample_results[ma_GT_plot_hist_i][4][ma_GT_DP_bin_i] ))
			ma_GT_DP_bin_width=ma_GT_sample_results[ma_GT_plot_hist_i][3][1]-ma_GT_sample_results[ma_GT_plot_hist_i][3][0]
			ma_GT_DP_hist_i=axs[ma_GT_plot_hist_i].bar(ma_GT_sample_results[ma_GT_plot_hist_i][3][:-1],ma_GT_sample_results[ma_GT_plot_hist_i][4],width=ma_GT_DP_bin_width,align="edge", log=True)
			axs[ma_GT_plot_hist_i].set_title(sample_list[ma_GT_plot_hist_i])		
		plt.gcf().set_size_inches(8, 2.5*n_samples)
		plt.savefig(output+"_GT_DP_multiallelic.pdf")

		ma_GQ_sample_results=[analyze_count_dict(i, 99, 99) for i in gt_ma_GQ_count_dict]
		fig, axs = plt.subplots(n_samples, layout="constrained")
		for ma_GQ_plot_hist_i in range(n_samples):
			for ma_GQ_DP_bin_i in range(len(ma_GQ_sample_results[ma_GQ_plot_hist_i][4])):
				output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "GT","GT_GQ",ma_GQ_sample_results[ma_GQ_plot_hist_i][3][ma_GQ_DP_bin_i] ,ma_GQ_sample_results[ma_GQ_plot_hist_i][3][ma_GQ_DP_bin_i+1] , sample_list[ma_GQ_plot_hist_i],ma_GQ_sample_results[ma_GQ_plot_hist_i][4][ma_GQ_DP_bin_i] ))
			#ma_GQ_bin_edges = ma_GQ_sample_results[ma_GQ_plot_hist_i][3]
			ma_GQ_DP_bin_width=ma_GQ_sample_results[ma_GQ_plot_hist_i][3][1]-ma_GQ_sample_results[ma_GQ_plot_hist_i][3][0]
			ma_GQ_DP_hist_i=axs[ma_GQ_plot_hist_i].bar(ma_GQ_sample_results[ma_GQ_plot_hist_i][3][:-1],ma_GQ_sample_results[ma_GQ_plot_hist_i][4],width=ma_GQ_DP_bin_width,align="edge", log=True)
			axs[ma_GQ_plot_hist_i].set_title(sample_list[ma_GQ_plot_hist_i])		
		plt.gcf().set_size_inches(8, 2.5*n_samples)
		plt.savefig(output+"_GT_GQ_multiallelic.pdf")

		#make violin plot for DP
		#subset DP array
#		ma_GT_depth_array=GT_depth_array[multiallelic_sites,:]
#		fig,ax=plt.subplots()
#		ax.violinplot(ma_GT_depth_array, showmeans=True, showmedians=True, showextrema=True, widths=0.7)
#		#ax.boxplot(ma_GT_depth_array)
#		ax.set_title("GT DP")
#		ax.set_xticks(numpy.arange(1, len(sample_list) + 1), labels=sample_list)
#		#ax.set_xlim(0.25, len(sample_list) + 0.75)
#		ax.set_xlabel('Sample')
#		plt.savefig(output+"_GT_DP_multiallelic.pdf")
#
#
#		#make violin plot for GQ
#		ma_GQ_array=GQ_array[multiallelic_sites,:]
#		fig,ax=plt.subplots()
#		ma_GQ_mask = ~numpy.isnan(ma_GQ_array)
#		ma_GQ_filtered_data = [d[m] for d, m in zip(ma_GQ_array.T, ma_GQ_mask.T)]
#		ax.violinplot(ma_GQ_filtered_data, showmeans=True, showmedians=True, showextrema=True, widths=0.7)
#		#ax.boxplot(ma_GT_depth_array)
#		ax.set_title("GT GQ")
#		ax.set_xticks(numpy.arange(1, len(sample_list) + 1), labels=sample_list)
#		ax.set_xlim(0.25, len(sample_list) + 0.75)
#		ax.set_xlabel('Sample')
#		plt.savefig(output+"_GT_GQ_multiallelic.pdf")
#		#plot 2d hist of DP and GQ
#
		#print("Done writing/plotting multiallelic sites in %s s" % str(time.time()-multiallelic_plot_time))

	if check_invariants==True:
		invariant_plot_time=time.time()
		#get array of invariant sites and count
		invariant_sites=n_alt_array==0
		invariant_sites[line_counter:]=False

		invariant_counter=sum(invariant_sites)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("invariant", "general", "site_count" ,"NA", "NA","all",invariant_counter))
		#Write/plot QUAL scores
		#1. Subset QUAL scores for invariant sites, mask inf and count them!
		#QUAL sites
		#count inf!
		invariant_inf_sites=sum(numpy.logical_and(invariant_sites, numpy.isinf(QUAL_array)))
		if invariant_inf_sites>0:
			QUAL_count_dict_inv.update({"Infinite":invariant_inf_sites})
		invariant_QUAL_float=numpy.logical_and(QUAL_not_na_inf, invariant_sites)
		QUAL_count_dict_inv.update({"float":sum(invariant_QUAL_float)})
		invariant_QUAL=QUAL_array[invariant_QUAL_float]
		#get mean, median, std
		inv_QUAL_mean=numpy.mean(invariant_QUAL)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("invariant", "general","QUAL_mean", "NA", "NA", "all", inv_QUAL_mean))
		inv_QUAL_median=numpy.median(invariant_QUAL)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("invariant", "general","QUAL_median", "NA", "NA", "all",inv_QUAL_median))
		inv_QUAL_std=numpy.std(invariant_QUAL)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("invariant", "general", "QUAL_std", "NA", "NA", "all",inv_QUAL_std))
		#get/plot histogram-> write output array of the function to file
		fig,ax=plt.subplots()
		inv_QUAL_hist=plt.hist(invariant_QUAL, hist_bins, log=True)
		ax.set_title("invariant QUAL")
		fig.tight_layout()
		plt.savefig(output+"_QUAL_invariant.pdf")
		#get hist from plot object and write to output tsv n-> counts for each bin (as array), bins -> start and end of each bin (length: n+1)
		for inv_QUAL_hist_bin_i in range(len(inv_QUAL_hist[0])):
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("invariant", "general", "QUAL_hist",inv_QUAL_hist[1][inv_QUAL_hist_bin_i] ,inv_QUAL_hist[1][inv_QUAL_hist_bin_i+1] , "all",inv_QUAL_hist[0][inv_QUAL_hist_bin_i]))
		#plot other QUAL counts
		fig,ax=plt.subplots()
		inv_qual_plot_cats=sorted(QUAL_count_dict_inv)
		inv_qual_plot_values=[QUAL_count_dict_inv[inv_qual_cat] for inv_qual_cat in inv_qual_plot_cats]
		plt.bar(inv_qual_plot_cats, inv_qual_plot_values, log=True)
		ax.set_title("Multiallelic QUAL categories")
		fig.tight_layout()
		plt.savefig(output+"_QUAL_categories_invariant.pdf")
		for inv_q_cat in inv_qual_plot_cats:
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("invariant", "general", "QUAL_category", inv_q_cat, "NA", "all",QUAL_count_dict_inv[inv_q_cat] ))
			
		#Write plot INFO stats
		inv_INFO_DP=info_DP_array[invariant_sites]
		inv_INFO_DP_mean=numpy.mean(inv_INFO_DP)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("invariant", "INFO","DP_mean", "NA", "NA", "all", inv_INFO_DP_mean))
		inv_INFO_DP_median=numpy.median(inv_INFO_DP)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("invariant", "INFO","DP_median", "NA", "NA", "all",inv_INFO_DP_median))
		inv_INFO_DP_std=numpy.std(inv_INFO_DP)
		output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("invariant", "INFO", "DP_std", "NA", "NA", "all",inv_INFO_DP_std))
		#Plot GATK stats
		fig,ax=plt.subplots()
		inv_INFO_DP_hist=ax.hist(inv_INFO_DP, hist_bins, log=True)
		for inv_INFO_DP_hist_bin_i in range(len(inv_INFO_DP_hist[0])):
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("invariant", "INFO", "DP_hist",inv_INFO_DP_hist[1][inv_INFO_DP_hist_bin_i] ,inv_INFO_DP_hist[1][inv_INFO_DP_hist_bin_i+1] , "all",inv_INFO_DP_hist[0][inv_INFO_DP_hist_bin_i]  ))
		ax.set_title("invariant INFO DP")		
		fig.tight_layout()
		plt.savefig(output+"_INFO_invariant.pdf")

		#Write/plot GT stats
		#plot GTs
		fig, ax = plt.subplots(layout="constrained")
		#x = numpy.arange(len(sample_list))
		#get all gts
		#gt_ma_gts=set()
		samples_gts=0
		for gt_inv_sub_dict in gt_inv_GT_dict:
		#	for gt_inv_sub_key in gt_inv_sub_dict:
		#		gt_inv_gts.add(gt_inv_sub_key)
			samples_gts+=1
		width = 1/(samples_gts+1+len(sample_list))  # the width of the bars
		#multiplier = 0
		#iterate over samples
		gt_counter=0
		sample_positions=[]
		for gt_dict in gt_inv_GT_dict:
			#for gt_type in sorted(list(gt_dict.keys()), key=len):
			#iterate over genotypes
			sample_counter=0
			for gt_type in gt_dict:
				output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("invariant", "GT","GT_count", gt_type, "NA", sample_list[len(sample_positions)], gt_dict[gt_type]))
				offset=width*gt_counter
				rects =ax.bar(offset, gt_dict[gt_type], width, label=gt_type, log=True)
				ax.bar_label(rects, labels=[gt_type], padding=3, fontsize=6, rotation='vertical')
				gt_counter += 1
				sample_counter+=1
			sample_positions.append((gt_counter-(sample_counter/2))*width)
			gt_counter+=1
		ax.set_ylabel('GT count')
		ax.set_title("Invariant Genotypes")
		ax.set_xticks(sample_positions, sample_list)
		plt.gcf().set_size_inches(len(sample_list)*3, 5)
		plt.savefig(output+"_GT_counts_invariant.pdf",dpi=200)

#Plot GT DP dicts

		inv_GT_sample_results=[analyze_count_dict(i, hist_bins, 0) for i in gt_inv_depth_count_dict]

		fig, axs = plt.subplots(n_samples, layout="constrained")
		for inv_GT_plot_hist_i in range(n_samples):
			for inv_GT_DP_bin_i in range(len(inv_GT_sample_results[inv_GT_plot_hist_i][4])):
				output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("invariant", "GT","GT_DP",inv_GT_sample_results[inv_GT_plot_hist_i][3][inv_GT_DP_bin_i] ,inv_GT_sample_results[inv_GT_plot_hist_i][3][inv_GT_DP_bin_i+1] , sample_list[inv_GT_plot_hist_i],inv_GT_sample_results[inv_GT_plot_hist_i][4][inv_GT_DP_bin_i] ))
			inv_GT_DP_bin_width=inv_GT_sample_results[inv_GT_plot_hist_i][3][1]-inv_GT_sample_results[inv_GT_plot_hist_i][3][0]
			inv_GT_DP_hist_i=axs[inv_GT_plot_hist_i].bar(inv_GT_sample_results[inv_GT_plot_hist_i][3][:-1],inv_GT_sample_results[inv_GT_plot_hist_i][4],width=inv_GT_DP_bin_width,align="edge", log=True)
			axs[inv_GT_plot_hist_i].set_title(sample_list[inv_GT_plot_hist_i])		
		plt.gcf().set_size_inches(8, 2.5*n_samples)
		plt.savefig(output+"_GT_DP_invariant.pdf")

#Plot RGQ dicts


		inv_RGQ_sample_results=[analyze_count_dict(i, 99, 99) for i in gt_inv_RGQ_count_dict]
		fig, axs = plt.subplots(n_samples, layout="constrained")
		for inv_RGQ_plot_hist_i in range(n_samples):
			for inv_RGQ_DP_bin_i in range(len(inv_RGQ_sample_results[inv_RGQ_plot_hist_i][4])):
				output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("invariant", "GT","GT_RGQ",inv_RGQ_sample_results[inv_RGQ_plot_hist_i][3][inv_RGQ_DP_bin_i] ,inv_RGQ_sample_results[inv_RGQ_plot_hist_i][3][inv_RGQ_DP_bin_i+1] , sample_list[inv_RGQ_plot_hist_i],inv_RGQ_sample_results[inv_RGQ_plot_hist_i][4][inv_RGQ_DP_bin_i] ))
			#inv_RGQ_bin_edges = inv_RGQ_sample_results[inv_RGQ_plot_hist_i][3]
			inv_RGQ_DP_bin_width=inv_RGQ_sample_results[inv_RGQ_plot_hist_i][3][1]-inv_RGQ_sample_results[inv_RGQ_plot_hist_i][3][0]
			inv_RGQ_DP_hist_i=axs[inv_RGQ_plot_hist_i].bar(inv_RGQ_sample_results[inv_RGQ_plot_hist_i][3][:-1],inv_RGQ_sample_results[inv_RGQ_plot_hist_i][4],width=inv_RGQ_DP_bin_width,align="edge", log=True)
			axs[inv_RGQ_plot_hist_i].set_title(sample_list[inv_RGQ_plot_hist_i])		
		plt.gcf().set_size_inches(8, 2.5*n_samples)
		plt.savefig(output+"_GT_RGQ_invariant.pdf")



		#make violin plot for DP
#		#subset DP array
#		inv_GT_depth_array=GT_depth_array[invariant_sites,:]
#		fig,ax=plt.subplots()
#		ax.violinplot(inv_GT_depth_array, showmeans=True, showmedians=True, showextrema=True, widths=0.7)
#		#ax.boxplot(inv_GT_depth_array)
#		ax.set_title("GT DP")
#		ax.set_xticks(numpy.arange(1, len(sample_list) + 1), labels=sample_list)
#		#ax.set_xlim(0.25, len(sample_list) + 0.75)
#		ax.set_xlabel('Sample')
#		plt.savefig(output+"_GT_DP_invariant.pdf")
#
#
#		#make violin plot for iRGQ
#		inv_RGQ_array=RGQ_array[invariant_sites,:]
#		fig,ax=plt.subplots()
#		inv_RGQ_mask = ~numpy.isnan(inv_RGQ_array)
#		inv_RGQ_filtered_data = [d[m] for d, m in zip(inv_RGQ_array.T, inv_RGQ_mask.T)]
#		ax.violinplot(inv_RGQ_filtered_data, showmeans=True, showmedians=True, showextrema=True, widths=0.7)
#		#ax.boxplot(inv_GT_depth_array)
#		ax.set_title("GT RGQ")
#		ax.set_xticks(numpy.arange(1, len(sample_list) + 1), labels=sample_list)
#		ax.set_xlim(0.25, len(sample_list) + 0.75)
#		ax.set_xlabel('Sample')
#		plt.savefig(output+"_GT_RGQ_invariant.pdf")
#		#plot 2d hist of DP and GQ

		#print("Done writing/plotting invariant sites in %s s" % str(time.time()-invariant_plot_time))

#currently only for biallelic sites
	if plot_filter:
		filter_dict_cats_list=[]
		for filter_dict_cat in filter_dict:
			output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("biallelic", "FILTER", "GATK_COUNT", filter_dict_cat,"NA", "all", filter_dict[filter_dict_cat]))
			#if filter_dict_cat!="PASS":
			#	filter_dict_cats_list=filter_dict_cats_list+filter_dict_cat.split(";")
		unique_filter_dict_cats=["MQRankSum","FS","ReadPosRankSum","SOR","MQ","QD"]
		filter_iter_dict={}
		for filter_comb in range(1,len(unique_filter_dict_cats)+1):#
			filter_iter_dict.update({i:0 for i in itertools.combinations(unique_filter_dict_cats, filter_comb)})
		for filter_dict_cat_again in filter_dict:
			filter_dict_cat_again_set=set(filter_dict_cat_again.split(";"))
			for filter_iter_cat in filter_iter_dict:
				if set(filter_iter_cat)==filter_dict_cat_again_set:
					filter_iter_dict[filter_iter_cat]=filter_dict[filter_dict_cat_again]
		fig, ax=plt.subplots(2,1, sharex=True, figsize=(16,9), gridspec_kw={'height_ratios': [3, 1]})
		filter_n_pass=0
		if "PASS" in filter_dict:
			filter_n_pass=filter_dict["PASS"]
		ax[0].bar(["PASS"]+[",".join(key) for key in filter_iter_dict.keys()], [filter_n_pass]+list(filter_iter_dict.values()), log=True)
		ax[0].tick_params('x', labelbottom=False, bottom=False)
		filter_set_array=[[0]+[255]*len(filter_iter_dict)]
		for filter_unique_cat in unique_filter_dict_cats:
			f_u_c_array=[255]
			for f_u in filter_iter_dict:
				if filter_unique_cat in set(f_u):
					f_u_c_array.append(0)
				else:	
					f_u_c_array.append(255)
			filter_set_array.append(f_u_c_array)
		ax[1].imshow(numpy.array(filter_set_array), cmap='gray', vmin=0, vmax=255)
		ax[1].set_yticks(numpy.arange(len(unique_filter_dict_cats)+1), labels=["PASS"]+unique_filter_dict_cats) 
		ax[1].tick_params('x', labelbottom=False, bottom=False)
		fig.tight_layout(h_pad=-1)
		plt.savefig(output+"_GATK_filters.pdf")
