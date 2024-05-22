import argparse
import sys
import numpy
import numpy.ma as ma
import matplotlib.pyplot as plt
import time
import pickle

parser=argparse.ArgumentParser()
parser.add_argument('--output')
parser.add_argument('--histogram_bins')
parser.add_argument('--invariants', action='store_true')
parser.add_argument('--biallelic', action='store_true')
parser.add_argument('--multiallelic', action='store_true')
parser.add_argument('--pickle_list')


args=parser.parse_args()

hist_bins=int(args.histogram_bins)
output=args.output
check_invariants=args.invariants
check_biallelic=args.biallelic
check_multiallelic=args.multiallelic
pickle_list=args.pickle_list

line_counter=0
sample_list=[]
n_alt_array=numpy.zeros(0, dtype=numpy.uint16)
QUAL_array=numpy.empty(0)
if check_biallelic==True:
	QUAL_count_dict_bi={}
	gt_bi_GT_dict=[]
	gt_bi_depth_count_dict=[]
	gt_bi_GQ_count_dict=[]


if check_invariants==True:
	QUAL_count_dict_inv={}
	gt_inv_GT_dict=[]
	gt_inv_depth_count_dict=[]
	gt_inv_RGQ_count_dict=[]
if check_multiallelic==True:
	QUAL_count_dict_ma={}
	gt_ma_GT_dict=[]
	gt_ma_depth_count_dict=[]
	gt_ma_GQ_count_dict=[]
if check_biallelic==True or check_multiallelic==True:
	info_array_GATK_bp=numpy.empty((0,6))
	GATK_index_array={"FS":0, "QD":1, "MQ":2, "MQRankSum":3, "ReadPosRankSum":4, "SOR":5}



info_DP_array=numpy.empty(0, dtype=numpy.uint32)


#load pickles!
with open(pickle_list) as pickle_input:
	for pickle_line in pickle_input:
		with open(pickle_line.strip(), 'rb') as f:
			pickle_data = pickle.load(f)
			line_counter+=pickle_data["line_counter"]
			sample_list=pickle_data["sample_list"]
			n_samples=len(sample_list)
			n_alt_array=numpy.concatenate((n_alt_array, pickle_data["n_alt_array"]))
			QUAL_array=numpy.concatenate((QUAL_array, pickle_data["QUAL_array"]))
			info_DP_array=numpy.concatenate((info_DP_array, pickle_data["info_DP_array"]))
			if check_biallelic==True or check_multiallelic==True:
				info_array_GATK_bp=numpy.concatenate((info_array_GATK_bp, pickle_data["info_array_GATK_bp"]), axis=0)
			if check_biallelic==True:
				QUAL_count_dict_bi={x: QUAL_count_dict_bi.get(x, 0) + pickle_data["QUAL_count_dict_bi"].get(x, 0) for x in set(QUAL_count_dict_bi).union(pickle_data["QUAL_count_dict_bi"])}
				if len(gt_bi_GT_dict)==0:
					#Check if arrays are empty, if yes just copy from pickle
					gt_bi_GT_dict=pickle_data["gt_bi_GT_dict"]
					gt_bi_depth_count_dict=pickle_data["gt_bi_depth_count_dict"]				
					gt_bi_GQ_count_dict=pickle_data["gt_bi_GQ_count_dict"]				
				else:
					#Otherwise iterate over samples and merge dicts
					for sample_i in range(len(sample_list)):
						gt_bi_GT_dict[sample_i]={x: gt_bi_GT_dict[sample_i].get(x, 0) + pickle_data["gt_bi_GT_dict"][sample_i].get(x, 0) for x in set(gt_bi_GT_dict[sample_i]).union(pickle_data["gt_bi_GT_dict"][sample_i])}
						gt_bi_depth_count_dict[sample_i]={x: gt_bi_depth_count_dict[sample_i].get(x, 0) + pickle_data["gt_bi_depth_count_dict"][sample_i].get(x, 0) for x in set(gt_bi_depth_count_dict[sample_i]).union(pickle_data["gt_bi_depth_count_dict"][sample_i])}
						gt_bi_GQ_count_dict[sample_i]={x: gt_bi_GQ_count_dict[sample_i].get(x, 0) + pickle_data["gt_bi_GQ_count_dict"][sample_i].get(x, 0) for x in set(gt_bi_GQ_count_dict[sample_i]).union(pickle_data["gt_bi_GQ_count_dict"][sample_i])}


			if check_invariants==True:
				QUAL_count_dict_inv={x: QUAL_count_dict_inv.get(x, 0) + pickle_data["QUAL_count_dict_inv"].get(x, 0) for x in set(QUAL_count_dict_inv).union(pickle_data["QUAL_count_dict_inv"])}
				if len(gt_inv_GT_dict)==0:
					#Check if arrays are empty, if yes just copy from pickle
					gt_inv_GT_dict=pickle_data["gt_inv_GT_dict"]
					gt_inv_depth_count_dict=pickle_data["gt_inv_depth_count_dict"]				
					gt_inv_RGQ_count_dict=pickle_data["gt_inv_RGQ_count_dict"]				
				else:
					#Otherwise iterate over samples and merge dicts
					for sample_i in range(len(sample_list)):
						gt_inv_GT_dict[sample_i]={x: gt_inv_GT_dict[sample_i].get(x, 0) + pickle_data["gt_inv_GT_dict"][sample_i].get(x, 0) for x in set(gt_inv_GT_dict[sample_i]).union(pickle_data["gt_inv_GT_dict"][sample_i])}
						gt_inv_depth_count_dict[sample_i]={x: gt_inv_depth_count_dict[sample_i].get(x, 0) + pickle_data["gt_inv_depth_count_dict"][sample_i].get(x, 0) for x in set(gt_inv_depth_count_dict[sample_i]).union(pickle_data["gt_inv_depth_count_dict"][sample_i])}
						gt_inv_RGQ_count_dict[sample_i]={x: gt_inv_RGQ_count_dict[sample_i].get(x, 0) + pickle_data["gt_inv_RGQ_count_dict"][sample_i].get(x, 0) for x in set(gt_inv_RGQ_count_dict[sample_i]).union(pickle_data["gt_inv_RGQ_count_dict"][sample_i])}
			if check_multiallelic==True:
				QUAL_count_dict_ma={x: QUAL_count_dict_ma.get(x, 0) + pickle_data["QUAL_count_dict_ma"].get(x, 0) for x in set(QUAL_count_dict_ma).union(pickle_data["QUAL_count_dict_ma"])}
				if len(gt_ma_GT_dict)==0:
					#Check if arrays are empty, if yes just copy from pickle
					gt_ma_GT_dict=pickle_data["gt_ma_GT_dict"]
					gt_ma_depth_count_dict=pickle_data["gt_ma_depth_count_dict"]				
					gt_ma_GQ_count_dict=pickle_data["gt_ma_GQ_count_dict"]				
				else:
					#Otherwise iterate over samples and merge dicts
					for sample_i in range(len(sample_list)):
						gt_ma_GT_dict[sample_i]={x: gt_ma_GT_dict[sample_i].get(x, 0) + pickle_data["gt_ma_GT_dict"][sample_i].get(x, 0) for x in set(gt_ma_GT_dict[sample_i]).union(pickle_data["gt_ma_GT_dict"][sample_i])}
						gt_ma_depth_count_dict[sample_i]={x: gt_ma_depth_count_dict[sample_i].get(x, 0) + pickle_data["gt_ma_depth_count_dict"][sample_i].get(x, 0) for x in set(gt_ma_depth_count_dict[sample_i]).union(pickle_data["gt_ma_depth_count_dict"][sample_i])}
						gt_ma_GQ_count_dict[sample_i]={x: gt_ma_GQ_count_dict[sample_i].get(x, 0) + pickle_data["gt_ma_GQ_count_dict"][sample_i].get(x, 0) for x in set(gt_ma_GQ_count_dict[sample_i]).union(pickle_data["gt_ma_GQ_count_dict"][sample_i])}
		print("Done reading pickle %s" % pickle_line)


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
#		fig, ax = plt.subplots(layout="constrained")
		#x = numpy.arange(len(sample_list))
		#get all gts
		#gt_ma_gts=set()
#		samples_gts=0
#		for gt_ma_sub_dict in gt_ma_GT_dict:
		#	for gt_ma_sub_key in gt_ma_sub_dict:
		#		gt_ma_gts.add(gt_ma_sub_key)
#			samples_gts+=1
#		width = 1/(samples_gts+1+len(sample_list))  # the width of the bars
		#multiplier = 0
		#iterate over samples
#		gt_counter=0
#		sample_positions=[]
		for gt_dict_i in range(len(gt_ma_GT_dict)):
			#for gt_type in sorted(list(gt_dict.keys()), key=len):
			#iterate over genotypes
#			sample_counter=0
			for gt_type in gt_ma_GT_dict[gt_dict_i]:
				output_table.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("multiallelic", "GT","GT_count", gt_type, "NA", sample_list[gt_dict_i], gt_ma_GT_dict[gt_dict_i][gt_type]))
#				offset=width*gt_counter
#				rects =ax.bar(offset, gt_dict[gt_type], width, label=gt_type, log=True)
#				ax.bar_label(rects, labels=[gt_type], padding=3, fontsize=6, rotation='vertical')
#				gt_counter += 1
#				sample_counter+=1
#			sample_positions.append((gt_counter-(sample_counter/2))*width)
#			gt_counter+=1
#		ax.set_ylabel('GT count')
#		ax.set_title("Multiallelic Genotypes")
#		ax.set_xticks(sample_positions, sample_list)
#		plt.gcf().set_size_inches(len(sample_list)*3, 5)
#		plt.savefig(output+"_GT_counts_multiallelic.pdf",dpi=200)



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





