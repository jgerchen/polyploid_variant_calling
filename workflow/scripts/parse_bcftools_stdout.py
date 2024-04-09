import argparse
import sys
import numpy
import numpy.ma as ma
import matplotlib.pyplot as plt
import time


parser=argparse.ArgumentParser()
parser.add_argument('--n_sites')
parser.add_argument('--output')
parser.add_argument('--histogram_bins')
parser.add_argument('--invariants', action='store_true')
parser.add_argument('--biallelic', action='store_true')
parser.add_argument('--multiallelic', action='store_true')
args=parser.parse_args()

n_sites=int(args.n_sites)
hist_bins=int(args.histogram_bins)
output=args.output
check_invariants=args.invariants
check_biallelic=args.biallelic
check_multiallelic=args.multiallelic
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
			GT_depth_array=numpy.zeros((n_sites, n_samples), dtype=numpy.uint32)
			#Try to put GQ/RGQ in arrays!
			if check_biallelic==True or check_invariants==True:
				GQ_array=numpy.empty((n_sites, n_samples))
				GQ_array.fill(numpy.nan)
			if check_invariants==True:
				RGQ_array=numpy.empty((n_sites, n_samples))
				RGQ_array.fill(numpy.nan)

			if check_biallelic==True:
			#	gt_bi_GQ_dict=[{i_gq:0 for i_gq in range(100)} for i_s in sample_list]
				gt_bi_GT_dict=[{} for i_s in sample_list]
			if check_multiallelic==True:
			#	gt_ma_GQ_dict=[{i_gq:0 for i_gq in range(100)} for i_s in sample_list]
				gt_ma_GT_dict=[{} for i_s in sample_list]
			if check_invariants==True:
			#	gt_inv_RGQ_dict=gt_bi_GQ_dict=[{i_gq:0 for i_gq in range(100)} for i_s in sample_list]
				gt_inv_GT_dict=[{} for i_s in sample_list]

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
				if qual in QUAL_count_dict_inv:
					QUAL_count_dict_inv[qual]+=1
				else:	
					QUAL_count_dict_inv.update({qual:1})
			elif n_alt_alleles==1:
				if qual in QUAL_count_dict_bi:
					QUAL_count_dict_bi[qual]+=1
				else:	
					QUAL_count_dict_bi.update({qual:1})
			else:
				if qual in QUAL_count_dict_ma:
					QUAL_count_dict_ma[qual]+=1
				else:	
					QUAL_count_dict_ma.update({qual:1})

		filter_cat=line_cats[6]
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
		gt_DP_i=gt_info.index("DP")
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
			GT_depth_array[line_counter][gt_i]=gt_depth	
			gt_GT=gt_all[0]
			#test if invariant, bi or ma
			if gt_RGQ_i:
				try:
					RGQ_array[line_counter][gt_i]=float(gt_all[gt_RGQ_i])
				except:
					pass
			if gt_GQ_i:
				try:
					GQ_array[line_counter][gt_i]=float(gt_all[gt_GQ_i])
				except:
					pass
			if n_alt_alleles==0:
				if gt_GT not in gt_inv_GT_dict[gt_i]:
					gt_inv_GT_dict[gt_i].update({gt_GT:1})
				else:
					gt_inv_GT_dict[gt_i][gt_GT]+=1
				#	try:
				#		gt_inv_RGQ_dict[gt_i]=int(gt_all[gt_RGQ_i])
						
				#	except:
				#		pass
			elif n_alt_alleles==1:
				if gt_GT not in gt_bi_GT_dict[gt_i]:
					gt_bi_GT_dict[gt_i].update({gt_GT:1})
				else:
					gt_bi_GT_dict[gt_i][gt_GT]+=1
				#if gt_GQ_i:
					#try:
					#	gt_bi_GQ_dict[gt_i]=int(gt_all[gt_GQ_i])
					#except:
					#	pass
			else:
				if gt_GT not in gt_ma_GT_dict[gt_i]:
					gt_ma_GT_dict[gt_i].update({gt_GT:1})
				else:
					gt_ma_GT_dict[gt_i][gt_GT]+=1
				#if gt_GQ_i:
				#	try:
				#		gt_ma_GQ_dict[gt_i]=int(gt_all[gt_GQ_i])
				#	except:
				#		pass
				
		line_counter+=1
	line = sys.stdin.readline()

#print("Done reading VCF in %s s" % str(time.time()-vcf_read_time))

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


		#make violin plot for DP
		#subset DP array
		bi_GT_depth_array=GT_depth_array[biallelic_sites,:]
		fig,ax=plt.subplots()
		ax.violinplot(bi_GT_depth_array, showmeans=True, showmedians=True, showextrema=True, widths=0.7)
		#ax.boxplot(bi_GT_depth_array)
		ax.set_title("GT DP")
		ax.set_xticks(numpy.arange(1, len(sample_list) + 1), labels=sample_list)
		#ax.set_xlim(0.25, len(sample_list) + 0.75)
		ax.set_xlabel('Sample')
		plt.savefig(output+"_GT_DP_biallelic.pdf")


		#make violin plot for GQ
		bi_GQ_array=GQ_array[biallelic_sites,:]
		fig,ax=plt.subplots()
		bi_GQ_mask = ~numpy.isnan(bi_GQ_array)
		bi_GQ_filtered_data = [d[m] for d, m in zip(bi_GQ_array.T, bi_GQ_mask.T)]
		ax.violinplot(bi_GQ_filtered_data, showmeans=True, showmedians=True, showextrema=True, widths=0.7)
		#ax.boxplot(bi_GT_depth_array)
		ax.set_title("GT GQ")
		ax.set_xticks(numpy.arange(1, len(sample_list) + 1), labels=sample_list)
		ax.set_xlim(0.25, len(sample_list) + 0.75)
		ax.set_xlabel('Sample')
		plt.savefig(output+"_GT_GQ_biallelic.pdf")
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


		#make violin plot for DP
		#subset DP array
		ma_GT_depth_array=GT_depth_array[multiallelic_sites,:]
		fig,ax=plt.subplots()
		ax.violinplot(ma_GT_depth_array, showmeans=True, showmedians=True, showextrema=True, widths=0.7)
		#ax.boxplot(ma_GT_depth_array)
		ax.set_title("GT DP")
		ax.set_xticks(numpy.arange(1, len(sample_list) + 1), labels=sample_list)
		#ax.set_xlim(0.25, len(sample_list) + 0.75)
		ax.set_xlabel('Sample')
		plt.savefig(output+"_GT_DP_multiallelic.pdf")


		#make violin plot for GQ
		ma_GQ_array=GQ_array[multiallelic_sites,:]
		fig,ax=plt.subplots()
		ma_GQ_mask = ~numpy.isnan(ma_GQ_array)
		ma_GQ_filtered_data = [d[m] for d, m in zip(ma_GQ_array.T, ma_GQ_mask.T)]
		ax.violinplot(ma_GQ_filtered_data, showmeans=True, showmedians=True, showextrema=True, widths=0.7)
		#ax.boxplot(ma_GT_depth_array)
		ax.set_title("GT GQ")
		ax.set_xticks(numpy.arange(1, len(sample_list) + 1), labels=sample_list)
		ax.set_xlim(0.25, len(sample_list) + 0.75)
		ax.set_xlabel('Sample')
		plt.savefig(output+"_GT_GQ_multiallelic.pdf")
		#plot 2d hist of DP and GQ

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


		#make violin plot for DP
		#subset DP array
		inv_GT_depth_array=GT_depth_array[invariant_sites,:]
		fig,ax=plt.subplots()
		ax.violinplot(inv_GT_depth_array, showmeans=True, showmedians=True, showextrema=True, widths=0.7)
		#ax.boxplot(inv_GT_depth_array)
		ax.set_title("GT DP")
		ax.set_xticks(numpy.arange(1, len(sample_list) + 1), labels=sample_list)
		#ax.set_xlim(0.25, len(sample_list) + 0.75)
		ax.set_xlabel('Sample')
		plt.savefig(output+"_GT_DP_invariant.pdf")


		#make violin plot for iRGQ
		inv_RGQ_array=RGQ_array[invariant_sites,:]
		fig,ax=plt.subplots()
		inv_RGQ_mask = ~numpy.isnan(inv_RGQ_array)
		inv_RGQ_filtered_data = [d[m] for d, m in zip(inv_RGQ_array.T, inv_RGQ_mask.T)]
		ax.violinplot(inv_RGQ_filtered_data, showmeans=True, showmedians=True, showextrema=True, widths=0.7)
		#ax.boxplot(inv_GT_depth_array)
		ax.set_title("GT RGQ")
		ax.set_xticks(numpy.arange(1, len(sample_list) + 1), labels=sample_list)
		ax.set_xlim(0.25, len(sample_list) + 0.75)
		ax.set_xlabel('Sample')
		plt.savefig(output+"_GT_RGQ_invariant.pdf")
		#plot 2d hist of DP and GQ

		#print("Done writing/plotting invariant sites in %s s" % str(time.time()-invariant_plot_time))






