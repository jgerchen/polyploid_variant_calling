import argparse
import sys
import numpy
import numpy.ma as ma
#import matplotlib.pyplot as plt
import time
import pickle

parser=argparse.ArgumentParser()
parser.add_argument('--n_sites')
parser.add_argument('--output')
parser.add_argument('--invariants', action='store_true')
parser.add_argument('--biallelic', action='store_true')
parser.add_argument('--multiallelic', action='store_true')
args=parser.parse_args()

n_sites=int(args.n_sites)
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

########Pickle stuff!!!###############

pickle_dict={"line_counter":line_counter, "sample_list":sample_list, "n_alt_array":n_alt_array, "QUAL_array":QUAL_array, "info_DP_array":info_DP_array}

if check_biallelic==True or check_multiallelic==True:
	pickle_dict.update({"info_array_GATK_bp":info_array_GATK_bp})
if check_biallelic==True: 
	pickle_dict.update({"QUAL_count_dict_bi":QUAL_count_dict_bi,"gt_bi_GT_dict":gt_bi_GT_dict, "gt_bi_depth_count_dict":gt_bi_depth_count_dict, "gt_bi_GQ_count_dict":gt_bi_GQ_count_dict})
if check_invariants==True:
	pickle_dict.update({"QUAL_count_dict_inv":QUAL_count_dict_inv,"gt_inv_GT_dict":gt_inv_GT_dict, "gt_inv_depth_count_dict":gt_inv_depth_count_dict, "gt_inv_RGQ_count_dict":gt_inv_RGQ_count_dict})
if check_multiallelic==True:
#Do we still count n alleles?
	pickle_dict.update({"QUAL_count_dict_ma":QUAL_count_dict_ma,"gt_ma_GT_dict":gt_ma_GT_dict, "gt_ma_depth_count_dict":gt_ma_depth_count_dict, "gt_ma_GQ_count_dict":gt_ma_GQ_count_dict})

with open(output+".pickle", "wb") as f:
	pickle.dump(pickle_dict, f, pickle.HIGHEST_PROTOCOL)


