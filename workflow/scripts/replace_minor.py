import sys

nucs={"A", "C", "G", "T"}
#line = sys.stdin.readline()
#while line:
for line in sys.stdin:
	if line[0]=="#":
		#automatically fix header
		out_line=line.strip().replace("ID=PL,Number=G", "ID=PL,Number=.")
		print(out_line)
	else:
		line_cats=line.strip().split("\t")
		ref=line_cats[3]
		alt=line_cats[4]
		alt_alleles=alt.split(",")
		info=line_cats[7]
		gt_cats=line_cats[8]
		gts=line_cats[9:]		
		info_cats={info_cat.split("=")[0]:info_cat.split("=")[1].split(",") for info_cat in info.split(";")}
		if "AC" in info_cats:
			if info_cats["AC"].count("1")==1:
				singleton_index=info_cats["AC"].index("1")
				remove_allele=alt_alleles[singleton_index]
				if remove_allele not in nucs:
					#parse genotypes and replace AC genotype with no-call					
					for gt_i in range(len(gts)):
						gt_cats=gts[gt_i].split(":")
						gt_gt=gt_cats[0].split("/")
						if len(gt_gt)==1:	
							gt_gt=gt_cats[0].split("|")
						if str(singleton_index+1) in gt_cats:
							gt_gt_replace="/".join(["."]*len(gt_gt))
							gt_replace=":".join([gt_gt_replace]+gt_cats[1:])
							gts[gt_i]=gt_replace
							break
					#also change alleles and AC, TODO: also AD?
					new_alt_allele=alt_alleles[singleton_index]
					info_cats["AC"]=[info_cats["AC"][singleton_index]]
					#TODO: turn modified info dict back into string
					info_out=
					print("\t".join(line_cats[0:9]+gts))
