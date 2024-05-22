import argparse
import gzip
import numpy

parser = argparse.ArgumentParser(description='Creates depth mask for loci with depth >mean+2*sd')
parser.add_argument('-v', '--vcf')
parser.add_argument('-o', '--output')
parser.add_argument('-c', '--countout')
parser.add_argument('-p', '--plothist')
parser.add_argument('-n', '--nsamples')
parser.add_argument('-l', '--locinumber')
parser.add_argument('-m', '--min_depth')

args = parser.parse_args()

vcf_input=args.vcf
output=args.output
nsamples=int(round(float(args.nsamples), 0))
nloci=int(args.locinumber)
counts=args.countout
histcounts=args.plothist
min_depth=int(args.min_depth)


if vcf_input[-3:]==".gz":
    inp_file=gzip.open(vcf_input, 'rt')
else:
    inp_file=open(vcf_input)

loci_counter=0
for line in inp_file:
    if line[0]=="#":
        if line[1]!="#":
            sample_cats=line.strip().split("\t")
            sample_length=len(sample_cats[9:])
            sample_ids=sample_cats[9:]
#initialize a large numpy arrays

#            data_array=numpy.zeros((nloci,sample_length),dtype=numpy.uint16)
            data_array = numpy.empty((nloci,sample_length),dtype=numpy.uint16)
         #   data_array.fill(numpy.nan)
          #  name_array=numpy.full(nloci, "", dtype=numpy.object)
            name_array=numpy.full(nloci, "", dtype=object)
    else:
        if loci_counter==nloci:
            print("Warning: --locinumber is smaller than the actual number of loci in the VCF file. Growing internal data array by 100000 lines. This is computationally inefficient. Avoid this by increasing the value for --locinumber")
            add_array=numpy.empty((100000, sample_length),dtype=numpy.uint16)
            add_array.fill(numpy.nan)
            data_array=numpy.vstack((data_array, add_array))
            name_array=numpy.concatenate((name_array, numpy.full(100000, "", dtype=object)))
        loc_cats=line.strip().split("\t")
        name_array[loci_counter]=loc_cats[0]+"\t"+loc_cats[1]
        #Do we have to do this for every locus? Maybe yes to be safe.
        loc_info_dp=loc_cats[8].split(":").index("DP")
        loc_snp_cats=loc_cats[9:]
        for i in range(len(loc_snp_cats)):
            loc_snp_fields=loc_snp_cats[i].split(":")
            if len(loc_snp_fields)>=loc_info_dp+1:
                if loc_snp_fields[loc_info_dp]!=".":
                    data_array[loci_counter, i]=numpy.uint16(loc_snp_fields[loc_info_dp])
                else:
                    data_array[loci_counter, i]=0

        loci_counter+=1

#final_array=data_array[:loci_counter]

data_mean=numpy.empty(len(sample_ids))
data_std=numpy.empty(len(sample_ids))
for sample_column_i in range(len(sample_ids)):
    column_array=data_array[:loci_counter,sample_column_i]
    data_mean[sample_column_i]=numpy.mean(column_array[column_array>=min_depth])
    data_std[sample_column_i]=numpy.std(column_array[column_array>=min_depth])
cutoff=data_mean+(2*data_std)

for sample_it in range(len(sample_ids)):
    print("Sample %s: Mean depth: %s, Standard deviation: %s, Depth cutoff: %s" %(sample_ids[sample_it],data_mean[sample_it], data_std[sample_it], cutoff[sample_it]))

hist_dict={d:0 for d in range(sample_length+1)}

with open(output,"w") as output_file:
    with open(counts, "w") as count_output:
        for comp_locus in range(loci_counter):
            depth_count=(data_array[comp_locus]>cutoff).sum()
            hist_dict[depth_count]+=1
            count_output.write(name_array[comp_locus]+"\t"+str(depth_count)+"\n")
            if depth_count>nsamples:
                output_file.write(name_array[comp_locus]+"\n")
with open(histcounts, "w") as hist_output:
    for hist_count in hist_dict:
        hist_output.write("%s\t%s\n" % (hist_count, hist_dict[hist_count]))
