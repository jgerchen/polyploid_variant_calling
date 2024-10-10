import argparse
import gzip

parser = argparse.ArgumentParser(description='Generate het mask')
parser.add_argument('-v', '--vcf')
parser.add_argument('-o', '--output')

args = parser.parse_args()

vcf_input=args.vcf
output=args.output

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
    else:
        loc_cats=line.strip().split("\t")
        loc_snp_cats=loc_cats[9:]
