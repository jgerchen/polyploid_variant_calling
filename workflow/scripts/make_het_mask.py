# script for counting fixed heterozygotes and total number of variable sites in given genes
import argparse
import gzip
import re

# get the library of populations with sample names
# assuming output in format sample ... ploidy and that population name is not explicitly given
# def get_populations(samples_file):
#     # returns {population: [sample1, sample2 ...]}
#     population_dict = {}
#     with open (samples_file, 'r') as populations:
#         for sample in populations:
#             values = sample.strip().split('\t')
#             ploidy = values[-1]
#             if ploidy == '2':
#                 # here i specify pop name explicitly
#                 pop_name = re.sub(r'(\d+)[a-zA-Z]*$', r'\1', values[0]) # AG001g -> AG001
#                 # pop_name = values[0][:1]
#                 sample_name = values[0]
#                 if pop_name in population_dict:
#                     population_dict[pop_name].append(sample_name)
#                 else:
#                     population_dict[pop_name] = [sample_name]
#     return population_dict

# if the input file is in the format: "sample population ploidy" (so population name is given)
def get_populations(samples_file):
    population_dict = {}
    with open (samples_file, 'r') as populations:
        for sample in populations:
            values = sample.strip().split('\t')
            sample = values[0]
            population = values[1]
            ploidy = values[2]
            if ploidy == '2':
                population = values[1]
                if population in population_dict:
                    population_dict[population].append(sample)
                else:
                    population_dict[population] = [sample]
    return population_dict

# returns contigs = {contig: [(gene1_start, gene1_end), (gene2_start, geen2_end) ...]} 
def get_genes_from_annotation(annotation_file):
    contigs = {}
    with (gzip.open(annotation_file, 'rt') if annotation_file.endswith(".gz") else open(annotation_file)) as annotation:
        current_contig = ""
        genes_of_current_contig = []
        for line in annotation:
            if is_gene_line(line):
                (gene_start, gene_end) = get_gene_range(line)
                genes_of_current_contig.append((gene_start, gene_end))
            elif is_region_line(line):
                if current_contig != "":
                    contigs[current_contig] = genes_of_current_contig
                current_contig = get_contig_alias(line.split()[-1])
                genes_of_current_contig = []
        if current_contig != "":
            contigs[current_contig] = genes_of_current_contig
    return contigs
            
def get_contig_alias(annotation_line):
    return annotation_line.split()[-1].split(';')[1].split('=')[1].split(',')[0] #1	dhAlnGlut1.1	region	1	53352176	.	.	.	ID=region:1;Alias=OY340898.1,NC_084886.1 -> OY340898.1

def get_gene_range(annotation_line):
    return (int(annotation_line.split()[3]), int(annotation_line.split()[4]))#1       ensembl gene    43621   46779   .       +       .       ID=gene:ENSGUIG00005000116;biotype=protein_coding;gene_id=ENSGUIG00005000116;version=1 -> (43621, 46779)

def is_gene_line(annotation_line):   #1	ensembl	ncRNA_gene	12883	15538	.	-	.	ID=gene:ENSGUIG00005000383;biotype=lncRNA;gene_id=ENSGUIG00005000383;version=1 -> true; 
    values = annotation_line.split()
    try:
        return len(values) == 9 and values[2] == "gene"
    except (IndexError, ValueError):
        return False

def is_region_line(annotation_line): #1	dhAlnGlut1.1	region	1	53352176	.	.	.	ID=region:1;Alias=OY340898.1,NC_084886.1 -> true
    values = annotation_line.split()
    try:
        return len(values) > 8 and values[-1].split(':')[0].split('=')[1] == "region"
    except (IndexError, ValueError):
        return False
 
def process_vcf(contigs_genes, populations, args):
    with open ("contig_genes", "w") as file:
        for contig_gene in contigs_genes:
            file.write((f"{contig_gene}: {contigs_genes[contig_gene]}\n\n\n"))
    current_contig = None
    previous_contig = None
    indices = None
    current_gene_info = {population: (0, 0) for population in populations} # {population : (total_number, heterozygotes)}
    current_region, prev_region = (-1, -1), (-1, -1)
    pointer = 0 # keeps track of the index of current region in current contig
    need_to_print_info = False

    with (gzip.open(args.vcf, 'rt') if args.vcf.endswith(".gz") else open(args.vcf)) as vcf_file, open(args.output, "w") as output:
        output.write("contig\tstart\tend\t"+ '\t'.join(population for population in populations) + '\n') # printing populations header to output 
        for line in vcf_file:
            if not line.startswith('#'):
                current_contig, position = line.split()[:2]
                position = int(position)
                if current_contig != previous_contig and previous_contig != None: # if we moved to a next contig
                    if need_to_print_info:
                        output.write(print_gene_info(previous_contig, prev_region, current_gene_info))
                        need_to_print_info = False
                    pointer = 0
                    current_region, prev_region = (-1, -1), (-1, -1)
                    current_gene_info = {population: (0, 0) for population in populations}
                    previous_contig = current_contig
                elif previous_contig == None: # if we moved to the first actual contig 
                    previous_contig = current_contig

                current_region, pointer = find_gene(current_contig, position, contigs_genes, pointer)
                if current_region == (-1, -1): # if the position is in unannotated region
                    if need_to_print_info:
                        output.write(print_gene_info(current_contig, prev_region, current_gene_info))
                        need_to_print_info = False
                        current_gene_info = {population: (0, 0) for population in populations}
                else:
                    if current_region != prev_region and prev_region != (-1, -1):
                        output.write(print_gene_info(current_contig, prev_region, current_gene_info))
                        current_gene_info = {population: (0, 0) for population in populations}
                    position_info = get_variant_info(line, indices, args.missing)
                    current_gene_info = { key: tuple(x + y for x, y in zip(current_gene_info[key], position_info[key])) for key in current_gene_info}
                    need_to_print_info = True
                prev_region = current_region

            elif line.startswith('#') and not line.startswith("##"):
                indices = get_indices(populations, line)
        # printing info about the last gene
        if need_to_print_info:
            output.write(print_gene_info(current_contig, prev_region, current_gene_info))

#returns list: population: [index of sample1 in header, index of sample2 in header ...] from vcf file
def get_indices(populations, header_line):
    indices = {}
    values = header_line.split()
    # population: {pop_name: [index1, index2...]}
    for population, samples in populations.items():
        indices[population] = []
        for sample in samples:
            try:
                indices[population].append(values.index(sample))
            except ValueError:
                raise Exception(f"sample {sample} is not found in vcf header")
    return indices

def print_gene_info(contig, gene_range, gene_info):
    return f"{contig}\t{gene_range[0]}\t{gene_range[1]}\t" + \
                '\t'.join(f"{info[0]},{info[1]}" for info in gene_info.values()) + \
                '\n'

# returns 
#   1. the region position is in, (-1, -1) if position is in intergenic region or if a contig is absent
#   2. pointer to the 
#       1. index of the current    region if the position is in annotated region or
#       2. index of the next valid region if the position is in unannotated region or
#       3. index of the last valid region if the position is after the last annotated region  
def find_gene(contig: str, position: int, contigs_genes: dict, pointer:int) -> tuple[tuple, int]:
    if contig in contigs_genes:
        contig_genes = contigs_genes[contig]
    else:
        return (-1, -1), 0
    for i in range(len(contig_genes[pointer:])):
        (X, Y) = contig_genes[i + pointer]
        if X <= position <= Y:  # Check if the position is within the range
            return (X, Y), i + pointer
        elif position < X:
            return (-1, -1), i + pointer
    return (-1, -1), len(contig_genes) - 1 # we reached the end of annotated regions 



# returns {population : (total_number, heterozygotes)} - we need total number in case of missing data  
# returns (0,0) for population that did't pass missing data threshold
# indices is a dictionary with indices of samples of a certain population 
# line is an entire variant line
def get_variant_info(line, indices, missing_threshold):
    values = line.split()
    variant_info_per_population = {}
    for population, sample_indices in indices.items():
        sample_infos = [values[sample_index].split(':')[0] for sample_index in sample_indices]
        missing = sum(1 for sample_info in sample_infos if '.' in sample_info)       
        if missing/len(sample_infos) >= float(missing_threshold):
            variant_info_per_population[population] = (0, 0)
        else:            
            variant_info_per_population[population] = (int(is_fixed_hetero(sample_infos)), int(is_variable_site(sample_infos)))
    return variant_info_per_population

 filtering_script
def is_fixed_hetero(sample_infos):
    return all(s == "0/1" or s == "0|1" for s in sample_infos if "." not in s) and any ('.' not in s for s in sample_infos)

def is_variable_site(sample_infos):
    return any("0" in s for s in sample_infos) and any("1" in s for s in sample_infos) 

def main():
    parser = argparse.ArgumentParser(description='Generate het mask')
    parser.add_argument('-v', '--vcf')
    parser.add_argument('-o', '--output')
    parser.add_argument('-s', '--samples')
    parser.add_argument('-a', '--annotation')
    parser.add_argument('-m', '--missing', default='0.33')
    args = parser.parse_args()
    populations = get_populations(args.samples)
    contigs_genes = get_genes_from_annotation(args.annotation)
    process_vcf(contigs_genes, populations, args)
    
if __name__ == "__main__":
    main()
def get_contig_alias(annotation_line):
    #1	dhAlnGlut1.1	region	1	53352176	.	.	.	ID=region:1;Alias=OY340898.1,NC_084886.1 -> OY340898.1
    return annotation_line.split()[-1].split(';')[1].split('=')[1].split(',')[0]

def get_contig_id(annotation_line):
    #1	dhAlnGlut1.1	region	1	53352176	.	.	.	ID=region:1;Alias=OY340898.1,NC_084886.1 -> 1
    return annotation_line.split()[-1].split(';')[0].split(':')[1]

def get_gene_id(annotation_line):
    #1	ensembl	ncRNA_gene	12883	15538	.	-	.	ID=gene:ENSGUIG00005000383;biotype=lncRNA;gene_id=ENSGUIG00005000383;version=1 -> ENSGUIG00005000383
    return annotation_line.split()[-1].split(';')[0].split(':')[1]

def get_gene_range(annotation_line):
    ##1	ensembl	ncRNA_gene	12883	15538	.	-	.	ID=gene:ENSGUIG00005000383;biotype=lncRNA;gene_id=ENSGUIG00005000383;version=1 -> (12883, 15538)
    return (int(annotation_line.split()[3]), int(annotation_line.split()[4]))


# returns contigs = {contig: [(gene1_start, gene1_end), (gene2_start, geen2_end) ...]} 
def get_genes_from_annotation(annotation_file):
    contigs = {}
    with (gzip.open(annotation_file, 'rt') if annotation_file.endswith(".gz") else open(annotation_file)) as annotation:
        current_contig = ""
        genes_of_current_contig = []
        for line in annotation:
            if is_gene_line(line):
                (gene_start, gene_end) = get_gene_range(line)
                genes_of_current_contig.append((gene_start, gene_end))
            elif is_region_line(line):
                if current_contig != "":
                    contigs[current_contig] = genes_of_current_contig
                current_contig = get_contig_alias(line.split()[-1])
                genes_of_current_contig = []
        if current_contig != "":
            contigs[current_contig] = genes_of_current_contig
    return contigs
            

def is_gene_line(annotation_line):
    #1	ensembl	ncRNA_gene	12883	15538	.	-	.	ID=gene:ENSGUIG00005000383;biotype=lncRNA;gene_id=ENSGUIG00005000383;version=1 -> true
    values = annotation_line.split()
    try:
        return len(values) > 8 and values[-1].split(':')[0].split('=')[1] == "gene"
    except (IndexError, ValueError):
        return False

def is_region_line(annotation_line):
    #1	dhAlnGlut1.1	region	1	53352176	.	.	.	ID=region:1;Alias=OY340898.1,NC_084886.1 -> true
    values = annotation_line.split()
    try:
        return len(values) > 8 and values[-1].split(':')[0].split('=')[1] == "region"
    except (IndexError, ValueError):
        return False


args = parser.parse_args()
populations = get_populations(args.samples)
contigs_genes = get_genes_from_annotation(args.annotation)
process_vcf(contigs_genes, populations)
main
