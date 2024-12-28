import argparse
import gzip
import re
parser = argparse.ArgumentParser(description='Generate het mask')
parser.add_argument('-v', '--vcf')
parser.add_argument('-o', '--output')
parser.add_argument('-s', '--samples')
parser.add_argument('-a', '--annotation')
parser.add_argument('-m', '--missing', default='5')

# parser.add_argument('-v', '--vcf', default='alnus.bigt.dp.m.bt.vcf.gz')
# parser.add_argument('-o', '--output', default = 'output')
# parser.add_argument('-s', '--samples', default='alnus_samples.tsv')
# parser.add_argument('-a', '--annotation', default='Alnus_glutinosa-GCA_958979055.1-2024_02-genes.gff3.gz')

#get the library of populations with sample names 
# def get_populations(samples_file):
#     # returns {population: [sample1, sample2 ...]}
#     population_dict = {}
#     with open (samples_file, 'r') as populations:
#         for sample in populations:
#             values = sample.strip().split('\t')
#             # values = sample.strip().split()
#             ploidy = values[-1]
#             if ploidy == '2':
#                 # here i specify pop name explicitly
#                 pop_name = re.sub(r'(\d+)[a-zA-Z]*$', r'\1', values[0])
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
                population = values[0]
                if population in population_dict:
                    population_dict[population].append(sample)
                else:
                    population_dict[population] = [sample]
    return population_dict

def process_vcf(contigs_genes, populations):
    current_contig = None
    previous_contig = None
    indices = None
    current_gene_info = {population: (0, 0) for population in populations} # {population : (total_number, heterozygotes)}
    current_region, prev_region = (-1, -1), (-1, -1)
    pointer = 0 # keeps track of the index of current region in current contig
    need_to_print_info = False

    with (gzip.open(args.vcf, 'rt') if args.vcf.endswith(".gz") else open(args.vcf)) as vcf_file, open(args.output, "w") as output:
        output.write('\t' * 7 + '\t'.join(population for population in populations) + '\n') # printing populations header to output 
        for line in vcf_file:
            if not line.startswith('#'):
                current_contig, position = line.split()[:2]
                position = int(position)
                if current_contig != previous_contig and previous_contig != None:
                    if need_to_print_info:
                        output.write(f"{previous_contig} ({current_region[0]},{current_region[1]})\t" + '\t'.join(f"{info[0]},{info[1]}" for info in current_gene_info.values()) + '\n')
                        need_to_print_info = False
                    pointer = 0
                    current_region, prev_region = (-1, -1), (-1, -1)
                    current_gene_info = {population: (0, 0) for population in populations}
                previous_contig = current_contig


                current_region, pointer = find_gene(current_contig, position, contigs_genes, pointer)
                if current_region == (-1, -1): # if the position is in unannotated region
                    if need_to_print_info:
                        output.write(f"{current_contig} ({prev_region[0]},{prev_region[1]})\t" + '\t'.join(f"{info[0]},{info[1]}" for info in current_gene_info.values()) + '\n')
                        need_to_print_info = False
                        current_gene_info = {population: (0, 0) for population in populations}

                else:
                    if current_region != prev_region and prev_region != (-1, -1):
                        output.write(f"{current_contig} ({prev_region[0]},{prev_region[1]})\t" + '\t'.join(f"{info[0]},{info[1]}" for info in current_gene_info.values()) + '\n')
                        current_gene_info = {population: (0, 0) for population in populations}
                    position_info = get_variant_info(line, indices)
                    current_gene_info = { key: tuple(x + y for x, y in zip(current_gene_info[key], position_info[key])) for key in current_gene_info}
                    need_to_print_info = True
                prev_region = current_region

            elif line.startswith('#') and not line.startswith("##"):
                indices = get_indices(populations, line)
        # printing info about the last gene
        if need_to_print_info:
            output.write(f"{current_contig} ({current_region[0]},{current_region[1]})\t" + '\t'.join(f"{info[0]},{info[1]}" for info in current_gene_info.values()) + '\n')
            pass

# returns 
#   1. the region position is in, (-1, -1) if position is in intergenic region
#   2. pointer to the 
#       1. index of the current    region if the position is in annotated region or
#       2. index of the next valid region if the position is in unannotated region or
#       3. index of the last valid region if the position is after the last annotated region  
def find_gene(contig: str, position: int, contigs_genes: dict, pointer:int) -> tuple[tuple, int]:
    contig_genes = contigs_genes[contig]
    for i in range(len(contig_genes[pointer:])):
        (X, Y) = contig_genes[i + pointer]
        if X <= position <= Y:  # Check if the position is within the range
            return (X, Y), i + pointer
        elif position < X:
            return (-1, -1), i + pointer
    return (-1, -1), len(contig_genes) - 1 # we reached the end of annotated regions 

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

# indices are indices of samples of a certain population 
# line is an entire variant line
def get_variant_info(line, indices):
    # returns {population : (total_number, heterozygotes)} - we need total number in case of missing data  
    # returns (0,0) for population that did't pass missing data threshold
    values = line.split()
    variant_info_per_population = {}
    for population, sample_indices in indices.items():
        sample_infos = [values[sample_index].split(':')[0] for sample_index in sample_indices]
        missing = sum(1 for sample_info in sample_infos if '.' in sample_info)
        hetero = sum(1 for sample_info in sample_infos if '0' in sample_info and '1' in sample_info)
        total = len(sample_infos) - missing
        if missing >= int(args.missing):
            hetero, total = 0, 0
        variant_info_per_population[population] = (hetero, total)
    return variant_info_per_population

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

if __name__ == '__main__':
    args = parser.parse_args()
    populations = get_populations(args.samples)
    for population, samples in populations.items():
        print(f"{population}: [{samples}]")
    contigs_genes = get_genes_from_annotation(args.annotation)
    process_vcf(contigs_genes, populations)
