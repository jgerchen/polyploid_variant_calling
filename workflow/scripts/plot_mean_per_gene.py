# plots Mean fraction of fixed heterozygotes per gene

import matplotlib.pyplot as plt
from itertools import chain

figure_path = "./mean_per_gene.png"

def plot_mean_per_gene(populations):
    populations_num = len(populations)
    result = [[gene[0] / gene[1] if gene[1] != 0 else 0 for gene in population] for population in populations]
    result = [sum(values) / populations_num for values in zip(*result)]
    plt.hist(result)
    plt.xlabel('Mean fraction of fixed heterozygotes per gene')
    plt.ylabel('Number of genes')
    plt.savefig(figure_path)
    plt.show()