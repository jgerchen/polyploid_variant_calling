# plots number of fixed heterozygotes over total number of variable sites per each gene per population

import matplotlib.pyplot as plt
from itertools import chain

figure_path = "./per_gene_per_pop.png"

#input in the format lists of lists (genes of a population) of lists ([fixed heterozygotes, number of variable sites])
def plot_per_gene_per_pop(populations):
    flattened = chain.from_iterable(populations)
    var_sites, fixed_hetero = zip(*flattened)
    plt.scatter(var_sites, fixed_hetero)
    plt.xlabel('Number of variable sites')
    plt.ylabel('Number of fixed heterozygotes')
    plt.savefig(figure_path)
    plt.show()
    