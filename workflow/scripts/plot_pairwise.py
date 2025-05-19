# takes two vectors - number of variable sites for each gene for each population and the corresponsing number of fixed heterozygotes in these genes 

import matplotlib.pyplot as plt
from itertools import chain

figure_path = "./pairwise.png"

def plot_pairwise(pop1, pop2):
    res1 = [[gene[0] / gene[1] if gene[1] != 0 else 0 for gene in pop1]]
    res2 = [[gene[0] / gene[1] if gene[1] != 0 else 0 for gene in pop2]]
    plt.scatter(res1, res2)
    plt.xlabel('Mean fraction of fixed heterozygotes per gene')
    plt.ylabel('Number of genes')
    plt.savefig(figure_path)
    plt.show()
