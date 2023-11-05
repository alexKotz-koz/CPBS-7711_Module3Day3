import pandas as pd
import numpy as np
import random

string = pd.read_csv("STRING 1.txt", delimiter="\t", names=["Gene1", "Gene2", "Weight"])


def fa_genes(file_name):
    FA_genes = {}
    for line in open(file_name, "r"):
        locus = line.split("\t")[0]
        print(f"locus: {locus}")
        genes = line.split("\t")[2:]
        FA_genes[locus] = genes
    return FA_genes


FA_genes = fa_genes("Input.gmt.txt")


def network(FA_genes, string, Gene1, Gene2, Weight):
    FA_genes = set(
        gene for genes in FA_genes.values() for gene in genes
    )  # Flatten list of lists and convert to set
    network = []
    for index, row in string.iterrows():
        if row[Gene1] in FA_genes or row[Gene2] in FA_genes:
            network.append((row[Gene1], row[Gene2], row[Weight]))
    return network


network = network(FA_genes, string, "Gene1", "Gene2", "Weight")


locus = FA_genes.keys()
print(f"locus2: {locus}")


# randomly select one gene per key in FA_genes and their weight from network
def random_gene(network):
    random_network = []
    random.shuffle(network)
    for locus in FA_genes.keys():
        found = False
        for gene in FA_genes[locus]:
            for row in network:
                if gene == row[0] or gene == row[1]:
                    random_network.append(row)
                    found = True
                    break
            if found:
                break

    return random_network


random_network = random_gene(network)
print("rand network", random_network)
