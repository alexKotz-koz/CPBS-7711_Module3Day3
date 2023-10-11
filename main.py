import random
import json
from components.bins import Bins
from components.fagenes import FaGenes
from components.nonfagenes import NonFaGenes
from components.stage1_subnetworks import Stage1_SubNetworks


def find_bin(gene, bins):
    binToReturn = []
    for bin in bins:
        if gene in bins[bin]:
            binToReturn = bins[bin]

    return binToReturn


# Mating or 1:1 (fa to nonfa) replacement
def create_secondary_subnetwork(parentSubnetworks, nonfaGenes, bins):
    newSubnetwork = {}

    newRandomSubnetListToWrite = []

    nonfaGenes = set(nonfaGenes.keys())

    # for each subnetwork in first round of subnetworks
    keyForNewDictionary = 0

    for subnet in parentSubnetworks:
        print(subnet["subnet"].values())
        gene = ""
        gene2 = ""
        newRandomSubnetList = []
        # for each gene in subnetwork
        for item in parentSubnetworks[subnet]:
            newGene = ""
            newGene2 = ""
            if isinstance(item, list):
                gene = item[0]
                gene2 = item[1]
                # print(f"gene: {gene} | gene2: {gene2}")

                geneBin = find_bin(gene, bins)
                geneBin2 = find_bin(gene2, bins)
                newGene = geneBin[random.randrange(0, len(geneBin))]
                newGene2 = geneBin2[random.randrange(0, len(geneBin2))]

                if newGene in nonfaGenes:
                    newRandomSubnetList.append(newGene)

                if newGene2 in nonfaGenes:
                    newRandomSubnetList.append(newGene2)

            elif isinstance(item, str):
                geneBin = find_bin(item, bins)

                if len(geneBin) == 0:
                    geneBin = [item]

                newGene = geneBin[random.randrange(0, len(geneBin))]

                if newGene in nonfaGenes:
                    newRandomSubnetList.append(newGene)

        newSubnetwork.update({keyForNewDictionary: [newRandomSubnetList]})

        keyForNewDictionary += 1

    with open("stage2_random_subnetworks.json", "w") as outputFile:
        json.dump(newSubnetwork, outputFile)


def main():
    ##TEST TO SEE HOW MANY ROWS ARE IN NONFAGENES
    """with open("STRING 1.txt", "r") as file:
    se = set()
    results = [row.split("\t")[:2] for row in file]
    for i in results:
        for j in i:
            se.add(j)
    print(len(se))"""

    binsInstance = Bins("STRING 1.txt")
    bins = binsInstance.create_bins()

    faGenesInstance = FaGenes("STRING 1.txt")
    faGenes = faGenesInstance.fanconi_anemia_genes()

    nonfaGenesInstance = NonFaGenes("STRING 1.txt", faGenes=faGenes)
    nonfaGenes = nonfaGenesInstance.extract_nonfa_genes()

    parentSubnetworksInstance = Stage1_SubNetworks(
        "results.txt", "Input.gmt.txt", "STRING 1.txt"
    )
    parentSubnetworks = parentSubnetworksInstance.create_random_subnetworks()

    """create_secondary_subnetwork(
        parentSubnetworks=parentSubnetworks, nonfaGenes=nonfaGenes, bins=bins
    )"""


if __name__ == "__main__":
    main()
