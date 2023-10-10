import random
import json
from components.bins import Bins
from components.fagenes import FaGenes
from components.nonfagenes import NonFaGenes
from components.stage1_subnetworks import Stage1_SubNetworks


# Mating or 1:1 (fa to nonfa) replacement
def create_secondary_subnetwork(parentSubnetwork, nonfaGenes, bins):
    bins = bins.create_bins()
    parentSubnetwork = parentSubnetwork.create_random_subnetworks()
    newSubnetwork = {}
    newRandomSubnetList = []
    newRandomSubnetListToWrite = []

    # for each subnetwork in first round of subnetworks
    for subnet in parentSubnetwork:
        gene = ""
        gene2 = ""

        # for each gene in subnetwork
        for item in parentSubnetwork[subnet]:
            newGene = ""
            # print(bins)
            if isinstance(item, list):
                gene = item[0]
                gene2 = item[1]
                # print(f"gene: {gene} | gene2: {gene2}")

                # gene = first element in the sublist of the subnetwork
                # for each bin in the bins object
                """for bin in bins:
                    if gene in bins[bin]:
                        if gene in nonfaGenes.keys():
                            newRandomSubnetList.append(gene)
                            print(
                                f"Gene from parent: {gene} | bin: {bins}|{bins[bin]}\n"
                            )
                            # print(newRandomSubnetList)

                        # print(bins[bin][random.randrange(0, len(bins[bin]))])
                        newGene = bins[bin][random.randrange(0, len(bins[bin]))]
                        if newGene in nonfaGenes.keys():
                            newRandomSubnetList.append(newGene)"""
                # print(gene + " " + newGene)
            else:
                gene = item


def main():
    ##TEST TO SEE HOW MANY ROWS ARE IN NONFAGENES
    """with open("STRING 1.txt", "r") as file:
    se = set()
    results = [row.split("\t")[:2] for row in file]
    for i in results:
        for j in i:
            se.add(j)
    print(len(se))"""

    bins = Bins("STRING 1.txt")
    faGenes = FaGenes("STRING 1.txt")
    nonfaGenes = NonFaGenes("STRING 1.txt", faGenes=faGenes)

    parentSubnetwork = Stage1_SubNetworks("results.txt", "Input.gmt.txt")
    print("main")

    create_secondary_subnetwork(
        parentSubnetwork=parentSubnetwork, nonfaGenes=nonfaGenes, bins=bins
    )


if __name__ == "__main__":
    main()
