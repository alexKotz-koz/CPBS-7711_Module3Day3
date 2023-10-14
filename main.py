import random
import json
from components.bins import Bins
from components.fagenes import FaGenes
from components.nonfagenes import NonFaGenes
from components.stage1_subnetworks import Stage1_SubNetworks


def find_bin(gene, bins):
    binToReturn = {}
    edgeCount = ""
    for bin in bins:
        if gene in bins[bin]:
            binToReturn[bin] = bins[bin]
            edgeCount = bin
            break
    # print(f"findbin | ogGene: {gene} binToReturn:{binToReturn} edgeCount: {edgeCount}")
    return binToReturn, edgeCount


def count_edges(subnetwork, stringNetwork):
    edgeCount = 0
    for row in stringNetwork:
        if row[0] in subnetwork:
            if row[1] in subnetwork:
                edgeCount += 1
        elif row[1] in subnetwork:
            if row[0] in subnetwork:
                edgeCount += 1
    return edgeCount


def create_individual_nonfa_subnetwork(stage1Subnetwork, nonfaBin, bins):
    tempFlattendSubnetwork = set()
    subnet = set()
    newNonFaGeneBin = {}
    binNotFoundFlag = False
    faGeneBinFlag = False
    faGeneBin = bins[0]

    for gene in stage1Subnetwork["subnet"]:
        if isinstance(gene, list):
            for subGene in gene:
                tempFlattendSubnetwork.add(subGene)
        else:
            tempFlattendSubnetwork.add(gene)
    print("Finding bins...")

    for gene in tempFlattendSubnetwork:
        geneBin, binEdgeCount = find_bin(gene, bins)

        # if the gene from tempFlattendSubnetwork lives in a bin with no nonfa genes, a random nonfa gene is used in place, not from the same bin. add flag to this subnetwork
        if binEdgeCount == 0:
            faGeneBinFlag = True

        if binEdgeCount == "":
            print("Bin Edge Count null")
            binNotFoundFlag = True
            continue

        if binEdgeCount not in nonfaBin:
            print("here")

        if binEdgeCount in nonfaBin:
            nonfaBinGenes = nonfaBin[binEdgeCount]
            newGene = random.choice(nonfaBinGenes)
            newNonFaGeneBin[gene] = newGene
            subnet.add(newGene)

    return list(subnet), binNotFoundFlag, faGeneBinFlag


"""
        # Double check before creating nonfa subnetwork
    for i in subnet:
        if i in faGenes:
            print(i)"""


def create_secondary_subnetwork(
    parentNetworkFile, stage1Subnetworks, nonfaBin, bins, faGenes
):
    print("Creating stage 2 random subnetworks")

    stage2Subnetwork = {}
    # nonfaGenes = set(nonfaGenes.keys())
    parentNetwork = []

    with open("STRING 1.txt", "r") as file:
        parentNetwork = [row.split("\t")[:2] for row in file]

    for index, subnet in stage1Subnetworks.items():
        subnetworksFromStage1 = subnet["subnet"]

        # flatten stage 1 subnetwork (some subnetworks contain sublists, representing an existing connection between two genes)

        subnet, binNotFoundFlag, faGeneBinFlag = create_individual_nonfa_subnetwork(
            subnet, nonfaBin, bins
        )
        subnetEdgeCount = count_edges(subnet, parentNetwork)

        stage2Subnetwork[index] = {
            "edgeCount": subnetEdgeCount,
            "subnet": subnet,
            "faGeneBinFlag": faGeneBinFlag,
            "binNotFoundFlag": binNotFoundFlag,
        }
        # subnet.append(newGene)

    with open("stage2_random_subnetworks.json", "w") as outputFile:
        json.dump(stage2Subnetwork, outputFile)
    print("Second 5,000 subnetworks created")


def t_test(stage1Subnetworks, stage2Subnetworks):
    pValue = 0

    return pValue


def main():
    ##TEST TO SEE HOW MANY ROWS ARE IN NONFAGENES
    """with open("STRING 1.txt", "r") as file:
    se = set()
    results = [row.split("\t")[:2] for row in file]
    for i in results:
        for j in i:
            se.add(j)
    print(len(se))"""

    faGenesInstance = FaGenes("Input.gmt.txt")
    faGenes = faGenesInstance.fanconi_anemia_genes()

    nonfaGenesInstance = NonFaGenes("STRING 1.txt", faGenes=faGenes)
    nonfaGenes = nonfaGenesInstance.extract_nonfa_genes()

    binsInstance = Bins("STRING 1.txt", "results.txt", faGenes, nonfaGenes)
    bins, nonfaBins = binsInstance.create_bins()

    stage1_subnetworksInstance = Stage1_SubNetworks(
        "results.txt", "Input.gmt.txt", "STRING 1.txt"
    )

    (
        stage1Subnetworks,
        edgeCount,
    ) = stage1_subnetworksInstance.create_random_subnetworks()
    stage2_subnetworks = create_secondary_subnetwork(
        "STRING 1.txt", stage1Subnetworks, nonfaBins, bins, faGenes
    )
    # print(edgeCount)


if __name__ == "__main__":
    main()
