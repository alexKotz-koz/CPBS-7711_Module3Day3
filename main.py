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


def count_edges(subnetwork, parentNetwork):
    edgeCount = 0
    for row in parentNetwork:
        if row[0] in subnetwork:
            if row[1] in subnetwork:
                edgeCount += 1
        elif row[1] in subnetwork:
            if row[0] in subnetwork:
                edgeCount += 1
    return edgeCount


"""def find_new_nonFA_gene(countFaNodeFlag, bins):
    for bin in bins:
        if any(
            "geneType" in geneDict and geneDict["geneType"] == "nonfaGene"
            for geneDict in bins[bin]
        ):
            print(bins[bin])"""


def create_individual_nonfa_subnetwork(stage1Subnetwork, parentNetwork, bins):
    tempFlattendSubnetwork = set()
    newGeneSet = []
    newGeneSetToWrite = []

    for gene in stage1Subnetwork["subnet"]:
        if isinstance(gene, list):
            for subGene in gene:
                tempFlattendSubnetwork.add(subGene)
        else:
            tempFlattendSubnetwork.add(gene)

    # print("Finding bins...")

    for gene in tempFlattendSubnetwork:
        tempNonFaGenes = []
        for bin in bins:
            # HOW LONG is ^ Running for
            binNumber = bin
            binObject = bins[bin]
            binGenes = [list(key.keys())[0] for key in bins[bin]]
            if gene in binGenes:
                for index, geneDict in enumerate(binObject):
                    geneKey = list(geneDict.keys())[0]
                    if geneDict[geneKey]["geneType"] == "nonfaGene":
                        tempNonFaGenes.append(geneKey)
        if len(tempNonFaGenes) == 0:
            newGeneSet.append("faNodeFlag")
            # print("single fa node in stage 1 subnetwork")
            continue
        newGeneNum = random.randrange(0, len(tempNonFaGenes))
        newGeneSet.append(tempNonFaGenes[newGeneNum])

    """newGeneSet = [
        "PLIN3",
        "PEX11G",
        "faNodeFlag",
        "E2F5",
        "RNF128",
        "STON2",
        "PRKD1",
        "CDK17",
        "DAZ2",
        "TMSB15B",
        "PARL",
        "MMP19",
    ]"""

    # add connections if exist
    for row in parentNetwork:
        if row[0] in newGeneSet:
            if row[1] in newGeneSet:
                newGeneSetToWrite.append(row)
        elif row[1] in newGeneSet:
            if row[0] in newGeneSet:
                newGeneSetToWrite.append(row)
    # print(newGeneSetToWrite)
    # sort and unduplicate newGeneSetToWrite
    newGeneSetToWrite = sorted(newGeneSetToWrite, key=lambda x: x[0], reverse=True)
    newGeneSetToWrite = list(
        set(tuple(sorted(subList)) for subList in newGeneSetToWrite)
    )
    newGeneSetToWrite = [list(subList) for subList in newGeneSetToWrite]

    # flatten newGeneSetToWrite and add the connected list(if exists) and the reset of the genes fron newGeneSet
    flattendNewGeneSetToWrite = []
    for row in newGeneSetToWrite:
        if row[0] not in flattendNewGeneSetToWrite:
            flattendNewGeneSetToWrite.append(row[0])
        if row[1] not in flattendNewGeneSetToWrite:
            flattendNewGeneSetToWrite.append(row[1])
    for row in newGeneSet:
        if row not in flattendNewGeneSetToWrite:
            newGeneSetToWrite.append(row)
    # print(newGeneSetToWrite)
    print(newGeneSetToWrite)

    return newGeneSetToWrite

    ### WORK-AROUND FOR CASE WHEN GENE IN SUBNET IS IN FA ONLY BIN
    # find faNodeFlag, replace with (a nonfa gene from a random bin)

    """countFaNodeFlag = 0
    for gene in newGeneSetToWrite:
        if gene == "faNodeFlag":
            countFaNodeFlag += 1
    find_new_nonFA_gene(countFaNodeFlag, bins)"""


def create_secondary_subnetwork(parentNetworkFile, stage1Subnetworks, bins, faGenes):
    # print("Creating stage 2 random subnetworks")

    stage2Subnetwork = {}
    # nonfaGenes = set(nonfaGenes.keys())
    parentNetwork = []

    with open("STRING 1.txt", "r") as file:
        parentNetwork = [row.split("\t")[:2] for row in file]

    for index, subnet in stage1Subnetworks.items():
        subnetworksFromStage1 = subnet["subnet"]

        # flatten stage 1 subnetwork (some subnetworks contain sublists, representing an existing connection between two genes)

        subnet = create_individual_nonfa_subnetwork(subnet, parentNetwork, bins)
        subnetEdgeCount = count_edges(subnet, parentNetwork)

        stage2Subnetwork[index] = {"edgeCount": subnetEdgeCount, "subnet": subnet}
        # subnet.append(newGene)

    with open("stage2_random_subnetworks.json", "w") as outputFile:
        json.dump(stage2Subnetwork, outputFile)
    # print("Second 5,000 subnetworks created")


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
    bins = binsInstance.create_bins()

    stage1_subnetworksInstance = Stage1_SubNetworks(
        "results.txt", "Input.gmt.txt", "STRING 1.txt"
    )

    (
        stage1Subnetworks,
        edgeCount,
    ) = stage1_subnetworksInstance.create_random_subnetworks()
    stage2_subnetworks = create_secondary_subnetwork(
        "STRING 1.txt", stage1Subnetworks, bins, faGenes
    )


if __name__ == "__main__":
    main()
