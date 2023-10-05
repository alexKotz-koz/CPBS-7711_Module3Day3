import random
import json


def fanconi_anemia_genes(inputFile):
    faData = []
    faGenes = []

    with open(inputFile, "r") as file:
        faData = [row.strip().split("\t")[2:] for row in file]

    for row in faData:
        faGenes.extend(row)

    return faGenes


def create_loci(inputFile):
    faGenes = dict()
    with open(inputFile) as file:
        listFaGenes = [row.split("\t") for row in file]

        for row in listFaGenes:
            key = row[0][-2:].strip()
            genes = row[2:]
            genes = [gene.strip() for gene in genes]
            faGenes.update({key: {"genes": genes}})

    return faGenes


def extract_fa_genes(prevFaSubnetworkFile):
    module1FASubnetwork = []
    with open(prevFaSubnetworkFile, "r") as file:
        for row in file:
            row = row.split("\t")
            row[2] = row[2].strip()
            module1FASubnetwork.append(row)
    return module1FASubnetwork


def extract_nonfa_genes(inputFile, faGenes):
    """nonfaGenes = set()
    faGenes = set(faGenes)
    with open(inputFile, "r") as file:
        parentNetwork = [row.strip().split("\t")[:1] for row in file]
        allGenes = set(
            [gene for row in parentNetwork for gene in row]
        )  # Create a set of all genes in the parent network
        print(len(allGenes))
        nonfaGenes = allGenes - faGenes"""

    nonfaGenes = set()
    faGenes = set(faGenes)
    with open(inputFile, "r") as file:
        parentNetwork = [row.split("\t") for row in file]
        print
        for row in parentNetwork:
            if row[0] not in faGenes:
                nonfaGenes.add(row[0])
            if row[1] not in faGenes:
                nonfaGenes.add(row[1])
        # HOW MANY LINES DID EVERYONE ELSE GET
        print(len(nonfaGenes))

    return nonfaGenes


def generate_12_genes():
    faGenes = create_loci("Input.gmt.txt")
    genesForSubnetwork = set()
    for index, item in enumerate(faGenes):
        random_int = random.randint(0, len(faGenes[item]["genes"]))
        try:
            genesForSubnetwork.add(faGenes[item]["genes"][random_int])
        except IndexError:
            random_int = random.randint(0, len(faGenes[item]["genes"]) / 2)
            genesForSubnetwork.add(faGenes[item]["genes"][random_int])
    return list(genesForSubnetwork)
    """[print(f"row{row}") for row in genesForSubnetwork]
    print(len(genesForSubnetwork'))"""


def create_individual_subnetwork(module1FASubnetwork):
    subnetworkToWrite = []
    flattenedSubnetwork = []
    """geneSet12 = [
        "PPM1D",
        "ERN1",
        "CD79B",
        "CD19",
        "XPO6",
        "EFCAB3",
        "MRPS31",
        "RPS15A",
        "KIF6",
        "METTL15",
        "SLX4",
        "CORO7",
    ]"""
    geneSet12 = generate_12_genes()

    for gene in geneSet12:
        for row in module1FASubnetwork:
            if gene == row[0] and row[1] in geneSet12:
                subnetworkToWrite.append(row)
            elif gene == row[1] and row[0] in geneSet12:
                subnetworkToWrite.append(row)

    for item in subnetworkToWrite:
        for gene in item:
            flattenedSubnetwork.append(gene)

    for gene in geneSet12:
        if gene not in flattenedSubnetwork:
            subnetworkToWrite.append(gene)
    return subnetworkToWrite


def create_random_subnetworks():
    module1FASubnetwork = extract_fa_genes("results.txt")

    finalList = []
    finalDictionary = {}

    print(len(module1FASubnetwork))

    i = 0

    while i < 5000:
        individualSubnetwork = []
        individualSubnetwork = create_individual_subnetwork(
            module1FASubnetwork=module1FASubnetwork
        )
        finalList.append(individualSubnetwork)
        i += 1

    for index, item in enumerate(finalList):
        index = str(index)
        finalDictionary.update({index: item})
    with open("random_subnetworks.json", "w") as outputFile:
        json.dump(finalDictionary, outputFile)

    return finalList


def create_secondary_subnetwork(parentSubnetwork, nonfaGenes):
    print(len(nonfaGenes))
    """for row in nonfaGenes:
        print(row)"""


def main():
    faGenes = fanconi_anemia_genes("STRING 1.txt")
    nonfaGenes = extract_nonfa_genes("STRING 1.txt", faGenes=faGenes)
    parentSubnetwork = create_random_subnetworks()
    create_secondary_subnetwork(
        parentSubnetwork=parentSubnetwork, nonfaGenes=nonfaGenes
    )


if __name__ == "__main__":
    main()
