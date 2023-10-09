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
    nonfaGenesSet = set()
    nonfaGenes = {}
    faGenes = set(faGenes)

    with open(inputFile, "r") as file:
        parentNetwork = [row.split("\t") for row in file]

    for row in parentNetwork:
        if row[0] not in faGenes:
            if row[0] not in nonfaGenes:
                nonfaGenes[row[0]] = 0
            nonfaGenes[row[0]] += 1
        if row[1] not in faGenes:
            if row[1] not in nonfaGenes:
                nonfaGenes[row[1]] = 0
            nonfaGenes[row[1]] += 1

    return nonfaGenes


def create_bins(inputFile):
    countPerGene = {}
    bins = {}
    with open(inputFile, "r") as file:
        results = [row.split("\t")[:2] for row in file]

    for row in results:
        if row[0] not in countPerGene:
            countPerGene[row[0]] = 0
        elif row[0] in countPerGene:
            countPerGene[row[0]] += 1

        if row[1] not in countPerGene:
            countPerGene[row[1]] = 0
        elif row[1] in countPerGene:
            countPerGene[row[1]] += 1

    sorted_dict = dict(sorted(countPerGene.items(), key=lambda item: item[1]))

    for item in sorted_dict:
        bins.setdefault(sorted_dict[item], []).append(item)

    with open("bins.json", "w") as outputFile:
        json.dump(bins, outputFile)

    return bins


def generate_12_genes():
    faGenes = create_loci("Input.gmt.txt")

    genesForSubnetwork = set()

    ##REFACTOR
    for index, item in enumerate(faGenes):
        random_int = random.randint(0, len(faGenes[item]["genes"]))
        try:
            genesForSubnetwork.add(faGenes[item]["genes"][random_int])
        except IndexError:
            random_int = random.randrange(0, len(faGenes[item]["genes"]))
            genesForSubnetwork.add(faGenes[item]["genes"][random_int])
    return list(genesForSubnetwork)


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

    subnetworkToWrite = [
        sublist
        for i, sublist in enumerate(subnetworkToWrite)
        if sublist not in subnetworkToWrite[:i]
    ]

    return subnetworkToWrite


def create_random_subnetworks():
    module1FASubnetwork = extract_fa_genes("results.txt")

    finalList = []
    finalDictionary = {}
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

    return finalDictionary


# Mating or 1:1 (fa to nonfa) replacement
def create_secondary_subnetwork(parentSubnetwork, nonfaGenes, bins):
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
            print(len(bins))
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

    faGenes = fanconi_anemia_genes("STRING 1.txt")
    nonfaGenes = extract_nonfa_genes("STRING 1.txt", faGenes=faGenes)
    parentSubnetwork = create_random_subnetworks()
    bins = create_bins("STRING 1.txt")
    create_secondary_subnetwork(
        parentSubnetwork=parentSubnetwork, nonfaGenes=nonfaGenes, bins=bins
    )


if __name__ == "__main__":
    main()
