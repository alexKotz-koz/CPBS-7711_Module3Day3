import random


def create_loci(inputFile):
    faGenes = dict()
    with open(inputFile) as file:
        listFaGenes = [row.split("\t") for row in file]

        for row in listFaGenes:
            key = row[0][-2:].strip()
            genes = row[2:]
            genes = [gene.strip() for gene in genes]
            faGenes.update({key: {"genes": genes}})

        """for key, value in faGenes.items():
            print(key, value)
            print("\n")"""
    return faGenes


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


def create_individual_subnetwork(prevFaSubnetwork):
    subnetworkToWrite = []
    module1FASubnetwork = []
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

    prevFaSubnetworkFile = open(prevFaSubnetwork, "r")
    # finalSubnetworks = open("output.txt", "w")

    geneSet12 = generate_12_genes()

    for row in prevFaSubnetworkFile:
        row = row.split("\t")
        row[2] = row[2].strip()
        module1FASubnetwork.append(row)
    prevFaSubnetworkFile.close()

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
    """for row in module1FASubnetwork:
        # A: fully inclusive
       if (row[0] in geneSet12) and (row[1] in geneSet12):
        print(f"Row from mod 1 in subnetwork {row}")
        subnetworkToWrite.append(row)
        # B
        if row[0] in subnetwork or row[1] in subnetwork:
            print(f"row 1: {row}")
            flattenedSubnetwork.append(row)
            # finalSubnetworks.write(row + "\n")
        print(len(flattenedSubnetwork))"""

    # A
    """for gene in geneSet12:
        if len(subnetworkToWrite) == 0:
            subnetworkToWrite.append(gene + "\n")
        else:
            if isinstance(gene, list):
                if (
                    gene[0] not in subnetworkToWrite
                    and gene[1] not in subnetworkToWrite
                ):
                    subnetworkToWrite.append(gene + "\n")

            elif isinstance(gene, str):
                if gene not in subnetworkToWrite:
                    subnetworkToWrite.append(gene + "\n")"""

    print(f"len of subnetwork: {len(subnetworkToWrite)} | {subnetworkToWrite}")
    return subnetworkToWrite


def create_random_subnetworks(prevFaSubnetwork, outputFile):
    print("here")


def main():
    create_individual_subnetwork("results.txt")


if __name__ == "__main__":
    main()
