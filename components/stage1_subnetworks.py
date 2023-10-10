import json
import random


class Stage1_SubNetworks:
    def __init__(self, prevFaSubnetworkFile, inputFile):
        self.prevFaSubnetworkFile = prevFaSubnetworkFile
        self.faInputFile = inputFile
        self.module1FASubnetwork = []
        self.faGenes = []

    def create_loci(self, faInputFile):
        faGenes = dict()
        with open(self.faInputFile) as file:
            listFaGenes = [row.split("\t") for row in file]

            for row in listFaGenes:
                key = row[0][-2:].strip()
                genes = row[2:]
                genes = [gene.strip() for gene in genes]
                faGenes.update({key: {"genes": genes}})

        return faGenes

    def generate_12_genes(self):
        faGenes = self.create_loci(self)

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

    def extract_fa_genes(self):
        module1FASubnetwork = []
        with open(self.prevFaSubnetworkFile, "r") as file:
            for row in file:
                row = row.split("\t")
                row[2] = row[2].strip()
                module1FASubnetwork.append(row)
        return module1FASubnetwork

    def create_individual_subnetwork(self, module1FASubnetwork):
        subnetworkToWrite = []
        flattenedSubnetwork = []

        geneSet12 = self.generate_12_genes()

        for gene in geneSet12:
            for row in self.module1FASubnetwork:
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

    def create_random_subnetworks(self):
        module1FASubnetwork = self.extract_fa_genes()

        finalList = []
        finalDictionary = {}
        i = 0

        while i < 5000:
            individualSubnetwork = []
            individualSubnetwork = self.create_individual_subnetwork(
                module1FASubnetwork
            )
            finalList.append(individualSubnetwork)
            i += 1

        for index, item in enumerate(finalList):
            index = str(index)
            finalDictionary.update({index: item})
        with open("random_subnetworks.json", "w") as outputFile:
            json.dump(finalDictionary, outputFile)

        return finalDictionary
