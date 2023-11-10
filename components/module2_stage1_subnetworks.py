import json
import random
from components.fa_utilities import FaUtilities


class Stage1_SubNetworks:
    def __init__(
        self, prevFaSubnetworkFile, faInputFile, stringInputFile, parentNetworkDict
    ):
        self.prevFaSubnetworkFile = prevFaSubnetworkFile
        self.faInputFile = faInputFile
        self.stringInputFile = stringInputFile
        self.parentNetworkDict = parentNetworkDict
        self.module1FASubnetwork = []
        self.faGenes = []
        self.sortedDictionary = {}

    # Input: Module 1 Day 3 fa subnetwork file
    # Output: An object containing sub-dictionaries for each fa locus
    def create_loci(self, faInputFile):
        faGenes = {}
        with open(self.faInputFile) as file:
            listFaGenes = [row.split("\t") for row in file]

            for row in listFaGenes:
                key = row[0][-2:].strip()
                genes = row[2:]
                genes = [gene.strip() for gene in genes]
                faGenes.update({key: {"genes": genes}})

        return faGenes

        def count_edges(self):
            sorted_dict = self.sortedDictionary
            edgeCountDict = {}
            totalEdgeCount = 0
            one = 0
            two = 0
            three = 0
            four = 0
            five = 0
            zero = 0
            six = 0
            for index, item in sorted_dict.items():
                if item["edgeCount"] == 0:
                    zero += 1
                elif item["edgeCount"] == 1:
                    one += 1
                elif item["edgeCount"] == 2:
                    two += 1
                elif item["edgeCount"] == 3:
                    three += 1
                elif item["edgeCount"] == 4:
                    four += 1
                elif item["edgeCount"] == 5:
                    five += 1
                elif item["edgeCount"] == 6:
                    six += 1
            edgeCountDict["zero"] = zero
            edgeCountDict["one"] = one
            edgeCountDict["two"] = two
            edgeCountDict["three"] = three
            edgeCountDict["four"] = four
            edgeCountDict["five"] = five
            edgeCountDict["six"] = six

            for index, item in sorted_dict.items():
                totalEdgeCount += item["edgeCount"]

            edgeCountDict["totalEdgeCount"] = totalEdgeCount

            return edgeCountDict

    # Input: loci dictionary (faGenes)
    def generate_12_genes(self):
        faGenes = self.create_loci(self)

        genesForSubnetwork = set()
        # for each locus in the faGenes dictionary, extract one gene at random, from each locus and store in a new list
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

                module1FASubnetwork.append(row)
        return module1FASubnetwork

    def create_individual_subnetwork(self, module1FASubnetwork):
        subnetworkToWrite = []
        flattenedSubnetwork = []
        geneSet12 = self.generate_12_genes()

        for gene in geneSet12:
            for row in module1FASubnetwork:
                if gene == row[0] and row[1] in geneSet12:
                    subnetworkToWrite.append(row[:2])
                elif gene == row[1] and row[0] in geneSet12:
                    subnetworkToWrite.append(row[:2])

        for item in subnetworkToWrite:
            for gene in item:
                flattenedSubnetwork.append(gene)

        for gene in geneSet12:
            if gene not in flattenedSubnetwork:
                subnetworkToWrite.append(gene)

        subnetworkToWrite = [
            sublist
            for index, sublist in enumerate(subnetworkToWrite)
            if sublist not in subnetworkToWrite[:index]
        ]

        return subnetworkToWrite

    def create_random_subnetworks(self):
        print("Creating stage 1 random subnetworks...")
        module1FASubnetwork = self.extract_fa_genes()

        finalList = []
        finalDictionary = {}

        count = 0
        while count < 5000:
            individualSubnetwork = []
            individualSubnetwork = self.create_individual_subnetwork(
                module1FASubnetwork
            )
            finalList.append(individualSubnetwork)
            count += 1

        # edge count
        # print(f"m2 parentnetworkdict: {type(self.parentNetworkDict)}")
        # print(f"m2 individualSubnet = {finalList}")
        ###REMOVE
        for index, item in enumerate(finalList):
            """faUtilitiesInstance = FaUtilities(
                parentNetworkFile=self.parentNetworkDF, individualSubnetwork=item
            )
            edgeCount = faUtilitiesInstance.count_edges()"""
            index = str(index)

            finalDictionary[index] = {"subnet": item}

        with open("stage1_random_subnetworks.json", "w") as outputFile:
            json.dump(finalDictionary, outputFile)
        print("First 5,000 subnetworks created")

        return finalDictionary
