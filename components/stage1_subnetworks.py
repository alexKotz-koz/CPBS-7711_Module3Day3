import json
import random


class Stage1_SubNetworks:
    def __init__(self, prevFaSubnetworkFile, faInputFile, stringInputFile):
        self.prevFaSubnetworkFile = prevFaSubnetworkFile
        self.faInputFile = faInputFile
        self.stringInputFile = stringInputFile
        self.module1FASubnetwork = []
        self.faGenes = []
        self.sortedDictionary = {}

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

        for item in sorted_dict.items():
            totalEdgeCount += item["edgeCount"]

        edgeCountDict["totalEdgeCount"] = totalEdgeCount

        return edgeCountDict

    def create_individual_subnetwork(self, module1FASubnetwork):
        subnetworkToWrite = []
        flattenedSubnetwork = []
        parentNetwork = []
        geneSet12 = self.generate_12_genes()

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

        for index, item in enumerate(finalList):
            index = str(index)
            edgeCount = 0
            for gene in item:
                if isinstance(gene, list):
                    edgeCount += 1
            finalDictionary[index] = {"edgeCount": edgeCount, "subnet": item}
        sortedDictionary = dict(
            sorted(finalDictionary.items(), key=lambda x: x[1]["edgeCount"])
        )

        sortedDictionary = dict(
            sorted(finalDictionary.items(), key=lambda x: x[1]["edgeCount"])
        )

        with open("stage1_random_subnetworks.json", "w") as outputFile:
            json.dump(sortedDictionary, outputFile)
        print("First 5,000 subnetworks created")

        self.sortedDictionary = sortedDictionary
        edgeCount = self.count_edges()

        return finalDictionary, edgeCount
