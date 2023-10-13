import json


class Bins:
    def __init__(self, stringInputFile, faInputFile, faGenes, nonfaGenes):
        self.stringInputFile = stringInputFile
        self.faInputFile = faInputFile
        self.faGenes = faGenes
        self.nonfaGenes = nonfaGenes

    def create_bins(self):
        print("Creating bins...")
        countPerGene = {}
        nonfaBin = {}
        faBin = {}
        bins = {}
        nfabins = {}
        seen = {}
        uniqueResults = []
        with open(self.stringInputFile, "r") as file:
            results = [row.split("\t")[:2] for row in file]
            results = set(tuple(row) for row in results)

            # Only add unique rows (a,b) (b,a) -- only add b,a to results
            for row in results:
                gene1, gene2 = row
                if (gene2, gene1) not in seen:
                    seen[row] = True
                    seen[(gene2, gene1)] = True
                    uniqueResults.append(tuple(row))

            """with open("temp.txt", "w") as file:
                for tup in uniqueResults:
                    file.write(str(tup) + "\n")"""

        ##REFACTOR lines 24-48

        # bin all genes from master network (STRING 1.txt)
        for row in uniqueResults:
            if row[0] not in countPerGene:
                countPerGene[row[0]] = 0
            elif row[0] in countPerGene:
                countPerGene[row[0]] += 1

            if row[1] not in countPerGene:
                countPerGene[row[1]] = 0
            elif row[1] in countPerGene:
                countPerGene[row[1]] += 1

        seen2 = {}
        uniqueResults2 = []
        # bin all genes from module 1 subnetwork (to include all faGenes in bins object)
        with open(self.faInputFile, "r") as file:
            results2 = [row.split("\t")[:2] for row in file]
            results2 = set(tuple(row) for row in results2)
            for row in results2:
                gene1, gene2 = row
                if (gene2, gene1) not in seen:
                    seen[row] = True
                    seen[(gene2, gene1)] = True
                    uniqueResults2.append(tuple(row))

        for row in uniqueResults2:
            if row[0] not in countPerGene:
                countPerGene[row[0]] = 0
            elif row[0] in countPerGene:
                countPerGene[row[0]] += 1

            if row[1] not in countPerGene:
                countPerGene[row[1]] = 0
            elif row[1] in countPerGene:
                countPerGene[row[1]] += 1

        # print(len(self.nonfaGenes))
        # print(len(self.faGenes))
        # print(self.nonfaGenes)
        for item in countPerGene:
            if item in self.nonfaGenes:
                if item == "OMP":
                    print(f"item{countPerGene[item]}")
                nonfaBin[item] = countPerGene[item]

        sortedNonFaDict = dict(sorted(nonfaBin.items(), key=lambda item: item[1]))

        sortedDict = dict(sorted(countPerGene.items(), key=lambda item: item[1]))

        for item in sortedNonFaDict:
            nfabins.setdefault(sortedNonFaDict[item], []).append(item)

        for item in sortedDict:
            bins.setdefault(sortedDict[item], []).append(item)

        """print(len(bins))
        print(len(nfabins))"""

        with open("nfabins.json", "w") as outputFile:
            json.dump(nfabins, outputFile)

        with open("bins.json", "w") as outputFile:
            json.dump(bins, outputFile)
        print("Bins created")
        return bins, nfabins
