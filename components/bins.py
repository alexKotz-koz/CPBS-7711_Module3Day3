import json
import pandas as pd
import numpy as np


class Bins:
    def __init__(self, parentNetwork, m1d3InputFile, faGenes, nonfaGenes):
        self.parentNetwork = pare
        self.m1d3InputFile = m1d3InputFile
        self.faGenes = faGenes
        self.nonfaGenes = nonfaGenes

    def create_bins(self):
        print("Creating bins...")
        # countPerGene = {GeneName:EdgeCount}
        countPerGene = {}
        nonfaBin = {}
        faBin = {}
        bins = {}
        nfabins = {}
        seen = {}
        uniqueResults = []

        # read in parent network and store in results set
        """with open(self.stringInputFile, "r") as file:
            results = [row.split("\t")[:2] for row in file]
            results = set(tuple(row) for row in results)

            # only add unique rows (a,b) (b,a) -- strictly add b,a to results
            for row in results:
                gene1, gene2 = row
                if (gene2, gene1) not in seen:
                    seen[row] = True
                    seen[(gene2, gene1)] = True
                    uniqueResults.append(tuple(row))"""

        # bin all genes from master network (STRING 1.txt)
        # ASSUMPTION: all genes from STRING 1.txt have at least one (node-node) connection
        for row in self.parentNetwork:
            if row[0] not in countPerGene:
                countPerGene[row[0]] = 1
            elif row[0] in countPerGene:
                countPerGene[row[0]] += 1

            if row[1] not in countPerGene:
                countPerGene[row[1]] = 1
            elif row[1] in countPerGene:
                countPerGene[row[1]] += 1

        seen2 = {}
        uniqueResults2 = []
        # bin all genes from module 1 subnetwork (to include all faGenes in bins object)
        with open(self.m1d3InputFile, "r") as file:
            results2 = [row.split("\t")[:2] for row in file]
            results2 = set(tuple(row) for row in results2)
            for row in results2:
                gene1, gene2 = row
                if (gene2, gene1) not in seen:
                    seen[row] = True
                    seen[(gene2, gene1)] = True
                    uniqueResults2.append(tuple(row))

        # ASSUMPTION: all genes from module 1 subnetwork have at least one (node-node) connection
        for row in uniqueResults2:
            if row[0] not in countPerGene:
                countPerGene[row[0]] = 1
            elif row[0] in countPerGene:
                countPerGene[row[0]] += 1

            if row[1] not in countPerGene:
                countPerGene[row[1]] = 1
            elif row[1] in countPerGene:
                countPerGene[row[1]] += 1

        for gene in self.faGenes:
            if gene not in countPerGene.keys():
                countPerGene[gene] = 0

        # create nonfaBin, based on the larger bin object
        for item in countPerGene:
            if item in self.nonfaGenes:
                nonfaBin[item] = countPerGene[item]

        sortedNonFaDict = dict(sorted(nonfaBin.items(), key=lambda item: item[1]))

        sortedDict = dict(sorted(countPerGene.items(), key=lambda item: item[1]))

        # add items sortedNonFaDict and sortedDict to corresponding objects
        for item in sortedNonFaDict:
            nfabins.setdefault(sortedNonFaDict[item], []).append(item)

        for item in sortedDict:
            bins.setdefault(sortedDict[item], []).append(item)

        # write the objects to files for use in creation of secondary subnetworks
        with open("nfabins.json", "w") as outputFile:
            json.dump(nfabins, outputFile)

        with open("bins.json", "w") as outputFile:
            json.dump(bins, outputFile)
        print("Bins created")
        return bins, nfabins
