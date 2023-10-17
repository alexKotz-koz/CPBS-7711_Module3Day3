import random
from bisect import bisect_left


class Create_Individual_Nonfa_Subnetwork_Thread:
    def __init__(self, nonfaBin, bins, parentNetwork, stage1Subnetwork):
        self.nonfaBin = nonfaBin
        self.bins = bins
        self.parentNetwork = parentNetwork
        self.stage1Subnetwork = stage1Subnetwork

        self.result = {}

    # Input: genes, stage 1 subnetwork genes; the bins object
    # Output: dictionary that contains twelve bins the correspond to the 12 genes from the Input
    def find_bins(self, genes, bins):
        binsToReturn = {}

        # for each of the 12 genes, find the bin which each gene is located in.
        for gene in genes:
            for binNumber, binGenes in bins.items():
                # sort binGenes for quicker access
                sortedGenes = sorted(binGenes)

                # using bisect_left() function from the bisect module to find the index where gene should be inserted into sortedGenes to maintain its sorted order.
                index = bisect_left(sortedGenes, gene)

                # if the legth of sortedGenes is not equal to index (the gene is not in sortedGenes) and if the gene is in sortedGenes at the specified index
                if index != len(sortedGenes) and sortedGenes[index] == gene:
                    binsToReturn[gene] = {binNumber: binGenes}
                    break
        return binsToReturn

    # Input: subnetwork, individual nonfa subnetwork
    # Output: total number of found connections between any of the 12 genes in the subnetwork
    def count_edges(self, subnetwork):
        edgeCount = 0
        for row in self.parentNetwork:
            if row[0] in subnetwork:
                if row[1] in subnetwork:
                    edgeCount += 1
            elif row[1] in subnetwork:
                if row[0] in subnetwork:
                    edgeCount += 1
        return edgeCount

    # Input: nonfaBin, used for gene replacement; bins, passed to find_bins(); stage1Subnetwork, a single stage 1 subnetwork
    # Output: individual randomized nonfa subnetwork
    def create_individual_nonfa_subnetwork(self, nonfaBin, bins, stage1Subnetwork):
        tempFlattendSubnetwork = set()
        subnet = []
        binNotFoundFlag = False
        faGeneBinFlag = False
        faGeneBin = bins[0]

        # flatten stage 1 subnetwork to get a set of genes, some of the stage 1 subnetworks have embedded lists which represent a connection between the two genes in the sublist
        tempFlattendSubnetwork = set(
            gene
            for sublist in stage1Subnetwork
            for gene in (sublist if isinstance(sublist, list) else [sublist])
        )

        # find the bins for all genes in the tempFlattendSubnetwork list
        binsToGenes = self.find_bins(tempFlattendSubnetwork, bins)

        # iterate over each genes bin, find a nonfa gene within each bin and add the 12 nonfa genes to a list
        for key, val in binsToGenes.items():
            newGeneList = []
            for name, genes in val.items():
                binNumber = name
                binGenes = genes

            # if the bin number (edge count) exists in the nonfaBin object, add the genes from that nonfaBin bin to a list
            if binNumber in nonfaBin:
                newGeneList = nonfaBin[binNumber]
            else:
                # if the bin number (edge count) does not exist in the nonfaBin object, add "faGene" to the new gene list
                # adding the "faGene" string to the list to 1) indicate no nonfa gene was found in the stage 1 genes bin, and
                # 2) since this fa gene will not be used in the final edge count for each subnetwork, adding the string will have no affect on the outcome.
                faGeneBinFlag = True
                newGeneList.append("faGene")

            # if the newGeneList contains nonfa genes, randomly select one of the nonfa genes from the list and add to subnet list
            if newGeneList == "faGene":
                subnet.append("faGene")
            else:
                subnet.append(random.choice(newGeneList))
        print(f"subnet created: {len(subnet)}")

        return list(subnet), binNotFoundFlag, faGeneBinFlag

    def run(self):
        (
            subnet,
            binNotFoundFlag,
            faGeneBinFlag,
        ) = self.create_individual_nonfa_subnetwork(
            self.nonfaBin, self.bins, self.stage1Subnetwork
        )
        subnetEdgeCount = self.count_edges(subnet)

        result = {
            "edgeCount": subnetEdgeCount,
            "subnet": subnet,
            "faGeneBinFlag": faGeneBinFlag,
            "binNotFoundFlag": binNotFoundFlag,
        }

        # print(f"result from class: {result}")
        return result
