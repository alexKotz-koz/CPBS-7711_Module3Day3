import random
from bisect import bisect_left


class Create_Individual_Nonfa_Subnetwork_Thread:
    def __init__(self, nonfaBin, bins, parentNetwork, stage1Subnetwork):
        self.nonfaBin = nonfaBin
        self.bins = bins
        self.parentNetwork = parentNetwork
        self.stage1Subnetwork = stage1Subnetwork

        self.result = {}

    # Input: stage 1 subnetwork genes and the bins object
    def find_bins(self, genes, bins):
        binsToReturn = {}
        count = 1
        for gene in genes:
            for binName, binGenes in bins.items():
                sortedGenes = sorted(binGenes)
                index = bisect_left(sortedGenes, gene)
                if index != len(sortedGenes) and sortedGenes[index] == gene:
                    count += 1
                    # binsToReturn.setdefault(binName, []).append(gene)
                    binsToReturn[gene] = {binName: binGenes}
                    break
        print(count)
        print(len(binsToReturn))
        return binsToReturn

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

    def create_individual_nonfa_subnetwork(self, nonfaBin, bins, stage1Subnetwork):
        tempFlattendSubnetwork = set()
        subnet = []
        binNotFoundFlag = False
        faGeneBinFlag = False
        faGeneBin = bins[0]

        tempFlattendSubnetwork = set(
            gene
            for sublist in stage1Subnetwork
            for gene in (sublist if isinstance(sublist, list) else [sublist])
        )

        # Find the bins for all genes in the tempFlattendSubnetwork list
        binsToGenes = self.find_bins(tempFlattendSubnetwork, bins)

        for key, val in binsToGenes.items():
            newGeneList = []
            for name, genes in val.items():
                binName = name
                binGenes = genes
            # Check if the bin has any nonfa genes
            if binName in nonfaBin:
                newGeneList = nonfaBin[binName]
            else:
                faGeneBinFlag = True
                newGeneList.append("faGene")

            # Choose a random gene from the gene list
            if newGeneList == "faGene":
                subnet.add("faGene")
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
        # print("run called")
        subnetEdgeCount = self.count_edges(subnet)

        result = {
            "edgeCount": subnetEdgeCount,
            "subnet": subnet,
            "faGeneBinFlag": faGeneBinFlag,
            "binNotFoundFlag": binNotFoundFlag,
        }

        # print(f"result from class: {result}")
        return result
