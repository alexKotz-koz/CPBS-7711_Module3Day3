import threading
import random


class Create_Individual_Nonfa_Subnetwork_Thread(threading.Thread):
    def __init__(self, subnet, nonfaBin, bins, parentNetwork, stage1Subnetwork):
        threading.Thread.__init__(self)
        self.subnet = subnet
        self.nonfaBin = nonfaBin
        self.bins = bins
        self.parentNetwork = parentNetwork
        self.stage1Subnetwork = stage1Subnetwork

    def find_bin(self, gene, bins):
        binToReturn = {}
        edgeCount = ""
        for bin in bins:
            if gene in bins[bin]:
                binToReturn[bin] = bins[bin]
                edgeCount = bin
                break
        # print(f"findbin | ogGene: {gene} binToReturn:{binToReturn} edgeCount: {edgeCount}")
        return binToReturn, edgeCount

    def count_edges(self, subnetwork):
        edgeCount = 0
        for row in self.stringNetwork:
            if row[0] in self.subnetwork:
                if row[1] in self.subnetwork:
                    edgeCount += 1
            elif row[1] in self.subnetwork:
                if row[0] in self.subnetwork:
                    edgeCount += 1
        return edgeCount

    def create_individual_nonfa_subnetwork(
        self, subnet, nonfaBin, bins, stage1Subnetwork
    ):
        tempFlattendSubnetwork = set()
        subnet = set()
        newNonFaGeneBin = {}
        binNotFoundFlag = False
        faGeneBinFlag = False
        faGeneBin = bins[0]
        ##########################################
        print(stage1Subnetwork)

        for gene in stage1Subnetwork["subnet"]:
            if isinstance(gene, list):
                for subGene in gene:
                    tempFlattendSubnetwork.add(subGene)
            else:
                tempFlattendSubnetwork.add(gene)
        print("Finding bins...")

        for gene in tempFlattendSubnetwork:
            geneBin, binEdgeCount = self.find_bin(gene, bins)

            # if the gene from tempFlattendSubnetwork lives in a bin with no nonfa genes, a random nonfa gene is used in place, not from the same bin. add flag to this subnetwork
            if binEdgeCount == 0:
                faGeneBinFlag = True

            if binEdgeCount == "":
                print("Bin Edge Count null")
                binNotFoundFlag = True
                continue

            if binEdgeCount in nonfaBin:
                nonfaBinGenes = nonfaBin[binEdgeCount]
                newGene = random.choice(nonfaBinGenes)
                newNonFaGeneBin[gene] = newGene
                subnet.add(newGene)

        return list(subnet), binNotFoundFlag, faGeneBinFlag

    def run(self):
        (
            subnet,
            binNotFoundFlag,
            faGeneBinFlag,
        ) = self.create_individual_nonfa_subnetwork(
            self.subnet, self.nonfaBin, self.bins, self.stage1Subnetwork
        )
        subnetEdgeCount = self.count_edges(subnet, self.parentNetwork)

        result = {
            "edgeCount": subnetEdgeCount,
            "subnet": subnet,
            "faGeneBinFlag": faGeneBinFlag,
            "binNotFoundFlag": binNotFoundFlag,
        }

        return result
