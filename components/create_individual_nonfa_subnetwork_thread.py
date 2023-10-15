import threading
import random
import multiprocessing


class Create_Individual_Nonfa_Subnetwork_Thread(threading.Thread):
    def __init__(self, subnet, nonfaBin, bins, parentNetwork, stage1Subnetwork):
        threading.Thread.__init__(self)
        self.subnet = subnet
        self.nonfaBin = nonfaBin
        self.bins = bins
        self.parentNetwork = parentNetwork
        self.stage1Subnetwork = stage1Subnetwork
        self.manager = multiprocessing.Manager()
        self.queue = self.manager.Queue()
        self.result = None

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
        for row in self.parentNetwork:
            if row[0] in subnetwork:
                if row[1] in subnetwork:
                    edgeCount += 1
            elif row[1] in subnetwork:
                if row[0] in subnetwork:
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
        # print(stage1Subnetwork)

        for gene in stage1Subnetwork:
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
        print("run called")
        subnetEdgeCount = self.count_edges(subnet)

        # print(f"subnetEdgeCount: {subnetEdgeCount}")
        # print(f"subnet: {subnet}")
        # print(f"faGeneBinFlag: {faGeneBinFlag}")
        # print(f"binNotFoundFlag: {binNotFoundFlag}")

        result = {
            "edgeCount": subnetEdgeCount,
            "subnet": subnet,
            "faGeneBinFlag": faGeneBinFlag,
            "binNotFoundFlag": binNotFoundFlag,
        }

        # print(f"result from class: {result}")

        self.queue.put(result)

    def get_result(self):
        self.result = self.queue.get()
        return self.result
