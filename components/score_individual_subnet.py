from components.fa_utilities import FaUtilities
import time
import cProfile


class ScoreIndividualSubnet:
    def __init__(self, individualSubnet, inputFile, parentNetwork):
        self.individualSubnet = individualSubnet
        self.inputFile = inputFile
        self.parentNetwork = parentNetwork
        faUtilitiesInstance = FaUtilities(inputFile=inputFile)
        self.loci = faUtilitiesInstance.extract_loci()
        self.canditdateGeneScores = {}

    def find_gene_locus(self, gene):
        for locus in self.loci:
            if gene in self.loci[locus]:
                return self.loci[locus], locus

    #
    def empty_locus_case(self, locus, subnet):
        subnet = subnet.copy()
        subnet.remove(locus)
        faUtilitiesInstance = FaUtilities(
            individualSubnetwork=subnet, parentNetworkFile=self.parentNetwork
        )
        edgeCount = faUtilitiesInstance.count_edges()
        return edgeCount

    def candidate_gene_score(self, locus, gene, subnet):
        # candidateGenesSCores = {"locusNumber":locusnumber, gene: , geneScore: }
        individualCandidateGeneScore = {}  # {gene,genescore}
        candidateGeneScores = []
        # {locusNumber: , genes:[individualCandidateGeneScore]}

        subnet = subnet.copy()
        locus = locus.copy()
        geneIndexInSubnet = subnet.index(gene)
        # print(f"GNE index: {geneIndexInSubnet}")

        faUtilitiesInstance = FaUtilities(
            individualSubnetwork=subnet, parentNetworkFile=self.parentNetwork
        )

        for item in locus:
            individualCandidateGeneScore = {}
            if item != gene:
                subnet[geneIndexInSubnet] = item
                # print(f"subnet after: {subnet}")
                faUtilitiesInstance.update_subnetwork(subnet)
                edgeCount = faUtilitiesInstance.count_edges()
                individualCandidateGeneScore["gene"] = item
                individualCandidateGeneScore["geneScore"] = edgeCount
                candidateGeneScores.append(individualCandidateGeneScore)
        return candidateGeneScores

    def gene_score(self):
        start = time.time()
        subnet = self.individualSubnet
        geneScores = {}
        print(f"score - genescore - subnet: {subnet}")
        faUtilitiesInstance = FaUtilities(
            individualSubnetwork=subnet, parentNetworkFile=self.parentNetwork
        )
        edgeCount = faUtilitiesInstance.count_edges()
        print(f"Subnet edgeCount: {edgeCount}")
        # 1 for each GENE in SUBNET
        for index, gene in enumerate(subnet):
            # 2 get empty locus score for each gene
            emptyLocusScore = self.empty_locus_case(gene, subnet)

            # print(f"emptyLocusScore: {emptyLocusScore}")

            # 3 find the locus and return list of locus genes
            locus, locusNumber = self.find_gene_locus(gene)
            # 4 create candidate gene scores for the x loops locus
            candidateGeneScores = self.candidate_gene_score(locus, gene, subnet)
            """for geneObj in candidateGeneScores:
                print(f"geneObj: {geneObj}")"""

            geneScores[locusNumber] = candidateGeneScores
            # print(f"candGeneScores: {candidateGeneScores}")
        end = time.time()
        print(f"gene score time: {end-start}")
        with open("gene.txt", "a") as file:
            for key, value in geneScores.items():
                file.write(str(key) + ": " + str(value) + "\n" + "\n")
