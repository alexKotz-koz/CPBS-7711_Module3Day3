from components.fa_utilities import FaUtilities
import time
import cProfile
import numpy as np


class ScoreIndividualSubnet:
    def __init__(self, individualSubnet, inputFile, parentNetwork):
        self.individualSubnet = individualSubnet
        self.inputFile = inputFile
        self.parentNetwork = parentNetwork
        faUtilitiesInstance = FaUtilities(inputFile=inputFile)
        self.loci = faUtilitiesInstance.extract_loci()
        self.canditdateGeneScores = {}

    def count_edges(self, subnet):
        # Convert subnetGenes to a set for faster membership tests
        start = time.time()
        subnetGenes = set(subnet)

        mask = self.parentNetwork["gene1"].isin(subnetGenes) & self.parentNetwork[
            "gene2"
        ].isin(subnetGenes)

        selectedRows = self.parentNetwork[mask].copy()

        # Use vectorized operations to create the sorted_genes column
        selectedRows["sorted_genes"] = np.sort(
            selectedRows[["gene1", "gene2"]], axis=1
        ).tolist()

        selectedRows.drop_duplicates(subset="sorted_genes", inplace=True)

        edgeCount = len(selectedRows)
        end = time.time()
        print(f"countEdgesTime: {end-start}")
        return edgeCount

    def find_gene_locus(self, gene):
        start = time.time()
        for locus in self.loci:
            if gene in self.loci[locus]:
                end = time.time()
                print(f"findgenelocus: {end-start}")
                return self.loci[locus], locus

    #
    def empty_locus_case(self, locus, subnet):
        start = time.time()
        subnet = subnet.copy()
        subnet.remove(locus)
        """faUtilitiesInstance = FaUtilities(
            individualSubnetwork=subnet, parentNetworkFile=self.parentNetwork
        )"""
        print("empty locus case count edges")
        edgeCount = self.count_edges(subnet)
        end = time.time()
        print(f"emptylocus: {end-start}")
        return edgeCount

    def candidate_gene_score(self, locus, gene, subnet):
        start = time.time()

        geneIndexInSubnet = subnet.index(gene)

        individualCandidateGeneScore = {"gene": None, "geneScore": None}
        candidateGeneScores = []

        for item in locus:
            if item != gene:
                subnet[geneIndexInSubnet] = item

                print("candidate gene score count edges")
                edgeCount = self.count_edges(subnet)
                individualCandidateGeneScore["gene"] = item
                individualCandidateGeneScore["geneScore"] = edgeCount
                candidateGeneScores.append(individualCandidateGeneScore.copy())

        end = time.time()
        print(f"candidategenescore: {end-start}")
        return candidateGeneScores

    def gene_score(self):
        start = time.time()
        subnet = self.individualSubnet
        geneScores = {}
        print(f"score - genescore - subnet: {subnet}")

        """faUtilitiesInstance = FaUtilities(
            individualSubnetwork=subnet, parentNetworkFile=self.parentNetwork
        )
        edgeCount = faUtilitiesInstance.count_edges()
        print(f"Subnet edgeCount: {edgeCount}")"""
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
