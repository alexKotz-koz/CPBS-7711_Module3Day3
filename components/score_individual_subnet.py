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

    def count_edges(self, subnets, emptyLocusScore, emptyLocusSubnet):
        # Convert subnetGenes to a set for faster membership tests
        candidateGeneScores = []
        print(f"Count Edges - emptyLocusScore {emptyLocusScore}")
        start = time.time()

        for subnet in subnets:
            individualCandidateGeneScore = {"gene": None, "geneScore": None}
            subnetGenes = [item for sublist in subnet.values() for item in sublist]
            gene = subnet.keys()

            # print(f"SUBNETGENES: {subnetGenes}")

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

            if emptyLocusScore > 0 and edgeCount < emptyLocusScore:
                print(
                    f"**********empty locus score: {emptyLocusScore} | edge count: {edgeCount} | subnet:{subnetGenes} | subnetEmpty:{emptyLocusSubnet}"
                )

            individualCandidateGeneScore["gene"] = gene
            individualCandidateGeneScore["geneScore"] = edgeCount - emptyLocusScore
            candidateGeneScores.append(individualCandidateGeneScore)
        end = time.time()
        print(f"countEdgesTime: {end-start}")
        return candidateGeneScores

    def find_gene_locus(self, gene):
        # start = time.time()
        for locus in self.loci:
            if gene in self.loci[locus]:
                # end = time.time()
                # print(f"findgenelocus: {end-start}")
                return self.loci[locus], locus

    #
    def empty_locus_case(self, geneLocus, subnet):
        # start = time.time()
        subnet = subnet.copy()
        subnet.remove(geneLocus)
        faUtilitiesInstance = FaUtilities(
            individualSubnetwork=subnet, parentNetworkFile=self.parentNetwork
        )
        # print("empty locus case count edges")
        edgeCount = faUtilitiesInstance.count_edges()
        emptyLocusSubnet = subnet
        # end = time.time()
        # print(f"emptylocus: {end-start}")
        return edgeCount, emptyLocusSubnet

    def candidate_gene_score(
        self, locus, gene, subnet, emptyLocusScore, emptyLocusSubnet
    ):
        # start = time.time()

        geneIndexInSubnet = subnet.index(gene)
        subnet = subnet.copy()
        swappedSubnets = []

        for item in locus:
            # test every gene in locus, including gi
            subnet[geneIndexInSubnet] = item
            tempSubnet = subnet.copy()
            # print(f"subnet: {tempSubnet}")
            swappedSubnets.append({item: tempSubnet})
        candidateGeneScores = self.count_edges(
            swappedSubnets, emptyLocusScore, emptyLocusSubnet
        )
        # print(f"subnets swapped: {candidateGeneScores}")

        # end = time.time()
        # print(f"candidategenescore: {end-start}")
        return candidateGeneScores

    def gene_score(self):
        start = time.time()
        subnet = self.individualSubnet
        geneScores = {}
        print(f"score - genescore - subnet: {subnet}")

        # 1 for each GENE in SUBNET
        for gene in subnet:
            # 2 get empty locus score for each gene
            emptyLocusScore, emptyLocusSubnet = self.empty_locus_case(gene, subnet)

            # print(f"emptyLocusScore: {emptyLocusScore}")

            # 3 find the locus and return list of locus genes
            locus, locusNumber = self.find_gene_locus(gene)
            # 4 create candidate gene scores for the x loops locus
            candidateGeneScores = self.candidate_gene_score(
                locus, gene, subnet, emptyLocusScore, emptyLocusSubnet
            )

            geneScores[locusNumber] = candidateGeneScores

        end = time.time()
        print(f"gene score time: {end-start}")
        with open("gene.txt", "a") as file:
            for key, value in geneScores.items():
                file.write(str(key) + ": " + str(value) + "\n")
