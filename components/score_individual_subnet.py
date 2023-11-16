from components.fa_utilities import FaUtilities
import time
import cProfile
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures


class ScoreIndividualSubnet:
    def __init__(self, individualSubnet, inputFile, parentNetwork, loci):
        self.individualSubnet = individualSubnet
        self.inputFile = inputFile
        self.parentNetwork = parentNetwork
        self.loci = loci
        self.canditdateGeneScores = {}

    # Stage 8
    def count_edges(self, subnets, emptyLocusScore, batch_size=60):
        candidateGeneScores = []

        with concurrent.futures.ThreadPoolExecutor() as executor:
            batched_subnets = []

            for i in range(0, len(subnets), batch_size):
                batched_subnets.append(subnets[i : i + batch_size])

            futures = {
                executor.submit(self.process_subnet_count_edges, batch, emptyLocusScore)
                for batch in batched_subnets
            }

        for future in concurrent.futures.as_completed(futures):
            try:
                individualCandidateGeneScore = future.result()
                candidateGeneScores.extend(individualCandidateGeneScore)
            except Exception as exc:
                print(f"Count Edges - Generated an exception: {exc} {emptyLocusScore}")
        return candidateGeneScores

    # Stage 9
    def process_subnet_count_edges(self, subnet, emptyLocusScore):
        batchScores = []
        for item in subnet:
            gene, subnetGenes = list(item.items())[0]

            mask = self.parentNetwork["gene1"].isin(subnetGenes) & self.parentNetwork[
                "gene2"
            ].isin(subnetGenes)

            selectedRows = self.parentNetwork[mask].copy()

            selectedRows["sorted_genes"] = np.sort(
                selectedRows[["gene1", "gene2"]], axis=1
            ).tolist()
            selectedRows.drop_duplicates(subset=["sorted_genes"], inplace=True)

            edgeCount = len(selectedRows)

            individualCandidateGeneScore = {
                "gene": gene,
                "geneScore": edgeCount - emptyLocusScore,
            }
            batchScores.append(individualCandidateGeneScore)
        return batchScores

    # Stage 3
    def empty_locus_case(self, geneLocus, subnet):
        subnet = subnet.copy()
        subnet.remove(geneLocus)

        edgeCount = self.process_empty_locus_case(geneLocus, subnet)

        return edgeCount

    # Stage 4
    def process_empty_locus_case(self, gene, subnet):
        mask = self.parentNetwork["gene1"].isin(subnet) & self.parentNetwork[
            "gene2"
        ].isin(subnet)
        selectedRows = self.parentNetwork[mask].copy()

        selectedRows["sorted_genes"] = np.sort(
            selectedRows[["gene1", "gene2"]], axis=1
        ).tolist()
        selectedRows.drop_duplicates(subset=["sorted_genes"], inplace=True)

        edgeCount = len(selectedRows)
        return edgeCount

    # Stage 5
    def find_gene_locus(self, gene):
        for locus in self.loci:
            if gene in self.loci[locus]:
                return self.loci[locus], locus

    # Stage 6
    def candidate_gene_score(self, locus, gene, subnet, emptyLocusScore):
        swappedSubnets = []
        # get the index of the gene in the subnet list
        geneIndexInSubnet = subnet.index(gene)

        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = {
                executor.submit(
                    self.process_item_candidate_gene_score,
                    item,
                    geneIndexInSubnet,
                    subnet,
                )
                for item in locus
            }

        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
                swappedSubnets.append(result)
            except Exception as exc:
                print(f"Candidate Gene Score - Generated an exception: {exc}")

        candidateGeneScores = self.count_edges(swappedSubnets, emptyLocusScore)
        return candidateGeneScores

    # Stage 7
    def process_item_candidate_gene_score(self, item, geneIndexInSubnet, subnet):
        subnet[geneIndexInSubnet] = item
        tempSubnet = subnet.copy()
        return {item: tempSubnet}

    # Stage 2
    # Input: Gene from subnet and subnet
    # Output:
    def process_gene_gene_score(self, gene, subnet):
        # get empty locus score for each gene
        estart = time.time()
        emptyLocusScore = self.empty_locus_case(gene, subnet)
        eend = time.time()
        print(f"Empty Locus Time: {eend - estart}")

        # find the locus and return list of locus genes
        lstart = time.time()
        locus, locusNumber = self.find_gene_locus(gene)
        lend = time.time()
        print(f"Find Locus Time: {lend-lstart}")

        # get candidate gene score for each gene in the subnet
        cstart = time.time()
        candidateGeneScores = self.candidate_gene_score(
            locus, gene, subnet, emptyLocusScore
        )
        cend = time.time()
        print(f"Candidate Gene Score Time: {cend-cstart}")

        return locusNumber, candidateGeneScores

    # Stage 1
    # Input: Individual subnet
    # Output: Average Gene Scores object
    def gene_score(self):
        print(f"Subnet Scoring Initialized for subnet: {self.individualSubnet}")
        start = time.time()
        subnet = self.individualSubnet
        geneScores = {}

        # for each gene in the subnet, call process_gene which handles the intialization logic of scoring the genes within the genes' locus
        with ProcessPoolExecutor() as executor:
            futures = {
                executor.submit(self.process_gene_gene_score, gene, subnet): gene
                for gene in subnet
            }
            # gather results of candidate_gene_scores, as they complete and store in geneScores
            for future in concurrent.futures.as_completed(futures):
                locusNumber, candidateGeneScores = future.result()
                geneScores[locusNumber] = candidateGeneScores

        end = time.time()
        print(f"gene score time: {end-start}")
        with open("gene.txt", "a") as file:
            for key, value in geneScores.items():
                file.write(str(key) + ": " + str(value) + "\n")
        return geneScores
