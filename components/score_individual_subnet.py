from components.fa_utilities import FaUtilities
import time
import cProfile
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures


class ScoreIndividualSubnet:
    def __init__(self, individualSubnet, inputFile, parentNetwork):
        self.individualSubnet = individualSubnet
        self.inputFile = inputFile
        self.parentNetwork = parentNetwork
        faUtilitiesInstance = FaUtilities(inputFile=inputFile)
        self.loci = faUtilitiesInstance.extract_loci()
        self.canditdateGeneScores = {}

    def count_edges(self, subnets, emptyLocusScore, batch_size=60):
        candidateGeneScores = []
        start = time.time()
        with concurrent.futures.ThreadPoolExecutor() as executor:
            batched_subnets = []

            for i in range(0, len(subnets), batch_size):
                batched_subnets.append(subnets[i : i + batch_size])

            futures = {
                executor.submit(self.process_subnet, batch, emptyLocusScore)
                for batch in batched_subnets
            }

        for future in concurrent.futures.as_completed(futures):
            try:
                individualCandidateGeneScore = future.result()
                candidateGeneScores.extend(individualCandidateGeneScore)
            except Exception as exc:
                print(f"Count Edges - Generated an exception: {exc} {emptyLocusScore}")
        end = time.time()
        # print(f"countEdgesTime: {end - start}")
        return candidateGeneScores

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

    def process_subnet(self, subnet, emptyLocusScore):
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

    def find_gene_locus(self, gene):
        for locus in self.loci:
            if gene in self.loci[locus]:
                return self.loci[locus], locus

    #
    def empty_locus_case(self, geneLocus, subnet):
        subnet = subnet.copy()
        subnet.remove(geneLocus)

        edgeCount = self.process_empty_locus_case(geneLocus, subnet)

        return edgeCount

    def candidate_gene_score(self, locus, gene, subnet, emptyLocusScore):
        swappedSubnets = []
        geneIndexInSubnet = subnet.index(gene)

        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = {
                executor.submit(
                    self.process_item, item, geneIndexInSubnet, subnet, emptyLocusScore
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

    def process_item(self, item, geneIndexInSubnet, subnet, emptyLocusScore):
        subnet[geneIndexInSubnet] = item
        tempSubnet = subnet.copy()
        return {item: tempSubnet}

    def process_gene(self, gene, subnet):
        # 2 get empty locus score for each gene
        emptyLocusScore = self.empty_locus_case(gene, subnet)

        # 3 find the locus and return list of locus genes
        locus, locusNumber = self.find_gene_locus(gene)

        candidateGeneScores = self.candidate_gene_score(
            locus, gene, subnet, emptyLocusScore
        )

        return locusNumber, candidateGeneScores

    def gene_score(self):
        start = time.time()
        subnet = self.individualSubnet
        geneScores = {}

        with ProcessPoolExecutor() as executor:
            futures = {
                executor.submit(self.process_gene, gene, subnet): gene
                for gene in subnet
            }
            for future in concurrent.futures.as_completed(futures):
                locusNumber, candidateGeneScores = future.result()
                geneScores[locusNumber] = candidateGeneScores

        end = time.time()
        print(f"gene score time: {end-start}")
        with open("gene.txt", "a") as file:
            for key, value in geneScores.items():
                file.write(str(key) + ": " + str(value) + "\n")
        return geneScores
