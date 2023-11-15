import pandas as pd
import numpy as np
import time


class FaUtilities:
    def __init__(
        self, parentNetworkFile=None, individualSubnetwork=None, inputFile=None
    ):
        self.individualSubnetwork = individualSubnetwork

        self.inputFile = inputFile

        if isinstance(parentNetworkFile, pd.DataFrame):
            self.parentNetwork = parentNetworkFile
        elif isinstance(parentNetworkFile, list):
            self.parentNetwork = parentNetworkFile
        elif isinstance(parentNetworkFile, str):
            self.parentNetworkFile = parentNetworkFile
            self.parentNetwork = pd.DataFrame()
        else:
            self.parentNetworkFile = parentNetworkFile

    def filter_parent_network(self):
        # Creating a filtered Parent Network: Only contains FA to FA connections. Limitation.
        print("Filtering Parent Network for FA Genes")
        faLoci = self.extract_loci()
        faGenes = [string for sublist in faLoci.values() for string in sublist]
        faNetwork = []
        with open(self.parentNetworkFile, "r") as file:
            for line in file:
                line = line.strip().split("\t")[:2]
                if line[0] in faGenes:
                    if line[1] in faGenes:
                        faNetwork.append(line)
                if line[1] in faGenes:
                    if line[0] in faGenes:
                        faNetwork.append(line)
        with open("faNetwork.txt", "w") as file:
            file.write("\n".join("\t".join(sublist) for sublist in faNetwork))
        print("FA Network Created")
        return faNetwork

    def create_parent_network(self):
        print("Creating Parent Network")
        start = time.time()
        self.filter_parent_network()
        # substitute for "faNetwork.txt"
        self.parentNetwork = pd.read_csv(
            "faNetwork.txt",
            sep="\t",
            header=None,
            names=["gene1", "gene2"],
            usecols=[0, 1],
        )

        # Convert gene1 and gene2 to strings
        self.parentNetwork["gene1"] = self.parentNetwork["gene1"].astype(str)
        self.parentNetwork["gene2"] = self.parentNetwork["gene2"].astype(str)

        # Create a set of sorted gene pairs
        sorted_gene_pairs = set(
            map(tuple, np.sort(self.parentNetwork[["gene1", "gene2"]].values, axis=1))
        )

        # Filter the DataFrame based on the set of sorted gene pairs
        self.parentNetwork = self.parentNetwork[
            self.parentNetwork.apply(
                lambda row: tuple(sorted([row["gene1"], row["gene2"]]))
                in sorted_gene_pairs,
                axis=1,
            )
        ]
        parentNetworkDict = self.parentNetwork.to_dict("records")
        end = time.time()
        ex = end - start
        print(f"Parent Network finished in : {ex}")
        return parentNetworkDict, self.parentNetwork

    def count_edges(self):
        """if isinstance(self.individualSubnetwork, dict):
            subnetGenes = self.individualSubnetwork[1]["subnet"]
        else:
            subnetGenes = self.individualSubnetwork"""

        """print(
            f"fa utilities - edge count - empty locus case - subnet: {len(subnetGenes)}{subnetGenes}"
        )"""
        # Convert subnetGenes to a set for faster membership tests
        subnetGenes = self.individualSubnetwork

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

        return edgeCount

    def extract_loci(self):
        loci = {}
        print("Extracting FA Loci")
        with open(self.inputFile, "r") as file:
            for line in file:
                # line = line.split()
                name = line.split()
                loci[name[3]] = line.strip().split("\t")[2:]
        print("Loci Extracted")
        return loci
        # print(f"loci: {loci} \n")
