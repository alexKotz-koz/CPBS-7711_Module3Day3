import pandas as pd
import numpy as np


class FaUtilities:
    def __init__(self, parentNetworkFile=None, individualSubnetwork=None):
        self.individualSubnetwork = individualSubnetwork
        if isinstance(parentNetworkFile, pd.DataFrame):
            self.parentNetwork = parentNetworkFile
        elif isinstance(parentNetworkFile, str):
            self.parentNetworkFile = parentNetworkFile
            self.parentNetwork = pd.DataFrame()
        else:
            raise TypeError("parentNetworkFile must be a DataFrame or a string")

    def create_parent_network(self):
        self.parentNetwork = pd.read_csv(
            self.parentNetworkFile,
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
        return self.parentNetwork

    def count_edges(self):
        subnetGenes = self.individualSubnetwork[1]["subnet"]
        flatSubnetGenes = []

        for item in subnetGenes:
            if isinstance(item, list):
                flatSubnetGenes.extend(item)
            else:
                flatSubnetGenes.append(item)

        print(self.individualSubnetwork)
        print(flatSubnetGenes)

        mask = self.parentNetwork["gene1"].isin(flatSubnetGenes) & self.parentNetwork[
            "gene2"
        ].isin(flatSubnetGenes)

        # Use the mask to select the rows where both conditions are true
        selected_rows = self.parentNetwork[mask].copy()

        # Create a new column that contains a sorted tuple of gene1 and gene2
        selected_rows["sorted_genes"] = selected_rows.apply(
            lambda row: tuple(sorted([row["gene1"], row["gene2"]])), axis=1
        )

        # Drop duplicates based on the new column
        selected_rows.drop_duplicates(subset="sorted_genes", inplace=True)

        # Count the number of unique edges
        edgeCount = len(selected_rows)

        return edgeCount
