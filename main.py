import numpy as np
import pandas as pd
import time
import json
import itertools
import cProfile
import pstats
import concurrent.futures

import matplotlib.pyplot as plt
import networkx as nx
from netgraph import Graph

from components.module2_stage1_subnetworks import Stage1_SubNetworks
from components.score_individual_subnet import ScoreIndividualSubnet
from components.fa_utilities import FaUtilities
from components.module1_subnetwork import Module1Subnetwork

# Globally used objects:
# bins: an object that contains a key value pair for each edge count (calculated from the parentNetwork). key = edge count # : value = list of genes
# nonfaBins: an object derived from the bins object containing the same structure, but the list of genes only contain nonfa genes
# parentNetwork: a list of lists that contain a sublist per row in the STRING 1.txt file.


def average_gene_scores(geneScoresFile):
    geneScoresFromFile = []
    scoresByGene = {}
    averageGeneScores = {}
    with open(geneScoresFile, "r") as file:
        for line in file:
            locusId = line.split()[0][:-1]
            locus_str = " ".join(line.split()[1:])  # Join all parts after the first one
            locus_str = locus_str.replace(
                "'", '"'
            )  # Replace single quotes with double quotes
            locus = json.loads(locus_str)  # Now this should work
            for gene in locus:
                geneScoresFromFile.append({locusId: gene})

    # Group scores by gene
    scoresByGene = {}
    for score in geneScoresFromFile:
        locusId = ",".join(score.keys())  # Convert keys to a string
        score = list(score.values())[0]
        gene = score["gene"]
        if gene not in scoresByGene:
            scoresByGene[gene] = {"locusId": locusId, "scores": []}
            scoresByGene[gene]["locusId"] = locusId
        scoresByGene[gene]["scores"].append(score["geneScore"])

    for item in scoresByGene.items():
        gene = item[0]
        subitem = item[1]

        locusId = subitem["locusId"]
        scores = subitem["scores"]

        averageGeneScores[gene] = {
            "averageScore": sum(scores) / len(scores),
            "locusId": locusId,
        }
    seen = []
    for item in averageGeneScores:
        if item in seen:
            print(f"duplicate: {item}")
        elif item not in seen:
            seen.append(item)
    return averageGeneScores


def visualize_gene_scores(averageGeneScores, faNetworkFile):
    """df = pd.DataFrame({"gene1": ["A", "B", "C", "A"], "gene2": ["D", "A", "E", "C"]})

    graph = nx.from_pandas_edgelist(df, "gene1", "gene2")
    nx.draw(graph, with_labels=True)
    Graph(graph)"""
    faNetwork = []

    with open(faNetworkFile, "r") as file:
        for row in file:
            row = row.split()
            faNetwork.append(row)

    G = nx.Graph()

    for gene, data in averageGeneScores.items():
        # print(data["locusId"])
        G.add_node(gene, averageScore=data["averageScore"], locusId=data["locusId"])

    for edge in faNetwork:
        edgeOne = edge[0]
        edgeTwo = edge[1]
        G.add_edge(edgeOne, edgeTwo)

    color_map_dict = {
        "0": "blue",
        "1": "green",
        "2": "red",
        "3": "cyan",
        "4": "magenta",
        "5": "yellow",
        "6": "black",
        "7": "pink",
        "8": "brown",
        "9": "orange",
        "10": "purple",
        "11": "grey",
        # Add more colors if you have more locusIds
    }

    # Create a color map and size list based on locusId and averageScore
    color_map = []
    size_list = []
    for node in G:
        locusId = G.nodes[node]["locusId"]
        averageScore = G.nodes[node]["averageScore"]

        # if averageScore > 0:
        #    print((averageScore + 1) * 100)
        color_map.append(
            color_map_dict[locusId]
        )  # Use color map dictionary for color mapping
        size_list.append(
            (averageScore + 1) * 100
        )  # Add 1 to averageScore and multiply by 100 for visibility

    # Draw the graph
    nx.draw(G, node_color=color_map, node_size=size_list, with_labels=True)
    """fig, ax = plt.subplots()
    Graph(
        g,
        node_color=node_color,  # indicates the community each belongs to
        node_edge_width=0,  # no black border around nodes
        edge_width=0.1,  # use thin edges, as they carry no information in this visualisation
        edge_alpha=0.5,  # low edge alpha values accentuates bundles as they appear darker than single edges
        node_layout="community",
        node_layout_kwargs=dict(node_to_community=node_to_community),
        ax=ax,
    )"""
    plt.show()


def main():
    start = time.time()
    # Create parent network
    faUtilitiesInstance = FaUtilities(
        parentNetworkFile="STRING 1.txt", inputFile="Input.gmt.txt"
    )
    parentNetworkDict, parentNetworkDF = faUtilitiesInstance.create_parent_network()
    loci = faUtilitiesInstance.extract_loci()
    """module1_subnetworkInstance = Module1Subnetwork(
        parentNetworkDict, faGenes, nonfaGenes
    )
    module1_subnetwork = module1_subnetworkInstance.create_subnetwork()"""

    stage1_subnetworksInstance = Stage1_SubNetworks(
        "results.txt", "Input.gmt.txt", "STRING 1.txt", parentNetworkDF
    )

    stage1Subnetworks = stage1_subnetworksInstance.create_random_subnetworks()

    # Test Dataset
    testData = dict(itertools.islice(stage1Subnetworks.items(), 3))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []
        for subnet in testData.items():
            subnet = subnet[1]["subnet"]
            scoreIndividualSubnetInstance = ScoreIndividualSubnet(
                subnet, "Input.gmt.txt", parentNetworkDF, loci
            )
            futures.append(executor.submit(scoreIndividualSubnetInstance.gene_score))

        for future in concurrent.futures.as_completed(futures):
            geneScores = future.result()

    end = time.time()
    print(f"total time: {end - start}")
    averageGeneScores = average_gene_scores("gene.txt")
    visualize_gene_scores(averageGeneScores, "faNetwork.txt")


if __name__ == "__main__":
    main()

    """# cPROFILE
    cProfile.run("main()", "output.pstats")

    # Open a new text file in write mode
    with open("output.txt", "w") as f:
        # Create a pstats.Stats object with the text file as the stream
        stats = pstats.Stats("output.pstats", stream=f)

        # Sort the statistics by the cumulative time spent in the function
        stats.sort_stats("cumulative")

        # Print the statistics to the text file
        stats.print_stats()"""
