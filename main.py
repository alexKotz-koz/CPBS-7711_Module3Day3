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
import community

from components.module2_stage1_subnetworks import Stage1_SubNetworks
from components.score_individual_subnet import ScoreIndividualSubnet
from components.fa_utilities import FaUtilities
from components.module1_subnetwork import Module1Subnetwork


def average_gene_scores(geneScoresFile):
    geneScoresFromFile = []
    scoresByGene = {}
    averageGeneScores = {}
    # read in averageGeneScores file
    with open(geneScoresFile, "r") as file:
        for line in file:
            locusId = line.split()[0][:-1]
            locus_str = " ".join(line.split()[1:])
            locus_str = locus_str.replace("'", '"')
            locus = json.loads(locus_str)
            # create gene score dictionary for each gene
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

        # print(f"LEN OF GENE SCORES: {gene}|{len(scores)}")

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
    faNetwork = []

    # get edges from filtered parent network
    with open(faNetworkFile, "r") as file:
        for row in file:
            row = row.strip().split("\t")[:2]
            faNetwork.append(row)

    G = nx.Graph()
    # add nodes to graph from averageGeneScores and add score and locusId as attributes
    for gene, data in averageGeneScores.items():
        for row in faNetwork:
            if gene in row:
                G.add_node(
                    gene, averageScore=data["averageScore"], locusId=data["locusId"]
                )
    # add edges to graph object from the filtered parent network object
    nodesList = list(G.nodes)
    for edge in faNetwork:
        edgeOne = edge[0]
        edgeTwo = edge[1]
        # print(f"edgeOne:{edgeOne}")
        if edgeOne in nodesList and edgeTwo in nodesList:
            G.add_edge(edgeOne, edgeTwo)

    commDict = {}
    for index, node in enumerate(list(G.nodes)):
        commDict[node] = G.nodes[node]["locusId"]

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
    nodeSize = {node: (G.nodes[node]["averageScore"]) for node in list(G.nodes)}
    nodeColor = {node: color_map_dict[locusId] for node, locusId in commDict.items()}
    nodeLabels = {node: node for node in G.nodes}
    fig, ax = plt.subplots(figsize=(10, 8))
    Graph(
        G,
        node_color=nodeColor,  # indicates the community each belongs to
        node_size=nodeSize,
        node_edge_width=0,  # no black border around nodes
        edge_width=0.1,  # use thin edges, as they carry no information in this visualisation
        edge_alpha=0.5,  # low edge alpha values accentuates bundles as they appear darker than single edges
        node_layout="community",
        node_layout_kwargs=dict(node_to_community=commDict),
        node_labels=nodeLabels,
        node_label_size=20,
        ax=ax,
    )
    # Draw the nodes
    """nx.draw_networkx_nodes(
        G,
        pos=nx.spring_layout(G),
        node_color=list(nodeColor.values()),
        node_size=list(nodeSize.values()),
        ax=ax,
    )

    # Draw the edges
    nx.draw_networkx_edges(G, pos=nx.spring_layout(G), width=0.1, alpha=0.5, ax=ax)

    # Draw the labels
    nx.draw_networkx_labels(
        G, pos=nx.spring_layout(G), labels=nodeLabels, font_size=20, ax=ax
    )"""

    plt.show()


def main():
    '''start = time.time()
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
    print(f"total time: {end - start}")'''
    averageGeneScores = average_gene_scores("gene.txt")
    visualize_gene_scores(averageGeneScores, "results.txt")


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
