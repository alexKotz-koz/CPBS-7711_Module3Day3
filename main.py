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

    G = Graph()
    # for each of the genes in the averageGeneScores object, create a node that contains the locusId and averageScore as attributes.
    for gene, data in averageGeneScores.items():
        G.add_node(gene, averageScore=data["averageScore"], locusId=data["locusId"])

    # add edges to graph object from the filtered parent network object
    nodesList = list(G.nodes)
    for edge in faNetwork:
        edgeOne = edge[0]
        edgeTwo = edge[1]
        # print(f"edgeOne:{edgeOne}")
        if edgeOne in nodesList:
            if edgeTwo in nodesList:
                G.add_edge(edgeOne, edgeTwo)
        elif edgeTwo in nodesList:
            if edgeOne in nodesList:
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
    for item in G:
        node = G.nodes[item]

        if node == {}:
            G.remove_node(item)
            break

        locusId = node["locusId"]
        averageScore = node["averageScore"]

        color_map.append(color_map_dict[locusId])
        size_list.append((averageScore + 1) * 100)

    # Step 1: Create a dictionary to store the positions of the nodes
    pos = {}

    # Step 2: Get a list of unique locusId values
    locusIds = set(nx.get_node_attributes(G, "locusId").values())

    # Step 3: For each locusId, create a subgraph containing only the nodes with that locusId
    for locusId in locusIds:
        # Get the nodes with this locusId
        nodes = [
            node for node, attr in G.nodes(data=True) if attr["locusId"] == locusId
        ]

        # Create a subgraph with these nodes
        subgraph = G.subgraph(nodes)

        # Step 4: Apply a layout algorithm to the subgraph
        subgraph_pos = nx.spring_layout(subgraph)

        # Step 5: Add the positions of the nodes in the subgraph to the positions dictionary
        pos.update(subgraph_pos)

    # Create a node_color dictionary
    node_color = {
        node: color_map_dict[data["locusId"]] for node, data in G.nodes(data=True)
    }
    nx.draw(
        G,
        pos=pos,
        node_color=color_map,
        node_size=size_list,
        with_labels=True,
    )

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
