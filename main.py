import numpy as np
import pandas as pd
import time
import json
import itertools
import concurrent.futures

import matplotlib.pyplot as plt
import networkx as nx
from netgraph import Graph
import community

from components.module2_stage1_subnetworks import Stage1_SubNetworks
from components.score_individual_subnet import ScoreIndividualSubnet
from components.fa_utilities import FaUtilities
from components.module1_subnetwork import Module1Subnetwork


# Input: geneScores text file (generated from score_individual_subnet.py), faNetwork text file (FA-FA connections)
# Output: averageGeneScores dictionary that contains one item per gene {geneScores list, locusId}
def average_gene_scores(geneScoresFile, faNetworkFile):
    geneScoresFromFile = []
    scoresByGene = {}
    averageGeneScores = {}
    faNetwork = []
    scoresByGene = {}

    # create fa network list
    with open(faNetworkFile, "r") as file:
        for line in file:
            line = line.split()
            faNetwork.append(line)

    # create gene scores dictionary
    with open(geneScoresFile, "r") as file:
        for line in file:
            locusId = line.split()[0][:-1]
            locus_str = " ".join(line.split()[1:])
            locus_str = locus_str.replace("'", '"')
            locus = json.loads(locus_str)
            for gene in locus:
                geneScoresFromFile.append({locusId: gene})

    # restructure geneScoresFromFile for use in calculating average gene scores
    for score in geneScoresFromFile:
        locusId = ",".join(score.keys())
        score = list(score.values())[0]
        gene = score["gene"]
        if gene not in scoresByGene:
            scoresByGene[gene] = {"locusId": locusId, "scores": []}
            scoresByGene[gene]["locusId"] = locusId
        scoresByGene[gene]["scores"].append(score["geneScore"])

    # calculate average gene scores and store in new dictionary
    for item in scoresByGene.items():
        gene = item[0]
        subitem = item[1]

        locusId = subitem["locusId"]
        scores = subitem["scores"]

        averageGeneScores[gene] = {
            "averageScore": sum(scores) / len(scores),
            "locusId": locusId,
        }

    # if a gene from averageGeneScores is not in the FA-FA network, mark genescore as "NA"
    for gene in averageGeneScores:
        if not any(gene in sublist for sublist in faNetwork):
            averageGeneScores[gene]["averageScore"] = "NA"

    return averageGeneScores


# Input: averageGeneScores dictionary, fa network text file
# Output: creation of nextworkx network graph
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
            if gene in row and data["averageScore"] != "NA":
                G.add_node(
                    gene, averageScore=data["averageScore"], locusId=data["locusId"]
                )

    # add edges to graph object from the filtered parent network object
    nodesList = list(G.nodes)
    for edge in faNetwork:
        edgeOne = edge[0]
        edgeTwo = edge[1]
        if edgeOne in nodesList and edgeTwo in nodesList:
            G.add_edge(edgeOne, edgeTwo)

    # create community dictionary to group nodes by locus id
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
        "12": "olive",
    }

    # create node size and color maps for visualization
    nodeSize = {node: (G.nodes[node]["averageScore"] * 10) for node in list(G.nodes)}
    nodeColor = {node: color_map_dict[locusId] for node, locusId in commDict.items()}
    fig, ax = plt.subplots(figsize=(10, 8))

    # trigger visualization
    Graph(
        G,
        node_color=nodeColor,
        node_size=nodeSize,
        node_edge_width=0.2,
        edge_width=0.1,
        edge_alpha=0.5,
        node_layout="community",
        node_layout_kwargs=dict(node_to_community=commDict),
        node_alpha=0.75,
        edge_layout="bundled",
        ax=ax,
    )

    plt.show()


def main():
    start = time.time()
    # create filtered parent network and fa loci
    faUtilitiesInstance = FaUtilities(
        parentNetworkFile="STRING 1.txt", inputFile="Input.gmt.txt"
    )
    parentNetworkDF = faUtilitiesInstance.create_parent_network()
    loci = faUtilitiesInstance.extract_loci()

    # create 5000 random fa subnetworks
    stage1_subnetworksInstance = Stage1_SubNetworks(
        "results.txt", "Input.gmt.txt", "STRING 1.txt", parentNetworkDF
    )
    stage1Subnetworks = stage1_subnetworksInstance.create_random_subnetworks()

    # Test Dataset: cut out the first 1000 fa subnetworks to reduce runtime...abs
    # THIS LINE TRUNCATES THE 5000 SUBNETWORKS TO REDUCE RUNTIME
    # TO TEST FULL FUNCTIONALITY: COMMENT THIS LINE AND REPLACE testData.items() WITH stage1Subnetworks.items() on line 170
    testData = dict(itertools.islice(stage1Subnetworks.items(), 3))

    # Process Pool to create gene scores for testData or stage1Subnetwork data
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []
        for subnet in testData.items():
            subnet = subnet[1]["subnet"]
            scoreIndividualSubnetInstance = ScoreIndividualSubnet(
                subnet, "Input.gmt.txt", parentNetworkDF, loci
            )
            futures.append(executor.submit(scoreIndividualSubnetInstance.gene_score))

        # as each process completes store the results in geneScores
        for future in concurrent.futures.as_completed(futures):
            geneScores = future.result()

    # calculate average gene scores and visualize the average gene scores for every FA gene that exists in the FA network
    averageGeneScores = average_gene_scores("gene.txt", "faNetwork.txt")
    visualize_gene_scores(averageGeneScores, "results.txt")

    end = time.time()
    print(f"total time: {end - start}")


if __name__ == "__main__":
    main()
