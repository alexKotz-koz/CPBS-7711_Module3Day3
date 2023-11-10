import random
import json
from math import sqrt
import numpy as np

import itertools
import cProfile
import pstats

from components.bins import Bins
from components.fagenes import FaGenes
from components.nonfagenes import NonFaGenes
from components.module2_stage1_subnetworks import Stage1_SubNetworks
from components.score_individual_subnet import ScoreIndividualSubnet
from components.fa_utilities import FaUtilities
from components.module1_subnetwork import Module1Subnetwork

# Globally used objects:
# bins: an object that contains a key value pair for each edge count (calculated from the parentNetwork). key = edge count # : value = list of genes
# nonfaBins: an object derived from the bins object containing the same structure, but the list of genes only contain nonfa genes
# parentNetwork: a list of lists that contain a sublist per row in the STRING 1.txt file.


def main():
    faGenesInstance = FaGenes("Input.gmt.txt")
    faGenes = faGenesInstance.fanconi_anemia_genes()

    nonfaGenesInstance = NonFaGenes("STRING 1.txt", faGenes=faGenes)
    nonfaGenes = nonfaGenesInstance.extract_nonfa_genes()

    # Create parent network
    faUtilitiesInstance = FaUtilities("STRING 1.txt")
    parentNetworkDict, parentNetworkDF = faUtilitiesInstance.create_parent_network()

    """module1_subnetworkInstance = Module1Subnetwork(
        parentNetworkDict, faGenes, nonfaGenes
    )
    module1_subnetwork = module1_subnetworkInstance.create_subnetwork()"""

    """binsInstance = Bins(parentNetworkDict, "results.txt", faGenes, nonfaGenes)
    bins, nonfaBins = binsInstance.create_bins()"""

    stage1_subnetworksInstance = Stage1_SubNetworks(
        "results.txt", "Input.gmt.txt", "STRING 1.txt", parentNetworkDF
    )

    stage1Subnetworks = stage1_subnetworksInstance.create_random_subnetworks()

    # Test Dataset
    testData = dict(itertools.islice(stage1Subnetworks.items(), 3))

    for subnet in testData.items():
        subnet = subnet[1]["subnet"]
        """faUtilitiesInstance = FaUtilities(parentNetworkDF, subnet, "Input.gmt.txt")
        edgeCount = faUtilitiesInstance.count_edges()
        print("edge: ", edgeCount)"""
        scoreIndividualSubnetInstance = ScoreIndividualSubnet(
            subnet, "Input.gmt.txt", parentNetworkDF
        )
        geneScores = scoreIndividualSubnetInstance.gene_score()


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
