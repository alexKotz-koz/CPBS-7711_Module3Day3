import numpy as np
import time

import itertools
import cProfile
import pstats
import concurrent.futures


from components.module2_stage1_subnetworks import Stage1_SubNetworks
from components.score_individual_subnet import ScoreIndividualSubnet
from components.fa_utilities import FaUtilities
from components.module1_subnetwork import Module1Subnetwork

# Globally used objects:
# bins: an object that contains a key value pair for each edge count (calculated from the parentNetwork). key = edge count # : value = list of genes
# nonfaBins: an object derived from the bins object containing the same structure, but the list of genes only contain nonfa genes
# parentNetwork: a list of lists that contain a sublist per row in the STRING 1.txt file.


def main():
    start = time.time()
    # Create parent network
    faUtilitiesInstance = FaUtilities("STRING 1.txt")
    parentNetworkDict, parentNetworkDF = faUtilitiesInstance.create_parent_network()

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
        for subnet in stage1Subnetworks.items():
            subnet = subnet[1]["subnet"]
            scoreIndividualSubnetInstance = ScoreIndividualSubnet(
                subnet, "Input.gmt.txt", parentNetworkDF
            )
            futures.append(executor.submit(scoreIndividualSubnetInstance.gene_score))

        for future in concurrent.futures.as_completed(futures):
            geneScores = future.result()

    print(geneScores)
    end = time.time()
    print(f"total time: {end - start}")


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
