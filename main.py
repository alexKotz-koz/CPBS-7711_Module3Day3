import random
import json
from math import sqrt
import numpy as np
import threading
import concurrent.futures
from components.bins import Bins
from components.fagenes import FaGenes
from components.nonfagenes import NonFaGenes
from components.stage1_subnetworks import Stage1_SubNetworks
from components.create_individual_nonfa_subnetwork_thread import (
    Create_Individual_Nonfa_Subnetwork_Thread,
)
from scipy.stats import permutation_test
from scipy.stats import norm


# constructor type function for use in create_secondary_subnetwork.
# this function initializes an instance of the Create_Individual_Nonfa_Subnetwork_Thread class as a thread and executes the classes run function
# Input: nonfaBin structure (containing a bin that represents the number of edge connections and the content of the bin are strictly nonfa genes)
# Ouput: the returned dictionary containing the individual random nonfa subnetwork, two flags used in the create_secondary_network function, and the edge count of the subnetwork
def create_individual_nonfa_subnetwork(nonfaBin, bins, parentNetwork, subnet):
    thread = Create_Individual_Nonfa_Subnetwork_Thread(
        nonfaBin, bins, parentNetwork, subnet
    )

    result = thread.run()
    return result


# this function generates 5000 nanfa subnetworks (supporting the null hypothesis)
# Input: parentNetworkFile = STRING 1.txt, nonfaBin structure, bins structure (similar to the nonfaBin structure, but includes all genes from parentNetwork), and the faGenes object
# Output: a structure containing all 5000 random nonfa subnetworks
def create_secondary_subnetwork(
    parentNetworkFile, stage1Subnetworks, nonfaBin, bins, faGenes
):
    print("Creating stage 2 random subnetworks")

    stage2Subnetworks = {}
    parentNetwork = []

    with open(parentNetworkFile, "r") as file:
        parentNetwork = [row.split("\t")[:2] for row in file]

    # test for individual subnetwork
    """for index, subnet in stage1Subnetworks.items():
        print(f"index{index}")

        instance = Create_Individual_Nonfa_Subnetwork_Thread(
            nonfaBin, bins, parentNetwork, subnet["subnet"]
        )
        result = instance.run()
        stage2Subnetworks[index] = result
        print(result)
    with open("stage2_random_subnetworks.json", "w") as outputFile:
        json.dump(stage2Subnetworks, outputFile)
    print("Second 5,000 subnetworks created")"""

    # create a list of threads for each subnet
    threads = [
        Create_Individual_Nonfa_Subnetwork_Thread(
            nonfaBin, bins, parentNetwork, subnet["subnet"]
        )
        for subnet in stage1Subnetworks.values()
    ]

    # create a thread pool to expidite the creation of individual null subnets
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []
        # test: chunk out 50 stage1 random fa subnetworks
        """stage1Subnetworks = dict(
            list(stage1Subnetworks.items())[: len(stage1Subnetworks) // 100]
        )"""

        # submit each subnetwork to the thread pool
        for index, subnet in stage1Subnetworks.items():
            future = executor.submit(
                create_individual_nonfa_subnetwork,
                nonfaBin,
                bins,
                parentNetwork,
                subnet["subnet"],
            )
            # add each running thread to a list for use in gathering results of all threads
            futures.append(future)

        # get the results from the thread pool and add to greater stage2Subnetworks object
        for index, future in enumerate(concurrent.futures.as_completed(futures)):
            result = future.result()
            stage2Subnetworks[index] = result

    # write stage2Subnetworks object for varification and raw visualization of data
    with open("stage2_random_subnetworks.json", "w") as outputFile:
        json.dump(stage2Subnetworks, outputFile)
    print("Second 5,000 subnetworks created")

    return stage2Subnetworks
    """for index, future in enumerate(
        concurrent.futures.as_completed(future for _, future in futures)
    ):"""
    """# Split the stage 1 subnetworks into batches
    subnetwork_groups = [
        dict(list(stage1Subnetworks.items())[i : i + 50])
        for i in range(0, len(stage1Subnetworks), 50)
    ]

    # Create a list of threads for each batch
    threads = []
    for subnetwork_group in subnetwork_groups:
        for subnet in subnetwork_group.values():
            subnet = subnet["subnet"]
            thread = Create_Individual_Nonfa_Subnetwork_Thread(
                nonfaBin, bins, parentNetwork, subnet
            )
            threads.append(thread)

    # Create a multiprocessing pool
    pool = multiprocessing.Pool()

    # Execute the threads in batches using the pool
    batch_size = 10
    batches = [threads[i : i + batch_size] for i in range(0, len(threads), batch_size)]
    batch_results = []
    for batch in batches:
        batch_result = pool.map(execute_thread, batch)
        batch_results.extend(batch_result)

    # Combine the results from each batch
    stage2Subnetworks = {}
    for result in batch_results:
        subnetIndex = result["index"]
        stage2Subnetworks[subnetIndex] = result
    # Close the pool
    pool.close()
    pool.join()"""


# set up function for p_test
def statistic(x, y, axis):
    return np.mean(x, axis=axis) - np.mean(y, axis=axis)


def p_test(stage1SubnetworksFile, stage2SubnetworksFile):
    stage1SubnetworksEdgeCount = []
    stage2SubnetworksEdgeCount = []
    alpha = 0.05

    with open(stage1SubnetworksFile, "r") as file:
        stage1Subnetworks = json.load(file)
    for bin in stage1Subnetworks:
        stage1SubnetworksEdgeCount.append(stage1Subnetworks[bin]["edgeCount"])
    stage1SubnetworksMean = sum(stage1SubnetworksEdgeCount) / 5000

    with open(stage2SubnetworksFile, "r") as file:
        stage2Subnetworks = json.load(file)

    totalEdgeCount = 0
    for bin in stage2Subnetworks:
        binEdgeCount = 0
        subnet = stage2Subnetworks[bin]["subnet"]
        for item in subnet:
            if isinstance(item, list):
                totalEdgeCount += 1
                binEdgeCount += 1
        if binEdgeCount >= 1:
            stage2SubnetworksEdgeCount.append(binEdgeCount)
        else:
            stage2SubnetworksEdgeCount.append(0)

    stage2SubnetworksMean = totalEdgeCount / 5000
    rng = np.random.default_rng()

    x = stage1SubnetworksEdgeCount
    y = stage2SubnetworksEdgeCount

    statistic(x, y, 0)

    res = permutation_test(
        (x, y), statistic, vectorized=True, n_resamples=9999, alternative="less"
    )
    print(f"pval: {res.pvalue}")


def main():
    faGenesInstance = FaGenes("Input.gmt.txt")
    faGenes = faGenesInstance.fanconi_anemia_genes()

    nonfaGenesInstance = NonFaGenes("STRING 1.txt", faGenes=faGenes)
    nonfaGenes = nonfaGenesInstance.extract_nonfa_genes()

    binsInstance = Bins("STRING 1.txt", "results.txt", faGenes, nonfaGenes)
    bins, nonfaBins = binsInstance.create_bins()

    stage1_subnetworksInstance = Stage1_SubNetworks(
        "results.txt", "Input.gmt.txt", "STRING 1.txt"
    )

    (
        stage1Subnetworks,
        edgeCount,
    ) = stage1_subnetworksInstance.create_random_subnetworks()
    stage2_subnetworks = create_secondary_subnetwork(
        "STRING 1.txt", stage1Subnetworks, nonfaBins, bins, faGenes
    )

    pVal = p_test("stage1_random_subnetworks.json", "stage2_random_subnetworks.json")
    print(pVal)


if __name__ == "__main__":
    main()
