import random
import json
from math import sqrt
import numpy as np
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

# Globally used objects:
# bins: an object that contains a key value pair for each edge count (calculated from the parentNetwork). key = edge count # : value = list of genes
# nonfaBins: an object derived from the bins object containing the same structure, but the list of genes only contain nonfa genes
# parentNetwork: a list of lists that contain a sublist per row in the STRING 1.txt file.


# constructor type function for use in create_secondary_subnetwork.
# this function initializes an instance of the Create_Individual_Nonfa_Subnetwork_Thread class as a thread and executes the classes run function
# Input: nonfaBin, bins, STRING 1.txt parentNetwork, stage1 subnetwork
# Ouput: the returned dictionary containing the individual random nonfa subnetwork, two flags used in the create_secondary_network function, and the edge count of the subnetwork
def create_individual_nonfa_subnetwork(nonfaBin, bins, parentNetwork, subnet):
    thread = Create_Individual_Nonfa_Subnetwork_Thread(
        nonfaBin, bins, parentNetwork, subnet
    )

    result = thread.run()
    return result


# this function generates 5000 nanfa subnetworks
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

    # create a list of threads, each thread is an instance of a nonfa subnet
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


# set up function for p_test
def statistic(x, y, axis):
    return np.mean(x, axis=axis) - np.mean(y, axis=axis)


def p_test(stage1SubnetworksFile, stage2SubnetworksFile):
    stage1SubnetworksEdgeCount = []
    stage2SubnetworksEdgeCount = []
    alpha = 0.05

    # find the mean edge count of both stage 1 subnetworks and stage 2 subnetworks
    with open(stage1SubnetworksFile, "r") as file:
        stage1Subnetworks = json.load(file)
    for bin in stage1Subnetworks:
        stage1SubnetworksEdgeCount.append(stage1Subnetworks[bin]["edgeCount"])
    stage1SubnetworksMean = sum(stage1SubnetworksEdgeCount) / 5000
    print(f"stage1Mean: {stage1SubnetworksMean}")

    with open(stage2SubnetworksFile, "r") as file:
        stage2Subnetworks = json.load(file)
    for bin in stage2Subnetworks:
        stage2SubnetworksEdgeCount.append(stage2Subnetworks[bin]["edgeCount"])
    stage2SubnetworksMean = sum(stage2SubnetworksEdgeCount) / 5000
    print(f"stage2Mean: {stage2SubnetworksMean}")

    # calculate difference in means
    observedStatistic = stage2SubnetworksMean - stage1SubnetworksMean

    numPermutations = 10000

    permStatistic = np.zeros(numPermutations)

    # perform the permutation test
    for i in range(numPermutations):
        # combine the data and shuffle it
        combinedSubnetworks = np.concatenate(
            (stage1SubnetworksEdgeCount, stage2SubnetworksEdgeCount)
        )
        np.random.shuffle(combinedSubnetworks)

        # calculate the test statistic for this permutation
        permGroup1 = combinedSubnetworks[: len(stage1SubnetworksEdgeCount)]
        permGroup2 = combinedSubnetworks[len(stage1SubnetworksEdgeCount) :]
        permStat = np.mean(permGroup2) - np.mean(permGroup1)

        # ctore the permuted statistic
        permStatistic[i] = permStat

    # calculate the p-value by comparing the observed statistic to the permuted statistics
    pVal = (np.abs(permStatistic) >= np.abs(observedStatistic)).mean()

    return pVal

    """    ic("Observed Statistic:", observedStatistic)
        print("abs of permStatistic", np.abs(permStatistic))
        ic(permStatistic >= observedStatistic)
        ic(p_value)
    """
    """rng = np.random.default_rng()

    x = stage1SubnetworksEdgeCount
    y = stage2SubnetworksEdgeCount
    x = norm.rvs(size=len(stage1SubnetworksEdgeCount), random_state=rng)
    y = norm.rvs(size=len(stage2SubnetworksEdgeCount), loc=3, random_state=rng)
    print(statistic(x, y, 0))

    res = permutation_test(
        (x, y),
        statistic,
        vectorized=True,
        n_resamples=9999,
        alternative="less",
        random_state=rng,
    )"""
    # print(f"pval: {res.pvalue}")


def main():
    """faGenesInstance = FaGenes("Input.gmt.txt")
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
    )"""

    pVal = p_test(
        "stage1_random_subnetworks.json",
        "stage2_random_subnetworks.json",
    )

    print(f"P-Value: {pVal}")


if __name__ == "__main__":
    main()
