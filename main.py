import random
import json
from math import sqrt
import numpy as np
import threading
import multiprocessing
from components.bins import Bins
from components.fagenes import FaGenes
from components.nonfagenes import NonFaGenes
from components.stage1_subnetworks import Stage1_SubNetworks
from components.create_individual_nonfa_subnetwork_thread import (
    Create_Individual_Nonfa_Subnetwork_Thread,
)
from scipy.stats import permutation_test
from scipy.stats import norm


def create_thread(subnet, nonfaBin, bins, parentNetwork):
    subnetworksFromStage1 = subnet["subnet"]
    # Create the thread object
    thread = Create_Individual_Nonfa_Subnetwork_Thread(
        subnet, nonfaBin, bins, parentNetwork, subnetworksFromStage1
    )
    # Return the thread object
    return thread


def execute_batch(batch):
    results = []
    print(batch)
    for thread in batch:
        print("executing thread")
        thread.start()
        thread.join()
        result = thread.get_result()
        results.append(result)
    return results


def create_secondary_subnetwork(
    parentNetworkFile, stage1Subnetworks, nonfaBin, bins, faGenes
):
    print("Creating stage 2 random subnetworks")

    stage2Subnetwork = {}
    parentNetwork = []

    with open(parentNetworkFile, "r") as file:
        parentNetwork = [row.split("\t")[:2] for row in file]
    #      Split the stage1Subnetworks dictionary into 5 groups
    subnetwork_groups = [
        dict(list(stage1Subnetworks.items())[i : i + 5])
        for i in range(0, len(stage1Subnetworks), 5)
    ]

    # Create a pool of processes
    num_processes = 5
    pool = multiprocessing.Pool(num_processes)

    # Submit each group of subnetworks to the pool
    for subnetwork_group in subnetwork_groups:
        results = []
        for subnet in subnetwork_group.values():
            # Create a new thread for each subnetwork
            thread = create_thread(subnet, nonfaBin, bins, parentNetwork)
            results.append(thread)
        # Execute the batch of threads
        batch_results = execute_batch(results)
        # Store the results in the stage2Subnetwork dictionary
        for index, batch_result in enumerate(batch_results):
            stage2Subnetwork[index] = batch_result
            edgeCount = batch_result["edgeCount"]
            subnet = batch_result["subnet"]
            faGeneBinFlag = batch_result["faGeneBinFlag"]
            binNotFoundFlag = batch_result["binNotFoundFlag"]
            print(edgeCount, subnet, faGeneBinFlag, binNotFoundFlag)

    # Close the pool
    pool.close()
    pool.join()

    with open("stage2_random_subnetworks.json", "w") as outputFile:
        json.dump(stage2Subnetwork, outputFile)
    print("Second 5,000 subnetworks created")


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

    print(f"Seed: {rng}")
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

    pVal = p_test(
        "stage1_random_subnetworks.json", "stage2_random_subnetworks copy.json"
    )
    print(pVal)


if __name__ == "__main__":
    main()
