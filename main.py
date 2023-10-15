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


def create_individual_nonfa_subnetwork(
    subnet, nonfaBin, bins, parentNetwork, subnetworksFromStage1
):
    return result


def create_secondary_subnetwork(
    parentNetworkFile, stage1Subnetworks, nonfaBin, bins, faGenes
):
    print("Creating stage 2 random subnetworks")

    stage2Subnetwork = {}
    parentNetwork = []

    with open(parentNetworkFile, "r") as file:
        parentNetwork = [row.split("\t")[:2] for row in file]

    try:
        results = []
        pool = multiprocessing.Pool(len(stage1Subnetworks.items()))
        for index, subnet in stage1Subnetworks.items():
            subnetworksFromStage1 = subnet["subnet"]
            print(len(subnetworksFromStage1))

            # Create an instance of the Create_Individual_Nonfa_Subnetwork_Thread class
            instance = Create_Individual_Nonfa_Subnetwork_Thread(
                subnet, nonfaBin, bins, parentNetwork, subnetworksFromStage1
            )
            thread = instance.run()
            # Start the thread
            result = pool.apply_async(thread)
            # Append the thread and AsyncResult object to the results list
            results.append((thread, result))

        for result in results:
            result.wait()

        # Close the pool
        pool.close()
        pool.join()
w
        for index, (thread, result) in enumerate(results):
            # Get the result from the thread object
            thread.get_result()
            stage2Subnetwork[index] = thread.result
            edgeCount = stage2Subnetwork[index]["edgeCount"]
            subnet = stage2Subnetwork[index]["subnet"]
            faGeneBinFlag = stage2Subnetwork[index]["faGeneBinFlag"]
            binNotFoundFlag = stage2Subnetwork[index]["binNotFoundFlag"]
            print(edgeCount, subnet)
    except Exception as e:
        print(f"Error:{e}")

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
