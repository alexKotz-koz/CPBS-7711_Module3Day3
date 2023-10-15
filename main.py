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


def execute_thread(thread):
    thread.start()
    thread.join()
    return thread.result


def create_secondary_subnetwork(
    parentNetworkFile, stage1Subnetworks, nonfaBin, bins, faGenes
):
    print("Creating stage 2 random subnetworks")

    stage2Subnetwork = {}
    parentNetwork = []

    with open("STRING 1.txt", "r") as file:
        parentNetwork = [row.split("\t")[:2] for row in file]

    threads = []
    for index, subnet in stage1Subnetworks.items():
        subnetworksFromStage1 = subnet["subnet"]

        # print(f"subnetworkFromStage1: {subnetworksFromStage1}")

        # Create the thread object
        thread = Create_Individual_Nonfa_Subnetwork_Thread(
            subnet, nonfaBin, bins, parentNetwork, subnetworksFromStage1
        )
        threads.append(thread)
    print(len(threads))

    thread_batches = [threads[i : i + 1000] for i in range(0, len(threads), 1000)]

    for result in results:
        stage2Subnetwork[index] = result
        edgeCount = result["edgeCount"]
        subnet = result["subnet"]
        faGeneBinFlag = result["faGeneBinFlag"]
        binNotFoundFlag = result["binNotFoundFlag"]
        print(
            edgeCount,
            subnet,
        )
    """try:
        threads = []
        for index, subnet in stage1Subnetworks.items():
            subnetworksFromStage1 = subnet["subnet"]
            # Create Thread instance here
            create_individual_nonfa_subnetwork_thread_instance = (
                create_individual_nonfa_subnetwork_thread(
                    subnet, nonfaBin, bins, parentNetwork, stage1Subnetworks
                )
            )
            t = threading.Thread(
                target=create_individual_nonfa_subnetwork_thread_instance.run
            )
            threads.append(t)

        for thread in threads:
            thread.start()

        for thread in threads:
            result = thread.join()
            stage2Subnetwork[index] = result
        print(stage2Subnetwork)
    except Exception as e:
        print("Error:", e)"""

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
