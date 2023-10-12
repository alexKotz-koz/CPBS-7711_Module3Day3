import random
import json
from components.bins import Bins
from components.fagenes import FaGenes
from components.nonfagenes import NonFaGenes
from components.stage1_subnetworks import Stage1_SubNetworks


def create_secondary_subnetwork(stage1Subnetworks, nonfaGenes, bins):
    print("Creating stage 2 random subnetworks")

    stage2Subnetwork = {}
    nonfaGenes = set(nonfaGenes.keys())

    for index, subnet in stage1Subnetworks.items():
        # print(subnet["subnet"])
        subnetworksFromStage1 = subnet["subnet"]
        tempFlattendSubnetwork = set()
        for gene in subnetworksFromStage1:
            if isinstance(gene, list):
                for subGene in gene:
                    tempFlattendSubnetwork.add(subGene)
            else:
                tempFlattendSubnetwork.add(gene)
        #print(tempFlattendSubnetwork)
        tempFlattendSubnetwork = list(tempFlattendSubnetwork)
        for gene in tempFlattendSubnetwork:
            
    """with open("stage2_random_subnetworks.json", "w") as outputFile:
        json.dump(stage2_randomSubnetworks, outputFile)
    print("Second 5,000 subnetworks created")"""


def t_test(stage1Subnetworks, stage2Subnetworks):
    pValue = 0

    return pValue


def main():
    ##TEST TO SEE HOW MANY ROWS ARE IN NONFAGENES
    """with open("STRING 1.txt", "r") as file:
    se = set()
    results = [row.split("\t")[:2] for row in file]
    for i in results:
        for j in i:
            se.add(j)
    print(len(se))"""

    binsInstance = Bins("STRING 1.txt")
    bins = binsInstance.create_bins()

    faGenesInstance = FaGenes("STRING 1.txt")
    faGenes = faGenesInstance.fanconi_anemia_genes()

    nonfaGenesInstance = NonFaGenes("STRING 1.txt", faGenes=faGenes)
    nonfaGenes = nonfaGenesInstance.extract_nonfa_genes()

    stage1_subnetworksInstance = Stage1_SubNetworks(
        "results.txt", "Input.gmt.txt", "STRING 1.txt"
    )

    (
        stage1Subnetworks,
        edgeCount,
    ) = stage1_subnetworksInstance.create_random_subnetworks()
    stage2_subnetworks = create_secondary_subnetwork(
        stage1Subnetworks, nonfaGenes, bins
    )
    # print(edgeCount)


if __name__ == "__main__":
    main()
