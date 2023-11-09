from components.fa_utilities import FaUtilities


class Module1Subnetwork:
    def __init__(self, parentNetwork, faGenes, nonFaGenes):
        self.parentNetwork = parentNetwork
        self.faGenes = faGenes
        self.nonFaGenes = nonFaGenes

    def find_possible_nonFa_connection(self, parentNetworkRow):
        faGenes = self.faGenes
        nonFaGenes = self.nonFaGenes

        pass

    ##PRE: subnetwork, generated from create_subnetwork()
    ##POST: subnetwork that only contains unique rows. ***this functionality is to remove the duplicated edge connection between any given node-node pair.
    def check_duplicate(self, results):
        # using set's, as this functionality is related to uniqueness and speed is crucial for effecient completion of the script.
        print("Checking parent subnetwork for duplicates")
        final = set()
        seen = set()
        duplicates = set()

        for row in results:
            # sort each row alphabetically and store as tuple to make membership testing more efficent.
            orderedRow = tuple(sorted(row))

            # if the row from the results file has already 'seen', indicating the contents of the row [gene1, gene2, edge] are identical to another row (order non-specific), add to dummy set.abs
            # else add to seen and final results set
            if orderedRow in seen:
                duplicates.add(tuple(orderedRow))
            else:
                seen.add(orderedRow)
                final.add(tuple(row))

        # write unique rows to final subnetwork
        with open("results.txt", "w") as outputFile:
            for row in final:
                outputFile.write("\t".join(row) + "\n")
        outputFile.close()

    def create_subnetwork(self):
        print("Creating parent subnetwork")
        parentNetwork = self.parentNetwork
        faGenes = self.faGenes
        results = []

        # Add FA
        for row in parentNetwork:
            # if both genes are in faGenes
            if row["gene1"] in faGenes:
                if row["gene2"] in faGenes:
                    results.append(row.values())

            # if the first gene is in faGenes, but not the second. call fx to find possible nonfagene connection
            elif row["gene1"] in faGenes:
                if row["gene2"] in nonFaGenes:
                    print(f"nonFAGene:{row['gene2']}")

            elif row["gene2"] in faGenes:
                if row["gene1"] not in faGenes:
                    pass

        with open("results.txt", "w") as outputFile:
            for row in results:
                outputFile.write("\t".join(row) + "\n")

        self.check_duplicate(results)
