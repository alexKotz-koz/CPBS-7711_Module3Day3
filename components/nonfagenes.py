from components.fagenes import FaGenes


class NonFaGenes:
    def __init__(self, inputFile, faGenes):
        self.inputFile = inputFile
        self.faGenes = faGenes

    # Input: faGenes object (all FA genes) and the STRING 1.txt file
    # Output: object containing all nonFA genes
    def extract_nonfa_genes(self):
        print("Creating non-FA genes...")
        nonfaGenesSet = set()
        nonfaGenes = {}
        faGenes = set(self.faGenes)

        with open(self.inputFile, "r") as file:
            parentNetwork = [row.split("\t") for row in file]

        # iterate over parentNetwork and add genes that do not exist in faGenes object, to nonfaGenes object
        for row in parentNetwork:
            if row[0] not in faGenes:
                if row[0] not in nonfaGenes:
                    nonfaGenes[row[0]] = 0
                nonfaGenes[row[0]] += 1
            if row[1] not in faGenes:
                if row[1] not in nonfaGenes:
                    nonfaGenes[row[1]] = 0
                nonfaGenes[row[1]] += 1
        print("Non-FA genes created")
        return nonfaGenes
