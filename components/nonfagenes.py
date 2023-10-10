from components.fagenes import FaGenes


class NonFaGenes:
    def __init__(self, inputFile, faGenes):
        self.inputFile = inputFile
        self.faGenes = faGenes

    def extract_nonfa_genes(self):
        nonfaGenesSet = set()
        nonfaGenes = {}
        faGenes = set(self.faGenes)

        with open(inputFile, "r") as file:
            parentNetwork = [row.split("\t") for row in file]

        for row in parentNetwork:
            if row[0] not in faGenes:
                if row[0] not in nonfaGenes:
                    nonfaGenes[row[0]] = 0
                nonfaGenes[row[0]] += 1
            if row[1] not in faGenes:
                if row[1] not in nonfaGenes:
                    nonfaGenes[row[1]] = 0
                nonfaGenes[row[1]] += 1

        return nonfaGenes
