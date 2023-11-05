class FaGenes:
    def __init__(self, inputFile):
        self.inputFile = inputFile

    # Input: results.txt from module 1 day 3 homework, a subnetwork of only connected fa genes
    # Output: faGenes object that contains fa genes
    def fanconi_anemia_genes(self):
        print("Creating FA genes object...")
        faData = []
        faGenes = []
        with open(self.inputFile, "r") as file:
            faData = [row.strip().split("\t")[2:] for row in file]
        for row in faData:
            faGenes.extend(row)
        print("FA genes created")
        return faGenes
