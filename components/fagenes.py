class FaGenes:
    def __init__(self, inputFile):
        self.inputFile = inputFile

    def fanconi_anemia_genes(self):
        faData = []
        faGenes = []
        with open(self.inputFile, "r") as file:
            faData = [row.strip().split("\t")[2:] for row in file]
        for row in faData:
            faGenes.extend(row)
        return faGenes
