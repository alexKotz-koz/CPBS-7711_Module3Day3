import json


class Bins:
    def __init__(self, inputFile):
        self.inputFile = inputFile

    def create_bins(self):
        countPerGene = {}
        bins = {}
        with open(self.inputFile, "r") as file:
            results = [row.split("\t")[:2] for row in file]

        for row in results:
            if row[0] not in countPerGene:
                countPerGene[row[0]] = 0
            elif row[0] in countPerGene:
                countPerGene[row[0]] += 1

            if row[1] not in countPerGene:
                countPerGene[row[1]] = 0
            elif row[1] in countPerGene:
                countPerGene[row[1]] += 1

        sorted_dict = dict(sorted(countPerGene.items(), key=lambda item: item[1]))

        for item in sorted_dict:
            bins.setdefault(sorted_dict[item], []).append(item)

        with open("bins.json", "w") as outputFile:
            json.dump(bins, outputFile)

        return bins
