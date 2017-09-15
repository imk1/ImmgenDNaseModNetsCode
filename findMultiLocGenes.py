def makeTFList(TFFileName):
	# Make a list of TFs in which all of the letters are capitalized
	TFFile = open(TFFileName)
	TFList = []
	for line in TFFile:
		# Iterate through TF file and make a list of all of the TFs in upper case
		TFCap = string.upper(line.strip())
		TFList.append(TFCap)
	TFFile.close()
	return TFList

def findMultiLocGenes(sortedGenesFileName, multiLocGenesFileName):
	# Find genes that occur multiple times in a sorted list of genes
	multiLocGeneFile = open(multiLocGenesFileName, 'w+')
	sortedGenes = makeTFList(sortedGenesFileName)
	currentGene = sortedGenes[0]
	currentIndex = 0
	while currentIndex < len(sortedGenes) - 1:
		# Iterate through genes and write every gene that occurs multiple times to the file
		if sortedGenes[currentIndex + 1] == currentGene:
			# The gene is repeated
			multiLocGeneFile.write(currentGene)
			multiLocGeneFile.write("\n")
			while sortedGenes[currentIndex + 1] == currentGene:
				# Iterate through genes until a new gene has been reached
				currentIndex = currentIndex + 1
		currentIndex = currentIndex + 1
		currentGene = sortedGenes[currentIndex]

if __name__=="__main__":
    import sys
    import string
    sortedGenesFileName = sys.argv[1]
    multiLocGenesFileName = sys.argv[2]
    findMultiLocGenes(sortedGenesFileName, multiLocGenesFileName)
