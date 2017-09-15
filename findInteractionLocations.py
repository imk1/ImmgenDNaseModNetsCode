def findInteractionLocations(geneListFileName, interactionFileName, interactionLocationFileName):
	# Finds the locations of the genes in the identified interactions in the gene list and creates a file with 2 columns, 1 with the TFs row number and the other with the other gene's row number
	geneListFile = open(geneListFileName)
	interactionFile = open(interactionFileName)
	interactionLocationFile = open(interactionLocationFileName, 'w+')
	geneList = []
	for line in geneListFile:
		# Iterate through genes and add them to the list
		geneList.append(line.strip())
	geneListFile.close()
	for line in interactionFile:
		# Iterate through interactions and find the locations of the corresponding genes in the list
		lineElements = line.split("\t")
		TFLocation = geneList.index(lineElements[0].strip())
		geneLocation = geneList.index(lineElements[1].strip())
		interactionLocationFile.write(str(TFLocation))
		interactionLocationFile.write("\t")
		interactionLocationFile.write(str(geneLocation))
		interactionLocationFile.write("\n")
	interactionFile.close()
	interactionLocationFile.close()

if __name__=="__main__":
    import sys
    geneListFileName = sys.argv[1]
    interactionFileName = sys.argv[2]
    interactionLocationFileName = sys.argv[3]
    findInteractionLocations(geneListFileName, interactionFileName, interactionLocationFileName)
