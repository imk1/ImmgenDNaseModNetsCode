def createDNaseFeatGeneMatrix(genesToDNaseFileName, genesDNaseFeatFileName, numDNaseSites, numGenes):
	# Create a d x g matrix, where d is the number of DNase sites and g is the number of genes, of DNase sites and the values for the features for the genes that they are near
	# -1 means that there is not feature for the gene, either because the DNase site is not near the gene's TSS or because the information for that feature was not available
	genesToDNaseFile = open(genesToDNaseFileName)
	genesDNaseFeatFile = open(genesDNaseFeatFileName)
	DNaseFeatGeneMat = []
	for i in range(numDNaseSites):
		# Add a row to the feature matrix for each DNase site
		if i % 10000 == 0:
			print i
		DNaseFeatGeneMatRow = []
		for j in range(numGenes):
			# DNaseFeatGeneMat all of the features to -1
			DNaseFeatGeneMatRow.append(-1)
		DNaseFeatGeneMat.append(DNaseFeatGeneMatRow)
	geneCount = 0
	for line in genesToDNaseFile:
		# Iterate through genes and modify their feature values
		if geneCount % 10000 == 0:
			print geneCount
		lineElements = line.split("\t")
		featLine = genesDNaseFeatFile.readline()
		featLineElements = featLine.split("\t")
		DNaseCount = 0
		for DNaseIndexStr in lineElements[4:len(lineElements)]:
			# Iterate through the DNase sites for the current gene and modify the features appropriately
			DNaseIndex = int(DNaseIndexStr)
			featVal = float(featLineElements[DNaseCount+1])
			DNaseFeatGeneMat[DNaseIndex][geneCount] = featVal
			DNaseCount = DNaseCount + 1
		geneCount = geneCount + 1
	genesToDNaseFile.close()
	genesDNaseFeatFile.close()
	return DNaseFeatGeneMat

def writeDNaseFeatGeneMat(DNaseFeatGeneMat, outputFileName):
	# Write the DNase-gene pairwise feature matrix to a file
	print "Writing to file!"
	outputFile = open(outputFileName, 'w+')
	for DNaseSiteFeats in DNaseFeatGeneMat:
		# Iterate through rows of DNase site features and write each row to the file
		for DNaseGeneFeat in DNaseSiteFeats:
			# Iterate through features for each gene and write each feature to the file
			outputFile.write(str(DNaseGeneFeat))
			outputFile.write("\t")
		outputFile.write("\n")
	outputFile.close()

if __name__=="__main__":
    import sys
    genesToDNaseFileName = sys.argv[1]
    genesDNaseFeatFileName = sys.argv[2]
    numDNaseSites = int(sys.argv[3])
    numGenes = int(sys.argv[4])
    outputFileName = sys.argv[5]
                        
    DNaseFeatGeneMat = createDNaseFeatGeneMatrix(genesToDNaseFileName, genesDNaseFeatFileName, numDNaseSites, numGenes)
    writeDNaseFeatGeneMat(DNaseFeatGeneMat, outputFileName)
