def findGenesWithDNaseSites(genesToDNaseFileName, outputFileName):
	# Finds the genes with near-by DNase sites and writes their indexs in the list to a file
	genesToDNaseFile = open(genesToDNaseFileName)
	outputFile = open(outputFileName, 'w+')
	lineCount = 0
	for line in genesToDNaseFile:
		# Iterate through the genes and record the indexes of those that have near-by DNase sites
		# LINES ARE 1-INDEXED because output will be used in Matlab
		lineCount = lineCount + 1
		lineElements = line.strip().split("\t")
		if len(lineElements) > 4:
			# There is at least 1 DNase hypersensitivity site near the current gene, so write the gene's index
			outputFile.write(str(lineCount))
			outputFile.write("\n")
	genesToDNaseFile.close()
	outputFile.close()

if __name__=="__main__":
    import sys
    genesToDNaseFileName = sys.argv[1]
    outputFileName = sys.argv[2]
    findGenesWithDNaseSites(genesToDNaseFileName, outputFileName)
