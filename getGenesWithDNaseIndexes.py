def getGenesWithDNaseIndexes(genesToDNaseFileName, genesWithDNaseFileName):
	# Get the indexes indicating which genes have at least 1 DNase site near-by
	# INDEXES ARE 1-INDEXED (NOT 0-INDEXED)
	genesToDNaseFile = open(genesToDNaseFileName)
	genesWithDNaseFile = open(genesWithDNaseFileName, 'w+')
	geneCount = 0
	for line in genesToDNaseFile:
		# Iterate through genes and record the indexes of those that have near-by DNase sites
		geneCount = geneCount + 1
		lineElements = line.strip().split("\t")
		if len(lineElements) > 4:
			# The gene is associated with at least 1 DNase site
			genesWithDNaseFile.write(str(geneCount))
			genesWithDNaseFile.write("\n")
	genesToDNaseFile.close()
	genesWithDNaseFile.close()

if __name__=="__main__":
    import sys
    genesToDNaseFileName = sys.argv[1] 
    genesWithDNaseFileName = sys.argv[2]                            
    getGenesWithDNaseIndexes(genesToDNaseFileName, genesWithDNaseFileName)
