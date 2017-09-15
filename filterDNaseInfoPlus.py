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

def filterDNaseInfoPlus(DNaseInfoFileName,nonRepGenesFileName, nonChangeGenesFileName, outFileName, probesetidGeneListFileName):
	# Iterate through genes, remove those that are non-reproducible and non-changing, and output DNase info.
	DNaseInfoFile = open(DNaseInfoFileName)
	probesetidGeneListFile = open(probesetidGeneListFileName)
	geneList = []
	probesetidGeneListFile.readline()
	for line in probesetidGeneListFile:
		# Create a list of the probeset IDs and a list of the genes
		lineElements = line.split("\t")
		geneList.append(lineElements[1].strip())
	nonRepGenesList = makeTFList(nonRepGenesFileName)
	nonChangeGenesList = makeTFList(nonChangeGenesFileName)
	outFile = open(outFileName, 'w+')
	for line in DNaseInfoFile:
		# Iterate through the DNase information for each gene and keep the information for genes that do not need to be filtered
		lineElements = line.split("\t")
		geneName = string.upper(lineElements[0].strip())
		if geneName not in geneList:
			# Immgen excluded this gene, so do not include it
			continue
		if (geneName in nonRepGenesList) or (geneName in nonChangeGenesList):
			# The gene is not reproducible, so do not include it
			continue
		outFile.write(line)
	DNaseInfoFile.close()
	outFile.close()
	probesetidGeneListFile.close()

if __name__=="__main__":
    import sys
    import string
    DNaseInfoFileName = sys.argv[1]
    nonRepGenesFileName = sys.argv[2]
    nonChangeGenesFileName = sys.argv[3]
    outFileName = sys.argv[4]
    probesetidGeneListFileName = sys.argv[5]
                        
    filterDNaseInfoPlus(DNaseInfoFileName,nonRepGenesFileName, nonChangeGenesFileName, outFileName, probesetidGeneListFileName)
