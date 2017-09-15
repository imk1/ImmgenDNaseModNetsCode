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

def filterDNaseInfoUCSC(DNaseInfoFileName,nonRepGenesFileName, multiLocGenesFileName, nonChangeGenesFileName, allOutFileName, numsAllFileName, probesetidGeneListFileName):
	# Iterate through genes, remove those that are non-reproducible, have multiple loci, or are non-changing, and output DNase info. and the numbers of DNase for each gene to a file
	DNaseInfoFile = open(DNaseInfoFileName)
	probesetidGeneListFile = open(probesetidGeneListFileName)
	geneList = []
	probesetidGeneListFile.readline()
	for line in probesetidGeneListFile:
		# Create a list of the probeset IDs and a list of the genes
		lineElements = line.split("\t")
		geneList.append(string.upper(lineElements[1].strip()))
	nonRepGenesList = makeTFList(nonRepGenesFileName)
	multiLocGenesList = makeTFList(multiLocGenesFileName)
	nonChangeGenesList = makeTFList(nonChangeGenesFileName)
	allOutFile = open(allOutFileName, 'w+')
	numsAllFile = open(numsAllFileName, 'w+')
	for line in DNaseInfoFile:
		# Iterate through the DNase information for each gene and keep the information for genes that do not need to be filtered
		lineElements = line.split("\t")
		geneName = string.upper(lineElements[0].strip())
		if geneName not in geneList:
			# Immgen excluded this gene, so do not include it
			continue
		if ((geneName in nonRepGenesList) or (geneName in multiLocGenesList)) or (geneName in nonChangeGenesList):
			# The gene is not reproducible, so do not include it
			continue
		numDNaseSites = len(lineElements) - 4
		allOutFile.write(line)
		numsAllFile.write(geneName)
		numsAllFile.write("\t")
		numsAllFile.write(str(numDNaseSites))
		numsAllFile.write("\n")
	DNaseInfoFile.close()
	allOutFile.close()
	numsAllFile.close()
	probesetidGeneListFile.close()

if __name__=="__main__":
    import sys
    import string
    DNaseInfoFileName = sys.argv[1]
    nonRepGenesFileName = sys.argv[2]
    multiLocGenesFileName = sys.argv[3]
    nonChangeGenesFileName = sys.argv[4]
    allOutFileName = sys.argv[5]
    numsAllFileName = sys.argv[6]
    probsetidGeneListFileName = sys.argv[7]
                        
    filterDNaseInfoUCSC(DNaseInfoFileName, nonRepGenesFileName, multiLocGenesFileName, nonChangeGenesFileName, allOutFileName, numsAllFileName, probsetidGeneListFileName)
