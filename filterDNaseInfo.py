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

def filterDNaseInfo(DNaseInfoFileName,nonRepGenesFileName, nonChangeGenesFileName, regulatorsFileName, nonTFsOutFileName, TFsOutFileName, numsNonTFsFileName, numsTFsFileName, probesetidGeneListFileName):
	# Iterate through genes, remove those that are non-reproducible and non-changing, separate TFs from non-TFs, and output DNase info. and the numbers of DNase for each gene to a file
	DNaseInfoFile = open(DNaseInfoFileName)
	probesetidGeneListFile = open(probesetidGeneListFileName)
	probesetidList = []
	geneList = []
	probesetidGeneListFile.readline()
	for line in probesetidGeneListFile:
		# Create a list of the probeset IDs and a list of the genes
		lineElements = line.split("\t")
		probesetidList.append(int(lineElements[0].strip()))
		geneList.append(lineElements[1].strip())
	nonRepGenesList = makeTFList(nonRepGenesFileName)
	nonChangeGenesList = makeTFList(nonChangeGenesFileName)
	regulatorsList = makeTFList(regulatorsFileName)
	nonTFsOutFile = open(nonTFsOutFileName, 'w+')
	TFsOutFile = open(TFsOutFileName, 'w+')
	numsNonTFsFile = open(numsNonTFsFileName, 'w+')
	numsTFsFile = open(numsTFsFileName, 'w+')
	for line in DNaseInfoFile:
		# Iterate through the DNase information for each gene and keep the information for genes that do not need to be filtered
		lineElements = line.split("\t")
		probesetID = int(lineElements[0])
		geneName = string.upper(lineElements[1].strip())
		if geneName not in geneList:
			# Immgen excluded this gene, so do not include it
			continue
		if probesetID not in probesetidList:
			# Immgen excluded this probeset ID, so do not include it
			print "No ID!"
			print probesetID
			print geneName
			continue
		probesetIDLocation = probesetidList.index(probesetID)
		if geneList[probesetIDLocation] != geneName:
			# Probeset ID and gene name do not match, so do not include them
			print "No match!"
			print probesetID
			print geneName
			continue
		if (geneName in nonRepGenesList) or (geneName in nonChangeGenesList):
			# The gene is not reproducible, so do not include it
			continue
		numDNaseSites = len(lineElements) - 5
		if geneName in regulatorsList:
			# The gene is a regulator, so put its information in the files with the regulators' information
			TFsOutFile.write(line)
			numsTFsFile.write(geneName)
			numsTFsFile.write("\t")
			numsTFsFile.write(str(numDNaseSites))
			numsTFsFile.write("\n")
		else:
			nonTFsOutFile.write(line)
			numsNonTFsFile.write(geneName)
			numsNonTFsFile.write("\t")
			numsNonTFsFile.write(str(numDNaseSites))
			numsNonTFsFile.write("\n")
	DNaseInfoFile.close()
	nonTFsOutFile.close()
	TFsOutFile.close()
	numsNonTFsFile.close()
	numsTFsFile.close()
	probesetidGeneListFile.close()

if __name__=="__main__":
    import sys
    import string
    DNaseInfoFileName = sys.argv[1]
    nonRepGenesFileName = sys.argv[2]
    nonChangeGenesFileName = sys.argv[3]
    regulatorsFileName = sys.argv[4]
    nonTFsOutFileName = sys.argv[5]
    TFsOutFileName = sys.argv[6]
    numsNonTFsFileName = sys.argv[7]
    numsTFsFileName = sys.argv[8]
    probsetidGeneListFileName = sys.argv[9]
                        
    filterDNaseInfo(DNaseInfoFileName, nonRepGenesFileName, nonChangeGenesFileName, regulatorsFileName, nonTFsOutFileName, TFsOutFileName, numsNonTFsFileName, numsTFsFileName, probsetidGeneListFileName)
