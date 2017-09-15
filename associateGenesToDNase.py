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

def writeGeneDNasesToFile(outputFile, probesetID, geneName, chrom, start, end, DNaseForGeneIndexes):
    # Write a gene, its location, and the indexes of the DNase hypersensitivity sites near its TSS to the output file
    outputFile.write(probesetID)
    outputFile.write("\t")
    outputFile.write(geneName)
    outputFile.write("\t")
    outputFile.write(chrom)
    outputFile.write("\t")
    outputFile.write(str(start))
    outputFile.write("\t")
    outputFile.write(str(end))
    for di in DNaseForGeneIndexes:
	# Write each DNase index to the file
	outputFile.write("\t")
	outputFile.write(str(di))
    outputFile.write("\n")

def getDNaseInfoSpecific(DNaseInfo, DNaseIndex):
    # Gets the information about a specific DNase hypersensitivity site
    DNaseElements = DNaseInfo[DNaseIndex].split("\t")
    DNaseChrom = DNaseElements[0]
    DNaseStart = int(DNaseElements[1])
    DNaseEnd = int(DNaseElements[2])
    DNaseMiddle = DNaseStart + (float(DNaseEnd - DNaseStart) / float(2))
    return [DNaseChrom, DNaseMiddle]

def associateGenesToDNase(genesFileName, DNaseFileName, outputFileName, cutoff):
    # Find all of the DNase hypersensitivity sites that are within the cutoff bases of each gene's TSS
    genesFile = open(genesFileName)
    outputFile = open(outputFileName, 'w+')
    genesInfo = genesFile.readlines()
    genesFile.close()
    DNaseFile = open(DNaseFileName)
    DNaseInfo = DNaseFile.readlines()
    DNaseFile.close()
    outputFile = open(outputFileName, 'w+')
    lastGenesIndex = -1
    chrom = ""
    start = 0
    end = 0
    probesetID = 0
    geneName = ""
    genesIndex = 0
    DNaseIndex = 0
    [DNaseChrom, DNaseMiddle] = getDNaseInfoSpecific(DNaseInfo, DNaseIndex)
    DNaseForGeneIndexes = []
    while genesIndex < len(genesInfo):
        # Find the DNase peaks associated with each gene
	if genesIndex > lastGenesIndex:
		# Get the information for the new gene
        	geneElements = genesInfo[genesIndex].split(",")
		probesetID = geneElements[0]
		if len(geneElements) <= 7:
			# No name is given for the gene
			genesIndex = genesIndex + 1
			continue
		chrom = geneElements[2][1:len(geneElements[2])-1]
		start = int(geneElements[4])
		end = int(geneElements[5])
        	genesElementsWithName = geneElements[7]
		genesInfoElements = genesElementsWithName.split("//")
		if len(genesInfoElements) == 1:
			# No name is given for the gene
			genesIndex = genesIndex + 1
			continue
		geneName = genesInfoElements[1]
		if len(DNaseForGeneIndexes) > 0:
			# Go back to the 1st DNase for the previous gene
			DNaseIndex = DNaseForGeneIndexes[0]
			[DNaseChrom, DNaseMiddle] = getDNaseInfoSpecific(DNaseInfo, DNaseIndex)
    			DNaseForGeneIndexes = []
		lastGenesIndex = genesIndex
	allDNasesUsed = False
	while DNaseChrom < chrom:
		# Iterate through DNase sites until a site on the correct chromosome has been found
		DNaseIndex = DNaseIndex + 1
		if DNaseIndex >= len(DNaseInfo):
			# All DNases have been used, so stop
			allDNasesUsed = True
			break
		[DNaseChrom, DNaseMiddle] = getDNaseInfoSpecific(DNaseInfo, DNaseIndex)
	if allDNasesUsed == True:
		# No more DNases remain, so write the current gene to the file
		writeGeneDNasesToFile(outputFile, probesetID, geneName, chrom, start, end, DNaseForGeneIndexes)
		genesIndex = genesIndex + 1
		continue			
	while (DNaseMiddle < start - cutoff) and (DNaseChrom == chrom):
		# Iterate through DNase sites until a site near the gene's TSS has been found
		DNaseIndex = DNaseIndex + 1
		if DNaseIndex >= len(DNaseInfo):
			# All DNases have been used, so stop
			allDNasesUsed = True
			break
		[DNaseChrom, DNaseMiddle] = getDNaseInfoSpecific(DNaseInfo, DNaseIndex)
	if allDNasesUsed == True:
		# No more DNases remain, so write the current gene to the file
		writeGeneDNasesToFile(outputFile, probesetID, geneName, chrom, start, end, DNaseForGeneIndexes)
		genesIndex = genesIndex + 1
		continue	
	if DNaseChrom > chrom:
		# There are no more DNase hypersensitivity sites on the chromosome with the gene
		writeGeneDNasesToFile(outputFile, probesetID, geneName, chrom, start, end, DNaseForGeneIndexes)
		genesIndex = genesIndex + 1
		continue
	if DNaseMiddle > start + cutoff:
		# There are no more DNase hypersensitivity sites near the gene's TSS
		writeGeneDNasesToFile(outputFile, probesetID, geneName, chrom, start, end, DNaseForGeneIndexes)
		genesIndex = genesIndex + 1
		continue
	DNaseForGeneIndexes.append(DNaseIndex)
	DNaseIndex = DNaseIndex + 1
	if DNaseIndex >= len(DNaseInfo):
		# All DNases have been used, so stop
		writeGeneDNasesToFile(outputFile, probesetID, geneName, chrom, start, end, DNaseForGeneIndexes)
		genesIndex = genesIndex + 1
		continue
	[DNaseChrom, DNaseMiddle] = getDNaseInfoSpecific(DNaseInfo, DNaseIndex)
    outputFile.close()

if __name__=="__main__":
    import sys
    import string
    genesFileName = sys.argv[1]
    DNaseFileName = sys.argv[2]
    outputFileName = sys.argv[3]
    cutoff = int(sys.argv[4])
                        
    associateGenesToDNase(genesFileName, DNaseFileName, outputFileName, cutoff)
