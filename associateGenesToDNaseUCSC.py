def writeGeneDNasesToFileUCSC(outputFile, geneName, chrom, start, end, DNaseForGeneIndexes):
    # Write a gene, its location, and the indexes of the DNase hypersensitivity sites near its TSS to the output file
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

def associateGenesToDNaseUCSC(genesFileName, DNaseFileName, outputFileName, cutoff):
    # Find all of the DNase hypersensitivity sites that are within the cutoff bases of each gene's TSS
    # ASSUMES THAT TRANSCRIPTS FROM DIFFERENT MULTI-TRANSCRIPT GENES DO NOT OVERLAP, THAT ALL TRANSCRIPTS IN A GENE ARE CONSECUTIVE, AND THAT THERE ARE NO GENES WITH THE SAME NAME ON DIFFERENT CHROMOSOMES
    # DNase hypersensitivity sites are 0-INDEXED
    genesFile = open(genesFileName)
    genesFile.readline()
    genesInfo = genesFile.readlines()
    genesFile.close()
    DNaseFile = open(DNaseFileName)
    DNaseInfo = DNaseFile.readlines()
    DNaseFile.close()
    outputFile = open(outputFileName, 'w+')
    lastGenesIndex = -1
    chrom = ""
    startList = []
    start = 0
    end = 0
    lastGenesName = ""
    geneName = ""
    genesIndex = 0
    DNaseIndex = 0
    [DNaseChrom, DNaseMiddle] = getDNaseInfoSpecific(DNaseInfo, DNaseIndex)
    DNaseForGeneIndexes = []
    while genesIndex < len(genesInfo):
        # Find the DNase peaks associated with each gene
	if genesIndex > lastGenesIndex:
		# Get the information for the new gene
        	geneElements = genesInfo[genesIndex].split("\t")
		geneName = geneElements[12].strip()[1:-1]
		if geneName != lastGenesName:
			# Change the starts and end if at a new gene
			if lastGenesName != "":
				# Write the information from the last gene to the file
				writeGeneDNasesToFileUCSC(outputFile, lastGenesName, chrom, startList[0], end, DNaseForGeneIndexes)
			start = int(geneElements[4].strip())
			startList = []
			startList.append(start)
			end = int(geneElements[5].strip())
			if len(DNaseForGeneIndexes) > 0:
				# Go back to the 1st DNase for the previous gene
				DNaseIndex = DNaseForGeneIndexes[0]
				[DNaseChrom, DNaseMiddle] = getDNaseInfoSpecific(DNaseInfo, DNaseIndex)
    				DNaseForGeneIndexes = []
		else:
			start = int(geneElements[4].strip())
			startList.append(start)
			end = max(end, int(geneElements[5].strip()))
		chrom = geneElements[2].strip()[1:-1]
		lastGenesIndex = genesIndex
		lastGenesName = geneName
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
		genesIndex = genesIndex + 1
		continue	
	if DNaseChrom > chrom:
		# There are no more DNase hypersensitivity sites on the chromosome with the gene
		genesIndex = genesIndex + 1
		continue
	if DNaseMiddle > start + cutoff:
		# There are no more DNase hypersensitivity sites near the gene's TSS
		genesIndex = genesIndex + 1
		continue
	if DNaseIndex not in DNaseForGeneIndexes:
		# The DNase hypersensitivity site has not yet been considered for the gene, so add it to the list
		DNaseForGeneIndexes.append(DNaseIndex)
	DNaseIndex = DNaseIndex + 1
	if DNaseIndex >= len(DNaseInfo):
		# All DNases have been used, so stop
		genesIndex = genesIndex + 1
		continue
	[DNaseChrom, DNaseMiddle] = getDNaseInfoSpecific(DNaseInfo, DNaseIndex)
    writeGeneDNasesToFileUCSC(outputFile, lastGenesName, chrom, startList[0], end, DNaseForGeneIndexes)
    outputFile.close()

if __name__=="__main__":
    import sys
    genesFileName = sys.argv[1]
    DNaseFileName = sys.argv[2]
    outputFileName = sys.argv[3]
    cutoff = int(sys.argv[4])
                        
    associateGenesToDNaseUCSC(genesFileName, DNaseFileName, outputFileName, cutoff)
