def makeGeneToTSS(geneInfo):
	# Makes a dictionary that maps each gene to its TSS
	# If there are multiple transcripts, ASSUMES THAT THE TRANSCRIPT WITH THE EARLIEST TSS IS LISTED FIRST
	geneToTSS = {}
	for currentGeneInfo in geneInfo:
		# Iterate through genes and enter the gene and its TSS into the dictionary
		geneElements = currentGeneInfo.split("\t")
		geneName = geneElements[12].strip()[1:-1]
		if geneName in geneToTSS:
			# The TSS for the gene has already been recorded because there are multiple transcripts
			continue
		start = int(geneElements[4].strip())
		geneToTSS[string.upper(geneName)] = start
	return geneToTSS

def getRandomRegions(genesFileName, DNaseFileName, genesToDNaseFileName, outputFileNamePrefix, cutoff, numRand):
	# Gets random regions for each DNase hypersensitivity site that is associated with each gene
	# ASSUMES THAT TRANSCRIPTS FROM DIFFERENT MULTI-TRANSCRIPT GENES DO NOT OVERLAP, THAT ALL TRANSCRIPTS IN A GENE ARE CONSECUTIVE, AND THAT THERE ARE NO GENES WITH THE SAME NAME ON DIFFERENT CHROMOSOMES
    	# DNase hypersensitivity sites are 0-INDEXED
	outputFileArray = []
	for i in range(numRand):
		# Open the output file for each random trial
		outputFileName = outputFileNamePrefix + str(i+1)
		outputFile = open(outputFileName, 'w+')
		outputFileArray.append(outputFile)
    	genesFile = open(genesFileName)
    	genesFile.readline()
    	genesInfo = genesFile.readlines()
    	genesFile.close()
	geneToTSS = makeGeneToTSS(genesInfo)
   	DNaseFile = open(DNaseFileName)
    	DNaseInfo = DNaseFile.readlines()
    	DNaseFile.close()
	genesToDNaseFile = open(genesToDNaseFileName)
	for line in genesToDNaseFile:
		# Iterate through genes and their associated DNase sites and choose random regions for each site
		# IF 2 GENES ARE ASSOCIATED WITH THE SAME SITE, THE RANDOM SITES CHOSEN FOR EACH GENE MIGHT BE DIFFERENT
		lineElements = line.split("\t")
		geneName = string.upper(lineElements[0])

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
