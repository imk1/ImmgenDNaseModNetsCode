def getChromLengths(chromosomeLengthsFileName):
	# Get the lengths of each chromosome from the file
	chromLengths = []
	chromosomeLengthsFile = open(chromosomeLengthsFileName)
	for line in chromosomeLengthsFile:
		# Iterate through chromosomes and add their lengths to the list
		chromLengths.append(int(line.strip()))
	chromosomeLengthsFile.close()
	return chromLengths

def getDNaseInfoSpecificPlus(DNaseInfo, DNaseIndex):
    # Gets the information about a specific DNase hypersensitivity site
    DNaseElements = DNaseInfo[DNaseIndex].split("\t")
    DNaseChrom = DNaseElements[0]
    DNaseStart = int(DNaseElements[1])
    DNaseEnd = int(DNaseElements[2])
    DNaseLength = DNaseEnd - DNaseStart
    return [DNaseChrom, DNaseLength]

def getRandomRegions(DNaseFileName, chromosomeLengthsFileName, genesToDNaseFileName, DNaseOutputFileNamePrefix, mappingOutputFileNamePrefix, cutoff, numRand):
	# Gets random regions for each DNase hypersensitivity site that is associated with each gene
    	# DNase hypersensitivity sites are 0-INDEXED
	DNaseOutputFileArray = []
	mappingOutputFileArray = []
	for i in range(numRand):
		# Open the output file for each random trial
		DNaseOutputFileName = DNaseOutputFileNamePrefix + str(i+1)
		DNaseOutputFile = open(DNaseOutputFileName, 'w+')
		DNaseOutputFileArray.append(DNaseOutputFile)
		mappingOutputFileName = mappingOutputFileNamePrefix + str(i+1)
		mappingOutputFile = open(mappingOutputFileName, 'w+')
		mappingOutputFileArray.append(mappingOutputFile)
   	DNaseFile = open(DNaseFileName)
    	DNaseInfo = DNaseFile.readlines()
    	DNaseFile.close()
	chromLengths = getChromLengths(chromosomeLengthsFileName)
	chrom = ""
    	chromCount = 0
    	lastChrom = "chr1"
    	start = 0
    	minusCutoff = 0
    	plusCutoff = 0
	genesToDNaseFile = open(genesToDNaseFileName)
	randDNaseCount = 0
	for line in genesToDNaseFile:
		# Iterate through genes and their associated DNase sites and choose random regions for each site
		# IF 2 GENES ARE ASSOCIATED WITH THE SAME SITE, THE RANDOM SITES CHOSEN FOR EACH GENE MIGHT BE DIFFERENT
		lineElements = line.strip().split("\t")
		start = int(lineElements[2].strip())
		chrom = lineElements[1].strip()
		if chrom != lastChrom:
			# A new chromosome has been reached
			print chrom
			chromCount = chromCount + 1
			print chromLengths[chromCount]
			lastChrom = chrom
		minusCutoffList = []
		minusCutoffList.append(1)
		minusCutoffList.append(start-cutoff)
		minusCutoff = max(minusCutoffList)
		plusCutoffList = []
		plusCutoffList.append(chromLengths[chromCount])
		plusCutoffList.append(start+cutoff)
		plusCutoff = min(plusCutoffList)
		if len(lineElements) <= 4:
			# The current gene has no near-by DNase sites
			for i in range(numRand):
				# Write the information of the gene to each file with the random regions
				mappingOutputFileArray[i].write(line)
			continue
		for i in range(numRand):
			# Write the gene name and location information for each gene to each random file
			mappingOutputFileArray[i].write(lineElements[0])
			mappingOutputFileArray[i].write("\t")
			mappingOutputFileArray[i].write(lineElements[1])
			mappingOutputFileArray[i].write("\t")
			mappingOutputFileArray[i].write(lineElements[2])
			mappingOutputFileArray[i].write("\t")
			mappingOutputFileArray[i].write(lineElements[3])
		for DNaseSite in lineElements[4:]:
			# Iterate through DNase sites and find the corresponding random DNase site for each
			DNaseSiteNum = int(DNaseSite)
			[DNaseChrom, DNaseLength] = getDNaseInfoSpecificPlus(DNaseInfo, DNaseSiteNum)
			for i in range(numRand):
				# Choose different random region for each random iteration
				randDNaseMiddle = random.randint(minusCutoff, plusCutoff)
				# ASSUMES THAT NO GENES ARE SO CLOSE TO CHROMOSOME ENDS THAT RANDOM SITES GO OFF THE ENDS
				randDNaseStart = randDNaseMiddle - (DNaseLength/2)
				randDNaseEnd = randDNaseMiddle + (DNaseLength/2)
				DNaseOutputFileArray[i].write(DNaseChrom)
				DNaseOutputFileArray[i].write("\t")
				DNaseOutputFileArray[i].write(str(randDNaseStart))
				DNaseOutputFileArray[i].write("\t")
				DNaseOutputFileArray[i].write(str(randDNaseEnd))
				DNaseOutputFileArray[i].write("\n")
				mappingOutputFileArray[i].write("\t")
				mappingOutputFileArray[i].write(str(randDNaseCount))
			randDNaseCount = randDNaseCount + 1
		for i in range(numRand):
			# Write a line break to each mapping file
			mappingOutputFileArray[i].write("\n")
	genesToDNaseFile.close()
	for i in range(numRand):
		# Close all of the output files
		DNaseOutputFileArray[i].close()
		mappingOutputFileArray[i].close()

if __name__=="__main__":
    import sys
    import random
    DNaseFileName = sys.argv[1]
    chromosomeLengthsFileName = sys.argv[2]
    genesToDNaseFileName = sys.argv[3]
    DNAseOutputFileNamePrefix = sys.argv[4]
    mappingOutputFileNamePrefix = sys.argv[5]
    cutoff = int(sys.argv[6])
    numRand = int(sys.argv[7])
                        
    getRandomRegions(DNaseFileName, chromosomeLengthsFileName, genesToDNaseFileName, DNAseOutputFileNamePrefix, mappingOutputFileNamePrefix, cutoff, numRand)
