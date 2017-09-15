def getChromLengthsDict(chromosomeLengthsFileName, chromList):
	# Get the lengths of each chromosome from the file and put them into a dictionary
	chromLengths = {}
	chromosomeLengthsFile = open(chromosomeLengthsFileName)
	chromCount = 0
	for line in chromosomeLengthsFile:
		# Iterate through chromosomes and add their lengths to the list
		chromLengths[chromList[chromCount]] = int(line.strip())
		chromCount = chromCount + 1
	chromosomeLengthsFile.close()
	return chromLengths

def getRandomRegionsNotGeneSpecific(DNaseFileName, chromosomeLengthsFileName, genesToDNaseFileName, DNaseOutputFileNamePrefix, cutoff, numRand):
	# Gets random regions for each DNase hypersensitivity site, and make sure the random regions are near a gene
	DNaseOutputFileArray = []
	for i in range(numRand):
		# Open the output file for each random trial
		DNaseOutputFileName = DNaseOutputFileNamePrefix + str(i+1)
		DNaseOutputFile = open(DNaseOutputFileName, 'w+')
		DNaseOutputFileArray.append(DNaseOutputFile)
	DNaseFile = open(DNaseFileName)
	chromList = []
	chromList.append("chr1")
	chromList.append("chr10")
	chromList.append("chr11")
	chromList.append("chr12")
	chromList.append("chr13")
	chromList.append("chr14")
	chromList.append("chr15")
	chromList.append("chr16")
	chromList.append("chr17")
	chromList.append("chr18")
	chromList.append("chr19")
	chromList.append("chr2")
	chromList.append("chr3")
	chromList.append("chr4")
	chromList.append("chr5")
	chromList.append("chr6")
	chromList.append("chr7")
	chromList.append("chr8")
	chromList.append("chr9")
	chromList.append("chrX")
	chromList.append("chrY")
	chromLengths = getChromLengthsDict(chromosomeLengthsFileName, chromList)
	genesToDNaseFile = open(genesToDNaseFileName)
	genesToDNaseInfo = genesToDNaseFile.readlines()
	genesToDNaseFile.close()
	for line in DNaseFile:
		# Iterate through DNase sites and choose a random region corresponding to each site
		DNaseElements = line.split("\t")
		DNaseStart = int(DNaseElements[1])
		DNaseEnd = int(DNaseElements[2])
		DNaseLength = DNaseEnd - DNaseStart
		for i in range(numRand):
			# Find a random site corresponding to each DNase site
			# ASSUMES THAT THE REGION SURROUNDING A GENE WILL NEVER BE SHORTER THAN A DNASE SITE
			randGeneIndex = random.randint(0, len(genesToDNaseInfo)-1)
			randGeneInfo = genesToDNaseInfo[randGeneIndex]
			randGeneInfoElements = randGeneInfo.split("\t")
			chrom = randGeneInfoElements[1]
			start = int(randGeneInfoElements[2])
			minusCutoffList = []
			minusCutoffList.append(1)
			minusCutoffList.append(start-cutoff)
			minusCutoff = max(minusCutoffList)
			plusCutoffList = []
			plusCutoffList.append(chromLengths[chrom])
			plusCutoffList.append(start+cutoff-DNaseLength)
			plusCutoff = min(plusCutoffList)
			randDNaseStart = random.randint(minusCutoff, plusCutoff)
			randDNaseEnd = randDNaseStart + DNaseLength
			DNaseOutputFileArray[i].write(chrom)
			DNaseOutputFileArray[i].write("\t")
			DNaseOutputFileArray[i].write(str(randDNaseStart))
			DNaseOutputFileArray[i].write("\t")
			DNaseOutputFileArray[i].write(str(randDNaseEnd))
			DNaseOutputFileArray[i].write("\n")
	DNaseFile.close()
	for i in range(numRand):
		# Close all of the output files
		DNaseOutputFileArray[i].close()

if __name__=="__main__":
    import sys
    import random
    DNaseFileName = sys.argv[1]
    chromosomeLengthsFileName = sys.argv[2]
    genesToDNaseFileName = sys.argv[3]
    DNAseOutputFileNamePrefix = sys.argv[4]
    cutoff = int(sys.argv[5])
    numRand = int(sys.argv[6])
                        
    getRandomRegionsNotGeneSpecific(DNaseFileName, chromosomeLengthsFileName, genesToDNaseFileName, DNAseOutputFileNamePrefix, cutoff, numRand)
