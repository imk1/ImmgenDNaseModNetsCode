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

def getChromLengths(chromosomeLengthsFileName):
	# Get the lengths of each chromosome from the file
	chromLengths = []
	chromosomeLengthsFile = open(chromosomeLengthsFileName)
	for line in chromosomeLengthsFile:
		# Iterate through chromosomes and add their lengths to the list
		chromLengths.append(int(line.strip()))
	chromosomeLengthsFile.close()
	return chromLengths

def getNearTSSInfo(genesFileName, finalGeneListFileName, chromosomeLengthsFileName, outputFileName, cutoff):
    # Find the region within cutoff kb of each of the genes on the final list's TSS
    genesFile = open(genesFileName)
    finalGeneList = makeTFList(finalGeneListFileName)
    chromLengths = getChromLengths(chromosomeLengthsFileName)
    outputFile = open(outputFileName, 'w+')
    genesFile.readline()
    genesInfo = genesFile.readlines()
    genesFile.close()
    outputFile = open(outputFileName, 'w+')
    chrom = ""
    chromCount = 0
    lastChrom = "chr1"
    start = 0
    minusCutoff = 0
    plusCutoff = 0
    lastGenesName = ""
    geneName = ""
    genesIndex = 0
    allGenesUsed = False
    while genesIndex < len(genesInfo):
        # Find the region within cutoff kb of each gene's TSS
        geneElements = genesInfo[genesIndex].split("\t")
	geneName = string.upper(geneElements[12].strip()[1:-1])
	while geneName not in finalGeneList:
		# Iterate through genes until a gene on the final list has been found
		genesIndex = genesIndex + 1
		if genesIndex == len(genesInfo):
			# All genes have been used, so stop
			allGenesUsed = True
			break
		geneElements = genesInfo[genesIndex].split("\t")
		geneName = string.upper(geneElements[12].strip()[1:-1])
	if allGenesUsed == True:
		# All genes have been seen, so stop
		break
	if geneName != lastGenesName:
		# Change the starts and end if at a new gene
		start = int(geneElements[4].strip())
		chrom = geneElements[2].strip()[1:-1]
		if chrom != lastChrom:
			# A new chromosome has been reached
			print chrom
			chromCount = chromCount + 1
			print chromLengths[chromCount]
			lastChrom = chrom
		# Write the information from the new gene to the file
		outputFile.write(chrom)
		minusCutoffList = []
		minusCutoffList.append(1)
		minusCutoffList.append(start-cutoff)
		minusCutoff = max(minusCutoffList)
		plusCutoffList = []
		plusCutoffList.append(chromLengths[chromCount])
		plusCutoffList.append(start+cutoff)
		plusCutoff = min(plusCutoffList) 
		outputFile.write("\t")
		outputFile.write(str(minusCutoff))
		outputFile.write("\t")
		outputFile.write(str(plusCutoff))
		outputFile.write("\n")
	lastGenesName = geneName
	genesIndex = genesIndex + 1
    outputFile.close()

if __name__=="__main__":
    import sys
    import string
    genesFileName = sys.argv[1]
    finalGeneListFileName = sys.argv[2]
    chromosomeLengthsFileName = sys.argv[3]
    outputFileName = sys.argv[4]
    cutoff = int(sys.argv[5])
                        
    getNearTSSInfo(genesFileName, finalGeneListFileName, chromosomeLengthsFileName, outputFileName, cutoff)
