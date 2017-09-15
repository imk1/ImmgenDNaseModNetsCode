def makeTFList(TFFileName):
    # Make a list of TFs in which all of the letters are capitalized
    TFFile = open(TFFileName)
    TFList = []
    for line in TFFile:
        # Iterate through TF file and make a list of all of the TFs in upper case
        TFCap = string.upper(line.strip())
        TFList.append(TFCap)
    return TFList

def makeInitialModules(moduleFileNameListFile, finalGeneListFileName):
	# Make the initial list of modules
	moduleFileNameList = open(moduleFileNameListFile)
	geneListFinal = makeTFList(finalGeneListFileName)
	assign = []
	for i in range(len(geneListFinal)):
		# Initialize the array of module assignments
		assign.append(0)
	for moduleFileName in moduleFileNameList:
		# Iterate through modules, find the location of each gene in each module, and set its assign value
		moduleFileNameElements = moduleFileName.split("/")
		moduleNumElements = moduleFileNameElements[5].split(".")
		moduleNum = int(moduleNumElements[0][4:])
		print moduleNum
		moduleFile = open(moduleFileName.strip())
		moduleFile.readline()
		for line in moduleFile:
			# Iterate through the genes in the module and record their modules in the assign array
			if line.strip() == ",":
				# All of the genes have been read, so stop
				break
			lineElements = line.split(",")
			geneName = lineElements[0]
			if geneName in geneListFinal:
				# The gene has not been filtered out, so record its module
				geneIndex = geneListFinal.index(geneName)
				assign[geneIndex] = moduleNum
		moduleFile.close()
	moduleFileNameList.close()
	return assign

def writeNumsToFile(numList, numListFileName):
	# Write a list of numbers to a file
	numListFile = open(numListFileName, 'w+')
	for num in numList:
		# Iterate through the numbers and write each number to its own line in the file
		numListFile.write(str(num))
		numListFile.write("\n")
	numListFile.close()

if __name__=="__main__":
    import sys
    import string
    moduleFileNameListFile = sys.argv[1] 
    finalGeneListFileName = sys.argv[2]
    numListFileName = sys.argv[3]
                            
    assign = makeInitialModules(moduleFileNameListFile, finalGeneListFileName)
    writeNumsToFile(assign, numListFileName)
