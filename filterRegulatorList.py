def makeTFList(TFFileName):
    # Make a list of TFs in which all of the letters are capitalized
    TFFile = open(TFFileName)
    TFList = []
    for line in TFFile:
        # Iterate through TF file and make a list of all of the TFs in upper case
        TFCap = string.upper(line.strip())
        TFList.append(TFCap)
    return TFList

def intersect(listOne, listTwo):
	# Find all of the members of listOne that are not in listTwo
	intersectList = []
	for element in listOne:
		# Iterate through the elements of listOne and add those that are not in listTwo to diffList
		if element in listTwo:
			# Element is not in listTwo, so add it to diffList
			intersectList.append(element)
	return intersectList

def filterRegulatorList(progFileName, finalTFListFileName):
	# Create a list of regulators for the module that are in the final list of regulators
	finalTFList = makeTFList(finalTFListFileName)
	progFile = open(progFileName)
	progFile.readline()
	TFListModule = []
	for line in progFile:
		# Iterate through the lines of the module file and make a list of the regulators
		lineElements = line.split(",")
		if len(lineElements[1].strip()) < 1:
			# The regulator list has ended, so stop
			break
		TFListModule.append(string.upper(lineElements[1].strip()))
	TFListModuleFiltered = intersect(TFListModule, finalTFList)
	print len(TFListModuleFiltered)
	return TFListModuleFiltered

def writeList(listForFile, fileName):
	# Writes a list to a file
	# ASSUMES THAT ALL ELEMENTS IN THE LIST ARE STRINGS
	listFile = open(fileName, 'w+')
	for element in listForFile:
		# Iterate through elements for a list and write each element to the file
		listFile.write(element)
		listFile.write("\n")
	listFile.close()

if __name__=="__main__":
    import sys
    import string
    progFileName = sys.argv[1]
    finalTFListFileName = sys.argv[2]
    outputFileName = sys.argv[3]

    TFListModuleFiltered = filterRegulatorList(progFileName, finalTFListFileName)
    writeList(TFListModuleFiltered, outputFileName)
