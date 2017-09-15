def separateGeneralCellType(dataFileName, cellTypeName, generalCellTypeFileName):
	# Create a file with the column numbers for a general cell type
	dataFile = open(dataFileName)
	cellTypeNameLine = dataFile.readline()
	cellTypeNameList = cellTypeNameLine.split("\t")
	dataFile.close()
	generalCellTypeFile = open(generalCellTypeFileName, 'w+')
	ctnCount = 0
	for ctn in cellTypeNameList:
		# Iterate through full cell type names to find the cell types with the general cell type name
		ctnElements = ctn.split(".")
		ctnGeneral = ctnElements[0]
		if ctnGeneral == cellTypeName:
			# A column with the general cell type has been found
			generalCellTypeFile.write(str(ctnCount))
			generalCellTypeFile.write("\t")
		ctnCount = ctnCount + 1
	generalCellTypeFile.close()

if __name__=="__main__":
    import sys
    dataFileName = sys.argv[1]
    cellTypeName = sys.argv[2]
    generalCellTypeFileName = sys.argv[3]
                        
    separateGeneralCellType(dataFileName, cellTypeName, generalCellTypeFileName)
