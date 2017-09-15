def filterCellTypeNames(cellTypeNamesFileName, outputFileName):
	# Filters the cell type names so that there is only 1 cell type name for each cell type
	cellTypes = []
	cellTypeNamesFile = open(cellTypeNamesFileName)
	outputFile = open(outputFileName, 'w+')
	for line in cellTypeNamesFile:
		# Iterate through cell types and write the unique ones to the output file
		lineElements = line.split("#")
		currentCellType = lineElements[0]
		if currentCellType not in cellTypes:
			# The current cell type has not yet been recorded, so record it and add it to the list of cell types
			outputFile.write(currentCellType)
			outputFile.write("\n")
			cellTypes.append(currentCellType)
	cellTypeNamesFile.close()
	outputFile.close()

if __name__=="__main__":
    import sys
    cellTypeNamesFileName = sys.argv[1] 
    outputFileName = sys.argv[2]                            
    filterCellTypeNames(cellTypeNamesFileName, outputFileName)
