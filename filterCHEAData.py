def filterCHEAData(CHEADataFileName, TFFileName, geneFileName, acceptedCellTypeFileName, CHEADataOutFileName, CHEALocationOutFileName):
	# Filters data from the CHEA site
	CHEADataFile = open(CHEADataFileName)
	TFFile = open(TFFileName)
	geneFile = open(geneFileName)
	acceptedCellTypeFile = open(acceptedCellTypeFileName)
	TFList = []
	for line in TFFile:
		# Add each TF to the list of TFs
		TFList.append(string.upper(line.strip()))
	geneList = []
	for line in geneFile:
		# Add each TF to the list of TFs
		geneList.append(string.upper(line.strip()))
	acceptedCellTypeList = []
	for line in acceptedCellTypeFile:
		# Add each accepted cell type to the list of accepted cell types
		acceptedCellTypeList.append(line.strip())
	TFFile.close()
	geneFile.close()
	acceptedCellTypeFile.close()
	CHEADataOutFile = open(CHEADataOutFileName, 'w+')
	CHEALocationOutFile = open(CHEALocationOutFileName, 'w+')
	for line in CHEADataFile:
		# Iterate through CHEA data and store the interactions from mouse and selected TFs
		lineElements = line.split(",")
		species = lineElements[7]
		#if species != "mouse":
		#	# Iteraction is not from mouse, so do not consider it
		#	continue
		cellType = lineElements[6]
		if cellType not in acceptedCellTypeList:
			# The cell type with the interaction is not a cell type with expression data, so do not consider interaction
			continue
		TF = string.upper(lineElements[1])
		TFLocation = -1
		if TF not in TFList:
			# TF was not selected, so do not consider its interaction
			continue
		else:
			TFLocation = TFList.index(TF) # TF locations are 0-INDEXED
		gene = string.upper(lineElements[3])
		geneLocation = -1
		if gene not in geneList:
			# Gene was not selected, so do not consider its interaction
			continue
		else:
			geneLocation = geneList.index(gene) # Gene locations are 0-INDEXED
		CHEADataOutFile.write(TF)
		CHEADataOutFile.write("\t")
		CHEADataOutFile.write(gene)
		CHEADataOutFile.write("\n")
		CHEALocationOutFile.write(str(TFLocation))
		CHEALocationOutFile.write("\t")
		CHEALocationOutFile.write(str(geneLocation))
		CHEALocationOutFile.write("\n")
	CHEADataFile.close()
	CHEADataOutFile.close()
	CHEALocationOutFile.close()

if __name__=="__main__":
    import sys
    import string
    CHEADataFileName = sys.argv[1]
    TFFileName = sys.argv[2]
    geneFileName = sys.argv[3]
    acceptedCellTypeFileName = sys.argv[4]
    CHEADataOutFileName = sys.argv[5]
    CHEALocationOutFileName = sys.argv[6]
    filterCHEAData(CHEADataFileName, TFFileName, geneFileName, acceptedCellTypeFileName, CHEADataOutFileName, CHEALocationOutFileName)
