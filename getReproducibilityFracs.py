def getCellTypeToCols(cellTypeInfoElements):
    # Makes a mapping from cell types to columns in the expression matrix
    cellTypeToColsIndex = -1
    cellTypesToCols = []
    currentCellType = ""
    for i in range(2, len(cellTypeInfoElements)):
        # Iterate through column headers and group the columns by cell types
        currentCellTypeInfo = cellTypeInfoElements[i].split("#")
        if currentCellTypeInfo[0] == currentCellType:
            # Cell type for current column is the same as cell type for previous column
            cellTypesToCols[cellTypeToColsIndex].append(i)
        else:
            cellTypeToColsIndex = cellTypeToColsIndex + 1
            cellTypesToCols.append([])
            cellTypesToCols[cellTypeToColsIndex].append(i)
            currentCellType = currentCellTypeInfo[0]
    return cellTypesToCols

def getReproducibilityFracs(expressionFileName, cutoff, repFracCutoff):
    # Gets the fraction of reproducible genes for each cell type
    expressionFile = open(expressionFileName)
    cellTypeInfo = expressionFile.readline()
    cellTypeInfoElements = cellTypeInfo.split("\t")
    cellTypesToCols = getCellTypeToCols(cellTypeInfoElements)
    numGenes = 0
    repCounts = []
    nonRepGenes = []
    for i in range(len(cellTypesToCols)):
        # Initialize reproducibility counts
        repCounts.append(0)
    for line in expressionFile:
        # Iterate through expressions and update counts
        lineElements = line.split("\t")
	gene = lineElements[1]
	numNonRep = 0
        for ct in range(len(cellTypesToCols)):
            # Iterate through cell types and find out whether the gene is reproducible for each
            expressions = []
            for col in cellTypesToCols[ct]:
                # Iterate through expressions and convert them into numbers
                expr = float(lineElements[col])
                expressions.append(expr)
            maxExpr = max(expressions)
            minExpr = min(expressions)
            if float(maxExpr)-float(minExpr) > cutoff:
                # Reproducibility problem
                repCounts[ct] = repCounts[ct] + 1
		if ct < len(cellTypesToCols) - 2:
			# Not a control
			numNonRep = numNonRep + 1
	if float(numNonRep)/float(len(cellTypesToCols) - 2) > repFracCutoff:
		# Gene is not reproducible for too high of a fraction of cell types
		nonRepGenes.append(gene)
        numGenes = numGenes + 1
    repFracs = []
    for rc in repCounts:
        # Use the counts to find the fraction of genes for which each cell type's expressions are not reproducible
        print rc
        frac = float(rc)/float(numGenes)
        repFracs.append(frac)
    print float(len(nonRepGenes))/float(numGenes)
    expressionFile.close()
    return [cellTypesToCols, repFracs, nonRepGenes]

def writeRepFracs(repFracs, outputFileName):
    # Write reproducibility fractions to a file
    outputFile = open(outputFileName, 'w+')
    for rf in repFracs:
        # Iterate through reproduciblity fractions and write them to the file
        outputFile.write(str(rf))
        outputFile.write("\n")

def writeCellTypesToCols(cellTypesToCols, cellTypesColsFileName):
    # Write the columns for each cell type to a file
    cellTypesColsFile = open(cellTypesColsFileName, 'w+')
    for cellTypeCols in cellTypesToCols:
        # Iterate through cell types and write the columns for each to the file
        for col in cellTypeCols:
            # Write each column number for the current cell type to the file
            cellTypesColsFile.write(str(col))
            cellTypesColsFile.write("\t")
        cellTypesColsFile.write("\n")
    cellTypesColsFile.close()

def writeList(geneList, geneListFileName):
	# Write a list of genes to a file
	geneListFile = open(geneListFileName, 'w+')
	for gene in geneList:
		# Write each gene to the file
		geneListFile.write(gene)
		geneListFile.write("\n")
	geneListFile.close()

if __name__=="__main__":
    import sys
    expressionFileName = sys.argv[1]
    cutoff = float(sys.argv[2])
    repFracCutoff = float(sys.argv[3])
    outputFileName = sys.argv[4]
    cellTypesColsFileName = sys.argv[5]
    nonRepGenesFileName = sys.argv[6]                    
    [cellTypesToCols, repFracs, nonRepGenes] = getReproducibilityFracs(expressionFileName, cutoff, repFracCutoff)
    writeRepFracs(repFracs, outputFileName)
    writeCellTypesToCols(cellTypesToCols, cellTypesColsFileName)
    writeList(nonRepGenes, nonRepGenesFileName)
