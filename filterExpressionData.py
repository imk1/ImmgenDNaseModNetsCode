def filterExpressionData(expressionInputFileName, geneListFileName):
	# Removes genes not in the final gene list from the expression data and creates a matrix with only expressions
	geneListFile = open(geneListFileName)
	geneList = []
	for line in geneListFile:
		# Iterate through genes and put them into the gene list
		geneList.append(string.upper(line.strip()))
	geneListFile.close()
	expressionInputFile = open(expressionInputFileName)
	expressionMat = []
	for i in range(len(geneList)):
		# Initialize expression matrix
		expressionMat.append([])
	for line in expressionInputFile:
		# Iterate through gene expression information, removes genes not in the list, and record only the expressions
		lineElements = line.split("\t")
		geneName = string.upper(lineElements[1].strip())
		if geneName in geneList:
			# The gene is in the gene list, so record the expressions in the output file
			geneLocation = geneList.index(geneName)
			for expr in lineElements[2:len(lineElements)]:
				# Iterate through the expressions and write each one to the output file
				if len(expr) > 1:
					# The expression is not whitespace
					expressionMat[geneLocation].append(expr)
	expressionInputFile.close()
	return expressionMat

def writeMatToFile(expressionMat, expressionOutputFileName):
	# Write a matrix to a file
	expressionOutputFile = open(expressionOutputFileName, 'w+')
	for exprRow in expressionMat:
		# Iterate through genes and write each gene's expressions to the file
		for expr in exprRow:
			# Write the expression from each cell type to the file
			expressionOutputFile.write(expr)
			expressionOutputFile.write("\t")
		expressionOutputFile.write("\n")
	expressionOutputFile.close()

if __name__=="__main__":
    import sys
    import string
    expressionInputFileName = sys.argv[1]
    geneListFileName = sys.argv[2]
    expressionOutputFileName = sys.argv[3]
                        
    expressionMat = filterExpressionData(expressionInputFileName, geneListFileName)
    writeMatToFile(expressionMat, expressionOutputFileName)
