def logGeneExpression(inputFileName, outputFileName):
	# logs the expression values for all genes in the input file and records the new expression in the output file
	# Also removes the 1st 2 lines of the input file, but keeps the 3rd line, which has the cell types
	inputFile = open(inputFileName)
	outputFile = open(outputFileName, 'w+')
	inputFile.readline()
	inputFile.readline()
	cellTypes = inputFile.readline()
	outputFile.write(cellTypes)
	for line in inputFile:
		# Iterates through gene expressions and takes the log of each for each cell type
		lineElements = line.split("\t")
		outputFile.write(lineElements[0])
		outputFile.write("\t")
		outputFile.write(lineElements[1])
		outputFile.write("\t")
		for exprString in lineElements[2:len(lineElements)]:
			# Iterate through the gene expressions for the different cell types and take the log2 of each
			expr = float(exprString.strip())
			logExpr = math.log(expr, 2)
			outputFile.write(str(logExpr))
			outputFile.write("\t")
		outputFile.write("\n")
	inputFile.close()
	outputFile.close()

if __name__=="__main__":
    import sys
    import math
    inputFileName = sys.argv[1]
    outputFileName = sys.argv[2]                 
    logGeneExpression(inputFileName, outputFileName)
