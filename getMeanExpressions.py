def makeCellTypesColsList(cellTypeColsFileName):
	# Make a list of columns for each cell type
	cellTypeColsFile = open(cellTypeColsFileName)
	cellTypesToCols = []
	lineCount = 0
	for line in cellTypeColsFile:
		# Iterate through cell types and make a list of the columns for each
		cellTypesToCols.append([])
		lineElements = line.split()
		for colNum in lineElements:
			# Put each column number into the array for this cell type
			cellTypesToCols[lineCount].append(int(colNum))
		lineCount = lineCount + 1
	cellTypeColsFile.close()
	return cellTypesToCols

def getAverage(floatList):
	# Fine the average of a list of floats
	listSum = sum(floatList)
	average = listSum/len(floatList)
	return average

def getMeanExpressions(expressionFileName, cellTypesColsFileName, expressionsMeansFileName):
	# Gets the mean value of the expressions for each cell type
	expressionFile = open(expressionFileName)
	expressionFile.readline()
	cellTypesToCols = makeCellTypesColsList(cellTypesColsFileName)
	expressionMeansFile = open(expressionsMeansFileName, 'w+')
	for line in expressionFile:
		# Iterate through expressions and average across cell types
	        lineElements = line.split("\t")
		expressionMeansFile.write(lineElements[0])
		expressionMeansFile.write("\t")
		expressionMeansFile.write(lineElements[1])
		expressionMeansFile.write("\t")
	        for ct in range(len(cellTypesToCols) - 2):
	            # Iterate through cell types and find the mean expression for each, EXCLUDE CONTROLS
	            expressions = []
	            for col in cellTypesToCols[ct]:
	                # Iterate through expressions and convert them into numbers
	                expr = float(lineElements[col])
	                expressions.append(expr)
		    meanExpressions = getAverage(expressions)
		    expressionMeansFile.write(str(meanExpressions))
		    expressionMeansFile.write("\t")
		expressionMeansFile.write("\n")
	expressionFile.close()
	expressionMeansFile.close()

if __name__=="__main__":
    import sys
    expressionFileName = sys.argv[1]
    cellTypesColsFileName = sys.argv[2]
    expressionsMeansFileName = sys.argv[3]
                        
    getMeanExpressions(expressionFileName, cellTypesColsFileName, expressionsMeansFileName)
