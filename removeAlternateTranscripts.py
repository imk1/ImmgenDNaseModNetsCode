def findAlternateTranscripts(expressionFileName):
	# Make an array with 4 columns: gene name, line number, average expression, isMax and a list of lines to remove
	exprInfoArray = []
	linesToRemove = []
	for i in range(4):
		# Initialize the array with the expression information
		exprInfoArray.append([])
	expressionFile = open(expressionFileName)
	expressionFile.readline()
	count = 0
	for line in expressionFile:
		# Iterate through genes, record necessary expression information, and put repeats' lines on the removal list
		if math.fmod(count, 1000) == 0:
			print count
		lineElements = line.split("\t")
		geneName = lineElements[1]
		expressions = []
		for exprString in lineElements[2:len(lineElements)-1]:
			# Iterate through expressions and put them into an array
			expressions.append(float(exprString.strip()))
		meanExpr = float(sum(expressions))/float(len(expressions))
		isMax = True
		if geneName in exprInfoArray[0]:
			# The gene has an alternate transcript
			maxFound = False
			geneLoc = 0
			while maxFound == False:
				# Iterate through the alternative transcripts to find the previous maximum expression
				geneLoc = geneLoc + exprInfoArray[0][geneLoc:len(exprInfoArray[0])].index(geneName)
				if exprInfoArray[3][geneLoc] == True:
					# The previous transcript with the maximum expression has been found
					maxFound = True
					break
				geneLoc = geneLoc + 1
			if exprInfoArray[2][geneLoc] >= meanExpr:
				# The earlier transcript will be used
				isMax = False
				linesToRemove.append(count)
			else:
				exprInfoArray[3][geneLoc] = False
				linesToRemove.append(exprInfoArray[1][geneLoc])
		exprInfoArray[0].append(geneName)
		exprInfoArray[1].append(count)
		exprInfoArray[2].append(meanExpr)
		exprInfoArray[3].append(isMax)
		count = count + 1
	expressionFile.close()
	return linesToRemove

def removeAlternateTranscripts(expressionFileName, linesToRemove, outputFileName):
	# Iterate through gene expression data and remove alternate transcripts
	expressionFile = open(expressionFileName)
	firstLine = expressionFile.readline()
	outputFile = open(outputFileName, 'w+')
	outputFile.write(firstLine)
	count = 0
	for line in expressionFile:
		# Iterate through lines with genes and record every line except for those in the removal list
		if count not in linesToRemove:
			# The current gene should not be removed, so record it
			outputFile.write(line)
		count = count + 1
	expressionFile.close()
	outputFile.close()

if __name__=="__main__":
    import sys
    import math
    expressionFileName = sys.argv[1]
    outputFileName = sys.argv[2]

    linesToRemove = findAlternateTranscripts(expressionFileName)                 
    removeAlternateTranscripts(expressionFileName, linesToRemove, outputFileName)
