def readMat(matFileName):
	# Create a matrix of floats from a file
	floatMat = []
	matFile = open(matFileName)
	lineCount = 0
	for line in matFile:
		# Iterate through the lines of the file and make each line a row of the matrix
		if lineCount % 100 == 0:
			print lineCount
		floatMat.append([])
		lineElements = line.strip().split("\t")
		for numStr in lineElements:
			# Iterate through floats in the line and make each float an entry in the matrix
			floatMat[lineCount].append(float(numStr))
		lineCount = lineCount + 1
	matFile.close()
	return floatMat

def getDNaseFeaturesWithPWMsNull(genesToDNaseFileName, PWMScoreMatFileName, numGenes, isBinary):
	# Create a g x r matrix, where g is the number of genes and r is the number of regulators, where each entry is sum_(DNase sites)[PWM Score_(DNase site, TF)]
	genesToDNaseFile = open(genesToDNaseFileName)
	PWMScoreMat = readMat(PWMScoreMatFileName) # Matrix that is r x d, where r is the number of regulators and d is the number of DNase sites, in which each entry has the corresponding PWM score
	geneTFFeatMatPWM = []
	for i in range(numGenes):
		# Add a row to the feature matrix for each gene
		if i % 5000 == 0:
			print i
		geneTFFeatMatPWMRow = []
		for j in range(len(PWMScoreMat)):
			# Initialize geneTFFeatMatPWM so that all of the features are 0
			geneTFFeatMatPWMRow.append(0)
		geneTFFeatMatPWM.append(geneTFFeatMatPWMRow)
	geneCount = 0
	for line in genesToDNaseFile:
		# Iterate through genes and compute their feature values with the PWM scores
		if geneCount % 5000 == 0:
			print str(geneCount)
		lineElements = line.strip().split("\t") # DNase indexes for each gene (with gene name and location information at beginning)
		for DNaseIndexStr in lineElements[4:len(lineElements)]:
			# Iterate through the DNase sites for the current gene and modify the features appropriately
			DNaseIndex = int(DNaseIndexStr)
			for i in range(len(PWMScoreMat)):
				# Iterate through TFs and add the product of the PWM score and the DNase feature for each to the appropriate entry of the gene TF feature value matrix
				if isBinary == 0:
					# Not binary PWM scores
					geneTFFeatMatPWM[geneCount][i] = geneTFFeatMatPWM[geneCount][i] + PWMScoreMat[i][DNaseIndex] # DNaseIndex is 0-INDEXED
				else:
					if PWMScoreMat[i][DNaseIndex] >= .5:
						# The TF is expected to bind at least once (rounding), so add 1 to the PWM score for the gene
						geneTFFeatMatPWM[geneCount][i] = geneTFFeatMatPWM[geneCount][i] + 1
		geneCount = geneCount + 1
	genesToDNaseFile.close()
	return geneTFFeatMatPWM

def writeFeatGeneMat(featGeneMat, outputFileName):
	# Write the DNase-gene pairwise feature matrix to a file
	print "Writing to file!"
	outputFile = open(outputFileName, 'w+')
	for featRow in featGeneMat:
		# Iterate through rows of matrix and write each row to the file
		for DNaseGeneFeat in featRow:
			# Iterate through features for each gene and write each feature to the file
			outputFile.write(str(DNaseGeneFeat))
			outputFile.write("\t")
		outputFile.write("\n")
	outputFile.close()

if __name__=="__main__":
    import sys
    genesToDNaseFileName = sys.argv[1]
    PWMScoreMatFileName = sys.argv[2]
    numGenes = int(sys.argv[3])
    outputFileName = sys.argv[4]
    isBinary = int(sys.argv[5]) # Should be 0 or 1
                        
    geneTFFeatMatPWM = getDNaseFeaturesWithPWMsNull(genesToDNaseFileName, PWMScoreMatFileName, numGenes, isBinary)
    writeFeatGeneMat(geneTFFeatMatPWM, outputFileName)
