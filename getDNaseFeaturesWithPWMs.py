def readMat(matFileName):
	# Create a matrix of floats from a file
	floatMat = []
	matFile = open(matFileName)
	lineCount = 0
	for line in matFile:
		# Iterate through the lines of the file and make each line a row of the matrix
		floatMat.append([])
		lineElements = line.strip().split("\t")
		for numStr in lineElements:
			# Iterate through floats in the line and make each float an entry in the matrix
			floatMat[lineCount].append(float(numStr))
		lineCount = lineCount + 1
	matFile.close()
	return floatMat

def getDNaseFeaturesWithPWMs(genesToDNaseFileName, genesDNaseFeatFileName, PWMScoreMatFileName, numGenes):
	# Create a g x r matrix, where g is the number of genes and r is the number of regulators, where each entry is sum_(DNase sites)[PWM Score_(DNase site, TF) * DNase Feature_(DNase site, gene)]
	genesToDNaseFile = open(genesToDNaseFileName)
	genesDNaseFeatFile = open(genesDNaseFeatFileName)
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
		featLine = genesDNaseFeatFile.readline()
		featLineElements = featLine.strip().split("\t") # DNase feature values for each gene (with gene name at beginning)
		DNaseCount = 0
		for DNaseIndexStr in lineElements[4:len(lineElements)]:
			# Iterate through the DNase sites for the current gene and modify the features appropriately
			DNaseIndex = int(DNaseIndexStr)
			featVal = float(featLineElements[DNaseCount+1]) # DNaseIndex is 0-INDEXED, but the first entry in each row is the gene and not the DNase site
			for i in range(len(PWMScoreMat)):
				# Iterate through TFs and add the product of the PWM score and the DNase feature for each to the appropriate entry of the gene TF feature value matrix
				geneTFFeatMatPWM[geneCount][i] = geneTFFeatMatPWM[geneCount][i] + (PWMScoreMat[i][DNaseIndex] * featVal) # DNaseIndex is 0-INDEXED
			DNaseCount = DNaseCount + 1
		geneCount = geneCount + 1
	genesToDNaseFile.close()
	genesDNaseFeatFile.close()
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
    genesDNaseFeatFileName = sys.argv[2]
    PWMScoreMatFileName = sys.argv[3]
    numGenes = int(sys.argv[4])
    outputFileName = sys.argv[5]
                        
    geneTFFeatMatPWM = getDNaseFeaturesWithPWMs(genesToDNaseFileName, genesDNaseFeatFileName, PWMScoreMatFileName, numGenes)
    writeFeatGeneMat(geneTFFeatMatPWM, outputFileName)
