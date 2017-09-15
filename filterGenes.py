def makeTFList(TFFileName):
    # Make a list of TFs in which all of the letters are capitalized
    TFFile = open(TFFileName)
    TFList = []
    for line in TFFile:
        # Iterate through TF file and make a list of all of the TFs in upper case
        TFCap = string.upper(line.strip())
        TFList.append(TFCap)
    return TFList

def filterGenes(expressionFileName, nonRepGenesFileName, multiLocGenesFileName, stdCutoff, TFListFileName, expressionFilteredFileName, TFExpressionFilteredFileName):
	# Removes non-reproducible genes and genes whose expression does not change from the list of genes and creates 2 files of remaining genes, 1 for everything and 1 for TFs
	expressionFile = open(expressionFileName)
    	nonRepGenes = makeTFList(nonRepGenesFileName)
	multiLocGenes = makeTFList(multiLocGenesFileName)
	TFList = makeTFList(TFListFileName)
	nonChangingGenes = []
	expressionFilteredFile = open(expressionFilteredFileName, 'w+')
	TFExpressionFilteredFile = open(TFExpressionFilteredFileName, 'w+')
	numGenes = 0
    	for line in expressionFile:
        	# Iterate through expressions and update counts
        	lineElements = line.split("\t")
		gene = string.upper(lineElements[1])
		if (gene in nonRepGenes) or (gene in multiLocGenes):
			# Gene is not reproducible or has multiple loci, so do not include it
			continue
		expressions = []
		#for exprString in lineElements[2:744]:
		for exprString in lineElements[2:len(lineElements)-1]: # Remove new line character
			# Make an array of expressions
			expressions.append(float(exprString))
		expressionsArray = array(expressions)
		if expressionsArray.std() <= stdCutoff:
			# Expression does not vary sufficiently
			nonChangingGenes.append(gene)
			continue
		if gene in TFList:
			# The gene is a TF
			TFExpressionFilteredFile.write(line)
		expressionFilteredFile.write(line)
		numGenes = numGenes + 1
	print numGenes
	return nonChangingGenes

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
    import string
    import math
    from numpy import *
    expressionFileName = sys.argv[1]
    nonRepGenesFileName = sys.argv[2]
    multiLocGenesFileName = sys.argv[3]
    stdCutoff = float(sys.argv[4])
    TFListFileName = sys.argv[5]
    expressionFilteredFileName = sys.argv[6]
    TFExpressionFilteredFileName = sys.argv[7]
    nonChangingGenesFileName = sys.argv[8]
                       
    nonChangingGenes = filterGenes(expressionFileName, nonRepGenesFileName, multiLocGenesFileName, stdCutoff, TFListFileName, expressionFilteredFileName, TFExpressionFilteredFileName)
    writeList(nonChangingGenes, nonChangingGenesFileName)
