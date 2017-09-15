def createGoDataMatrix(GOAnnotationsFileName):
	# Create a list of all genes with GO annotations, a list of all GO terms, and a gene x GO term matrix with a 1 for each (gene, GO term) pair and a 0 in all other locations
	GOAnnotationsFile = open(GOAnnotationsFileName)
	geneList = []
	GOList = []
	geneGO = []
	geneGOLocation = -1
	for i in range(31):
		# Iterate through the lines of the file with no GO annotations
		GOAnnotationsFile.readline()
	for line in GOAnnotationsFile:
		# Read GO annotation and, if necessary, add to the appropriate lists
		lineElements = line.split("\t")
		gene = lineElements[2]
		if gene not in geneList:
			# Put in the information for the new gene
			geneList.append(gene)
			geneGOLocation = geneGOLocation + 1
			geneGO.append([])
			for i in range(len(GOList)):
				# Initialize all GO entries to 0
				geneGO[geneGOLocation].append(0)
		geneLocation = geneList.index(gene)
		GO = lineElements[4]
		if GO not in GOList:
			# Put in the information for the new GO term
			GOList.append(GO)
			for geneGOInfo in geneGO:
				# Iterate through genes and add a 0 to each gene's array for the new GO term
				geneGOInfo.append(0)
		GOLocation = GOList.index(GO)
		geneGO[geneLocation][GOLocation] = 1
	GOAnnotationsFile.close()
	return [geneList, GOList, geneGO]

def writeList(nameList, listFileName):
	# Write a list of strings to a file
	listFile = open(listFileName, 'w+')
	for name in nameList:
		# Write each string to the file
		listFile.write(name)
		listFile.write("\n")
	listFile.close()

def writeTable(table, tableFileName):
	# Write a table of numbers to a file
	tableFile = open(tableFileName, 'w+')
	for row in table:
		# Iterate through the rows of the table
		for entry in row:
			# Iterate through the numbers in the row
			tableFile.write(str(entry))
			tableFile.write("\t")
		tableFile.write("\n")
	tableFile.close()

if __name__=="__main__":
    import sys
    GOAnnotationsFileName = sys.argv[1]
    geneFileName = sys.argv[2]
    GOFileName = sys.argv[3]
    geneGOFileName = sys.argv[4]
                        
    [geneList, GOList, geneGO] = createGoDataMatrix(GOAnnotationsFileName)
    writeList(geneList, geneFileName)
    writeList(GOList, GOFileName)
    writeTable(geneGO, geneGOFileName)
