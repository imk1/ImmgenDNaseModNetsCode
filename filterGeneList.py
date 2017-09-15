def filterGeneList(geneListOneFileName, geneListTwoFileName, geneListFiltFileName):
	# Find all genes in gene list one that are also in gene list two
	geneListTwoFile = open(geneListTwoFileName)
	geneListTwo = []
	for line in geneListTwoFile:
		# Iterate through lines of the file and put the gene in each line into geneListTwo
		geneListTwo.append(line.strip())
	geneListOneFile = open(geneListOneFileName)
	geneListFiltFile = open(geneListFiltFileName, 'w+')
	for line in geneListOneFile:
		# Iterate through genes in gene list one and write them to a the new file if they are also in geneListTwo
		gene = line.strip()
		if gene in geneListTwo:
			# Write the gene to the new file
			geneListFiltFile.write(gene)
			geneListFiltFile.write("\n")
	geneListOneFile.close()
	geneListFiltFile.close()

if __name__=="__main__":
    import sys
    geneListOneFileName = sys.argv[1]
    geneListTwoFileName = sys.argv[2]
    geneListFiltFileName = sys.argv[3]
    filterGeneList(geneListOneFileName, geneListTwoFileName, geneListFiltFileName)
