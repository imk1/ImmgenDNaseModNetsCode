def generateAllKMers(k):
	# Generates all k-mers of nucleotides
	bases = []
	bases.append('A')
	bases.append('C')
	bases.append('G')
	bases.append('T')
	kMers = []
	if k > 1:
		# The mers are longer than 1 base, so need to generate all k-1 mers
		merSuffixList = generateAllKMers(k - 1)
	for b in bases:
		# Iterate through bases and generate all k-1 mers for each
		if k == 1:
			# k = 1, so each k-mer contains only one base
			kMers.append(b)
		else:
			merPrefix = b
			for merSuffix in merSuffixList:
				# Append all k-1 mers to the current base
				currentMer = merPrefix + merSuffix
				kMers.append(currentMer)
	return kMers


def writeListToFile(stringList, outputFileName):
	# Writes a list of strings to a file
	outputFile = open(outputFileName, 'w+')
	for s in stringList:
		# Iterate through strings and write each string to its own line in the file
		outputFile.write(s)
		outputFile.write("\n")
	outputFile.close()


if __name__=="__main__":
    import sys
    k = int(sys.argv[1]) 
    outputFileName = sys.argv[2]  
                          
    kMers = generateAllKMers(k)
    writeListToFile(kMers, outputFileName)
