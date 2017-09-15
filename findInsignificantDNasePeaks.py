def findInsignificantDNasePeaks(DNaseFileName, cutoff, DNasesToRemoveFileName):
	# Find DNase peaks whose signals are not significant
	# INDEXES OF PEAKS ARE 1-INDEXED (NOT 0-INDEXED)
	DNaseFile = open(DNaseFileName)
	DNasesToRemoveFile = open(DNasesToRemoveFileName, 'w+')
	lineCount = 0
	for line in DNaseFile:
		# Iterate through the lines of the DNase file and make a list of the lines whose signals are not significant
		lineCount = lineCount + 1
		lineElements = line.split("\t")
		pVal = float(lineElements[7])
		if pVal <= cutoff:
			# The p-value is not sufficiently low, so designiate the peak for removal
			# Doing <= instead of > because p-values are -log10ed
			DNasesToRemoveFile.write(str(lineCount))
			DNasesToRemoveFile.write("\n")
	DNaseFile.close()
	DNasesToRemoveFile.close()

if __name__=="__main__":
    import sys
    DNaseFileName = sys.argv[1]
    cutoff = float(sys.argv[2])
    DNasesToRemoveFileName = sys.argv[3]
                        
    findInsignificantDNasePeaks(DNaseFileName, cutoff, DNasesToRemoveFileName)
