def writeOut(lineElements, lastChrom, lastStart, lastEnd, lastScore, lastBestpVal, outputFile):
	# Write a line of DNase hypersensitivity information to the file
	outputFile.write(lastChrom)
	outputFile.write("\t")
	outputFile.write(str(lastStart))
	outputFile.write("\t")
	outputFile.write(str(lastEnd))
	outputFile.write("\t")
	for i in range(3, 6):
		# Write the next 3 columns to the file
		outputFile.write(lineElements[i])
		outputFile.write("\t")
	outputFile.write(str(lastScore))
	outputFile.write("\t")
	outputFile.write(str(lastBestpVal))
	outputFile.write("\t")
	outputFile.write(lineElements[8])
	outputFile.write("\t")
	outputFile.write(lineElements[9])

def combineDNaseRegions(DNaseFileName, outputFileName):
	# Iterate through DNase regions, combine overlapping regions, and write everything to a file
	DNaseFile = open(DNaseFileName)
	outputFile = open(outputFileName, 'w+')
	lastChrom = ""
	lastStart = 0
	lastEnd = 0
	lastScore = 0
	lastBestpVal = 0
	for line in DNaseFile:
		# Iterate through DNase regions and combine the overlapping regions
		# Change the start to the start of the 1st overlapping region and the end to the end of the last
		# Change the score to the score with the lowest (highest because in negative log10 space) p-value
		# Change the p-value to the lowest (highest because in negative log10 space) p-value
		lineElements = line.split("\t")
		chrom = lineElements[0]
		start = int(lineElements[1])
		end = int(lineElements[2])
		score = float(lineElements[6])
		pVal = float(lineElements[7])
		if chrom != lastChrom:
			# The chromosome has changed, so write the last line to the file
			if lastChrom != "":
				# Not at the first line
				writeOut(lineElements, lastChrom, lastStart, lastEnd, lastScore, lastBestpVal, outputFile)
			lastChrom = chrom
			lastStart = start
			lastEnd = end
			lastScore = score
			lastBestpVal = pVal
			continue
		if start >= lastEnd:
			# The current region does not overlap with the last region, so write the last line to the file
			writeOut(lineElements, lastChrom, lastStart, lastEnd, lastScore, lastBestpVal, outputFile)
			lastChrom = chrom
			lastStart = start
			lastEnd = end
			lastScore = score
			lastBestpVal = pVal
			continue
		lastEnd = end
		if pVal > lastBestpVal:
			# The current peak is more significant, so replace the last peak's score with the current peak's
			lastBestpVal = pVal
			lastScore = score
			continue
		if pVal == lastBestpVal:
			# Average the last and current score because they are equally significant
			# ASSUMES THAT NO MORE THAN 2 PEAKS PER REGION WILL HAVE THE SAME P-VALUE
			lastScore = float(lastScore + score)/float(2)
	writeOut(lineElements, lastChrom, lastStart, lastEnd, lastScore, lastBestpVal, outputFile)
	DNaseFile.close()
	outputFile.close()

if __name__=="__main__":
    import sys
    DNaseFileName = sys.argv[1] 
    outputFileName = sys.argv[2]                            
    combineDNaseRegions(DNaseFileName, outputFileName)
