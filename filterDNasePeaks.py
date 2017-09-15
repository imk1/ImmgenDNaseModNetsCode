def filterDNasePeaks(DNasePeaksFileName, DNasesToRemoveFileName, DNasePeaksFilteredFileName):
	# Remove the list of DNase peaks for removal
	DNasesToRemoveFile = open(DNasesToRemoveFileName)
	DNasesToRemove = []
	for line in DNasesToRemoveFile:
		# Iterate through the DNases to remove and put them into an array
		index = int(line.strip())
		DNasesToRemove.append(index)
	DNasesToRemoveFile.close()
	DNasePeaksFile = open(DNasePeaksFileName)
	DNasePeaksFilteredFile = open(DNasePeaksFilteredFileName, 'w+')
	count = 0
	for line in DNasePeaksFile:
		# Record the line if its line number is not on the removal list
		if count % 10000 == 0:
			print count
		count = count + 1
		if count not in DNasesToRemove:
			# Write the current peak to the filtered file
			DNasePeaksFilteredFile.write(line)
	DNasePeaksFile.close()
	DNasePeaksFilteredFile.close()

if __name__=="__main__":
    import sys
    DNasePeaksFileName = sys.argv[1]
    DNasesToRemoveFileName = sys.argv[2]
    DNasePeaksFilteredFileName = sys.argv[3]
                        
    filterDNasePeaks(DNasePeaksFileName, DNasesToRemoveFileName, DNasePeaksFilteredFileName)
