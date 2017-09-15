def getMaxPWMScores(maxPWMScoreFileName, meanPWMScoreFileName, numDNaseSites, inputPWMScoreFileNames):
	# Gets the maximum and mean PWM scores for TFs for which there are mutliple motifs
	maxPWMScoreFile = open(maxPWMScoreFileName, 'w+')
	meanPWMScoreFile = open(meanPWMScoreFileName, 'w+')
	inputPWMScoreFiles = []
	for ipsfn in inputPWMScoreFileNames:
		# Open each file with PWM scores
		ipsf = open(ipsfn)
		inputPWMScoreFiles.append(ipsf)
	for i in range(0, numDNaseSites):
		# Iterate through DNase sites and find the maximum and minimum PWM scores for each
		PWMScoreList = []
		for j in range(0, len(inputPWMScoreFiles)):
			# Iterate through PWM score files and find the scores for the current DNase site
			currentPWMLine = inputPWMScoreFiles[j].readline()
			currentPWMLineElements = currentPWMLine.split("\t")
			if j == 0:
				# Write the DNase location information to the output files
				maxPWMScoreFile.write(currentPWMLineElements[0])
				maxPWMScoreFile.write("\t")
				meanPWMScoreFile.write(currentPWMLineElements[0])
				meanPWMScoreFile.write("\t")
			PWMScoreList.append(float(currentPWMLineElements[1]))
		maxPWMScore = max(PWMScoreList)
		maxPWMScoreFile.write(str(maxPWMScore))
		maxPWMScoreFile.write("\n")
		meanPWMScore = float(sum(PWMScoreList))/float(len(PWMScoreList))
		meanPWMScoreFile.write(str(meanPWMScore))
		meanPWMScoreFile.write("\n")
	maxPWMScoreFile.close()
	meanPWMScoreFile.close()
	for ipsf in inputPWMScoreFiles:
		# Close each of the input files
		ipsf.close()

if __name__=="__main__":
   import sys
   maxPWMScoreFileName = sys.argv[1] 
   meanPWMScoreFileName = sys.argv[2]
   numDNaseSites = int(sys.argv[3]) # Number of DNASE sites (NOT genes)
   inputPWMScoreFileNames = []
   for ipsfn in sys.argv[4:len(sys.argv)]:
	# Iterate through input file names and add them to the list
	inputPWMScoreFileNames.append(ipsfn)

   getMaxPWMScores(maxPWMScoreFileName, meanPWMScoreFileName, numDNaseSites, inputPWMScoreFileNames)
		
