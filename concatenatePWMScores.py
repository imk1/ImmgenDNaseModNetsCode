def concatenatePWMScores(PWMScoreFileNameList, PWMScoresFileName):
	# Concatenate all of the PWM scores for all of the TFs into 1 file
	PsfnlFile = open(PWMScoreFileNameList)
	PWMScoresFile = open(PWMScoresFileName, 'w+')
	for line in PsfnlFile:
		# Iterate through TFs and add each one's PWM scores for all of the DNase sites to the output file
		print line
		PWMScoreFileSingleTF = open(line.strip())
		for PWMScore in PWMScoreFileSingleTF:
			# Iterate through PWM scores for the different DNase sites and write each to the output file
			PWMScoreElements = PWMScore.split("\t")
			PWMScoreVal = PWMScoreElements[1].strip()
			PWMScoresFile.write(PWMScoreVal)
			PWMScoresFile.write("\t")
		PWMScoresFile.write("\n")
		PWMScoreFileSingleTF.close()
	PsfnlFile.close()
	PWMScoresFile.close()

if __name__=="__main__":
   import sys
   PWMScoreFileNameList = sys.argv[1] 
   PWMScoresFileName = sys.argv[2]
   concatenatePWMScores(PWMScoreFileNameList, PWMScoresFileName)
