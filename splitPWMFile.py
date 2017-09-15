def splitPWMFile(PWMFileName, outputFileNamePrefix):
	# Divides a file with multiple PWMs into multiple files in which each file has 1 PWM
	PWMFile = open(PWMFileName)
	currentOutputFile = ""
	firstTFLine = ""
	secondTFLine = ""
	thirdTFLine = ""
	for line in PWMFile:
		# Iterate through the lines of the PWM file, write the non-empty lines to an output file, and create a new output file for every PWM
		if line[0:2] == "TF" and line[3:7] != "Name":
			# At the first line for a new motif, so close the previous file if it exists
			if currentOutputFile != "":
				# Close the file for the previous TF
				currentOutputFile.close()
			firstTFLine = line
		elif line[0:2] == "TF" and line[3:7] == "Name":
			# At the second line for a new motif
			secondTFLine = line
		elif line[0:4] == "Gene":
			# At the third line for a new motif
			thirdTFLine = line
		elif line[0:5] == "Motif":
			# At the fourth line for a new motif, so create a new file for the motif and write the first 4 lines to the file
			currentOutputFileName = outputFileNamePrefix + secondTFLine[8:len(secondTFLine)].strip() + "_" + line[6:len(line)].strip()
			currentOutputFile = open(currentOutputFileName, 'w+')
			currentOutputFile.write(firstTFLine)
			currentOutputFile.write(secondTFLine)
			currentOutputFile.write(thirdTFLine)
			currentOutputFile.write(line)
		elif len(line) > 1:
			# Write the line to the current PWM file
			currentOutputFile.write(line)
	currentOutputFile.close()
	PWMFile.close()

if __name__=="__main__":
   import sys
   PWMFileName = sys.argv[1] 
   outputFileNamePrefix = sys.argv[2]
   splitPWMFile(PWMFileName, outputFileNamePrefix)
