def getConservation(DNaseFileName, PhastconsFileName, conservationFileName):
	# Gets the average conservation score of each DNase region
	DNaseFile = open(DNaseFileName)
	PhastconsFile = open(PhastconsFileName)
	conservationFile = open(conservationFileName, 'w+')
	PhastconsIndex = 0
	PhastconsChrom = "chr0"
	PhastconsScores = []
	PhastconsLine = ""
	lastChrom = "chr0"
	for line in DNaseFile:
		# Iterate through DNase regions and find the average conservation score for each
		lineElements = line.split("\t")
		chrom = lineElements[0]
		if chrom != lastChrom:
			print chrom
			lastChrom = chrom
		start = int(lineElements[1])
		end = int(lineElements[2])
		stopReached = False
		while chrom > PhastconsChrom:
			# At a new chromosome, so read through Phastcons file until the next chromosome is reached
			newChromReached = False
			while newChromReached == False:
				# Have not reached a new chromosome, so keep reading the Phastcons file
				PhastconsLine = PhastconsFile.readline()
				PhastconsLineElements = PhastconsLine.split(" ")
				if len(PhastconsLine) == 0:
					# There is no more conservation information, so do not record any conservation information for the region
					conservationFile.write(str(-1))
					conservationFile.write("\n")
					stopReached = True
					break
				if len(PhastconsLineElements) > 1:
					# A new chromosome has been reached
					newChromReached = True
					PhastconsChromInfo = PhastconsLineElements[1].split("=")
					PhastconsChrom = PhastconsChromInfo[1]
					PhastconsStartInfo = PhastconsLineElements[2].split("=")
					PhastconsIndex = int(PhastconsStartInfo[1])
					PhastconsLine = PhastconsFile.readline() #ADDED
					PhastconsLineElements = PhastconsLine.split(" ") #ADDED
				else:
					PhastconsIndex = PhastconsIndex + 1
			if stopReached == True:
				# The conservation information for this region cannot be obtained, so continue
				break
		if stopReached == True:
			# The conservation information for this region cannot be obtained, so continue
			continue
		if PhastconsChrom > chrom:
			# DNase conservation information is not in Phastcons file because a new chromosome has been reached in the Phastcons file
			conservationFile.write(str(-1))
			conservationFile.write("\n")
			continue
		while PhastconsIndex < start:
			# Go through bases until the start is reached
			PhastconsLine = PhastconsFile.readline()
			PhastconsLineElements = PhastconsLine.split(" ")
			if len(PhastconsLine) == 0:
				# There is no more conservation information, so do not record any conservation information for the region
				conservationFile.write(str(-1))
				conservationFile.write("\n")
				stopReached = True
				break
			if len(PhastconsLineElements) > 1:
				# At a new starting index for conservation scores
				PhastconsChromInfo = PhastconsLineElements[1].split("=")
				PhastconsChrom = PhastconsChromInfo[1]
				PhastconsStartInfo = PhastconsLineElements[2].split("=")
				if PhastconsChrom != chrom:
					# There is no more conservation information for the rest of the region, so do not compute its conservation
					conservationFile.write(str(-1))
					conservationFile.write("\n")
					PhastconsIndex = int(PhastconsStartInfo[1])
					stopReached = True
					break
				if int(PhastconsStartInfo[1]) != PhastconsIndex:
					# Modify PhastconsIndex appropriately
					PhastconsIndex = int(PhastconsStartInfo[1])
				# ASSUMES THAT THERE IS AT LEAST 1 CONSERVATION SCORE PER HEADING
				PhastconsLine = PhastconsFile.readline()
				PhastconsLineElements = PhastconsLine.split(" ")
			PhastconsIndex = PhastconsIndex + 1
		if stopReached == True:
			# The conservation information for this region cannot be obtained, so continue
			continue
		if PhastconsIndex > start:
			# The DNase region before the last conservation score
			if len(PhastconsScores) < PhastconsIndex - start:
				# The conservation read is after the DNase region and there is no sufficiently early conservation read
				conservationFile.write(str(-1))
				conservationFile.write("\n")
				continue
			firstPhastconsScore = len(PhastconsScores) - (PhastconsIndex - start)
			PhastconsScores = PhastconsScores[firstPhastconsScore:len(PhastconsScores)]
		else:
			PhastconsScores = []
		while PhastconsIndex < end:
			# ASSUMES THAT THERE ARE NO DNASE REGIONS THAT ARE COMPLETELY OVERLAPPING OTHER DNASE REGIONS
			# Gets conservation scores of every base in the DNase region
			PhastconsLine = PhastconsFile.readline()
			PhastconsLineElements = PhastconsLine.split(" ")
			if len(PhastconsLine) == 0:
				# There is no more conservation information, so do not record any conservation information for the region
				conservationFile.write(str(-1))
				conservationFile.write("\n")
				stopReached = True
				break
			if len(PhastconsLineElements) > 1:
				# At a new starting index for conservation scores
				PhastconsChromInfo = PhastconsLineElements[1].split("=")
				PhastconsChrom = PhastconsChromInfo[1]
				PhastconsStartInfo = PhastconsLineElements[2].split("=")
				if (PhastconsChrom != chrom) or (int(PhastconsStartInfo[1]) != PhastconsIndex):
					# There is no more conservation information for the rest of the region, so do not compute its conservation
					conservationFile.write(str(-1))
					conservationFile.write("\n")
					PhastconsIndex = int(PhastconsStartInfo[1])
					stopReached = True
					break
				PhastconsLine = PhastconsFile.readline() #ADDED
				PhastconsLineElements = PhastconsLine.split(" ") #ADDED
				PhastconsScores.append(float(PhastconsLineElements[0])) #ADDED
			else:
				PhastconsScores.append(float(PhastconsLineElements[0]))
				PhastconsIndex = PhastconsIndex + 1
		if stopReached == True:
			# The conservation information for this region cannot be obtained, so continue
			continue
		DNaseLength = end - start
		if len(PhastconsScores) == end - start:
			# There are conservation scores for every base in the DNase region
			averageScore = float(sum(PhastconsScores))/float(len(PhastconsScores))
			conservationFile.write(str(averageScore))
			conservationFile.write("\n")
		else:
			conservationFile.write(str(-1))
			conservationFile.write("\n")
	DNaseFile.close()
	PhastconsFile.close()
	conservationFile.close()

if __name__=="__main__":
    import sys
    DNaseFileName = sys.argv[1]
    PhastconsFileName = sys.argv[2]
    conservationFileName = sys.argv[3]
    getConservation(DNaseFileName, PhastconsFileName, conservationFileName)
