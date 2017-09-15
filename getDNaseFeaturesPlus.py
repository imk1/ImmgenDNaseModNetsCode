def makeTFList(TFFileName):
    # Make a list of TFs in which all of the letters are capitalized
    TFFile = open(TFFileName)
    TFList = []
    for line in TFFile:
        # Iterate through TF file and make a list of all of the TFs in upper case
        TFCap = string.upper(line.strip())
        TFList.append(TFCap)
    return TFList

def getDNaseFeaturesPlus(genesFileName, finalGeneListFileName, DNaseFileName, genesToDNaseFileName, hypersensFileName, conservFileName, genesToHypersensFileName, genesToConservFileName, genesToDistFileName, genesToUpFileName, distCutoff):
    # Find the DNase features for each DNase site that is associated witih each gene
    # ASSUMES THAT CONSERVATION IS THE ONLY FEATURE THAT CAN HAVE A MISSING VALUE
    # Assumes that the genes-to-DNase file is 0-INDEXED
    genesFile = open(genesFileName)
    finalGeneList = makeTFList(finalGeneListFileName)
    genesFile.readline()
    startListArray = []
    lastGeneName = ""
    geneIndex = -1
    for line in genesFile:
	# Iterate through genes and make a list of each gene's TSS's
	geneElements = line.split("\t")
	geneName = string.upper(geneElements[12].strip()[1:-1])
	if geneName not in finalGeneList:
		# The gene was removed, so do not consider it
		continue
	if geneName != lastGeneName:
		# At a new gene
		startListArray.append([])
		geneIndex = geneIndex + 1
		lastGeneName = geneName
	start = int(geneElements[4].strip())
        startListArray[geneIndex].append(start)
    genesFile.close()
    print len(startListArray)

    DNaseFile = open(DNaseFileName)
    DNaseMiddleArray = []
    for line in DNaseFile:
	# Iterate through DNase hypersensitivity regions and make a list of the middle of each
	DNaseElements = line.split("\t")
	DNaseStart = int(DNaseElements[1])
    	DNaseEnd = int(DNaseElements[2])
    	DNaseMiddle = float(DNaseEnd + DNaseStart + 1) / float(2) # The middle of the DNase hypersensitivity region is the average of the start and the end (the true start is 1 base passed the listed start)
	DNaseMiddleArray.append(DNaseMiddle)
    DNaseFile.close()
    
    hypersensFile = open(hypersensFileName)
    hypersensArray = []
    for line in hypersensFile:
	# Make an array of hypersensitivities
	hypersens = float(line.strip())
	if hypersens >= 5000:
		# Hypersensitivity is very large, so set it to 1
		hypersensArray.append(1)
	else:
		hypersensArray.append(hypersens/float(5000))
    hypersensFile.close()

    conservFile = open(conservFileName)
    conservArray = []
    for line in conservFile:
	# Make an array of conservation scores
	conservArray.append(float(line.strip()))
    conservFile.close()
    

    genesToDNaseFile = open(genesToDNaseFileName)
    genesToHypersensFile = open(genesToHypersensFileName, 'w+')
    genesToConservFile = open(genesToConservFileName, 'w+')
    genesToDistFile = open(genesToDistFileName, 'w+')
    genesToUpFile = open(genesToUpFileName, 'w+')
    geneCount = 0

    for line in genesToDNaseFile:
	# Iterate through genes and get the features for each corresponding DNase
	lineElements = line.split("\t")
	geneName = lineElements[0]
	genesToHypersensFile.write(geneName)
	genesToHypersensFile.write("\t")
	genesToConservFile.write(geneName)
	genesToConservFile.write("\t")
	genesToDistFile.write(geneName)
	genesToDistFile.write("\t")
	genesToUpFile.write(geneName)
	genesToUpFile.write("\t")

	if len(lineElements) < 5:
		# No DNase for the current gene
		genesToHypersensFile.write("\n")
		genesToConservFile.write("\n")
		genesToDistFile.write("\n")
		genesToUpFile.write("\n")
		geneCount = geneCount + 1
		continue

	for i in range(4, len(lineElements)):
		# Iterate through DNase sites and compute their features
		DNaseIndex = int(lineElements[i])
		DNaseConserv = conservArray[DNaseIndex]
		if DNaseConserv != -1:
			# The conservation is not missing
			genesToHypersensFile.write(str(hypersensArray[DNaseIndex]))
			genesToConservFile.write(str(DNaseConserv))
			DNaseMiddle = DNaseMiddleArray[DNaseIndex]
			distArray = []
			for start in startListArray[geneCount]:
				# Find the distances from each TSS to the current DNase
				dist = abs(start - DNaseMiddle)
				distArray.append(dist)
			minDist = min(distArray)
			minDistFeat = 1 - float(minDist)/float(distCutoff) # Do 1 - normalized distance so that closer DNase hypersensitivity sites have higher distance feature values
			genesToDistFile.write(str(minDistFeat))
			minDistIndex = distArray.index(minDist)
			trueDist = startListArray[geneCount][minDistIndex] - DNaseMiddle
			isUp = 1
			if trueDist < 0:
				# DNase hypersensitivity site is downstream
				isUp = 0
			genesToUpFile.write(str(isUp))
			
		else:
			# Conservation information is missing, so make all of the feature values 0
			genesToHypersensFile.write(str(0))
			genesToConservFile.write(str(0))
			genesToDistFile.write(str(0))
			genesToUpFile.write(str(0))

		genesToHypersensFile.write("\t")
		genesToConservFile.write("\t")
		genesToDistFile.write("\t")
		genesToUpFile.write("\t")

	genesToHypersensFile.write("\n")
	genesToConservFile.write("\n")
	genesToDistFile.write("\n")
	genesToUpFile.write("\n")
	geneCount = geneCount + 1

    genesToHypersensFile.close()
    genesToConservFile.close()
    genesToDistFile.close()
    genesToUpFile.close()
    genesToDNaseFile.close()
	

if __name__=="__main__":
    import sys
    import string
    genesFileName = sys.argv[1]
    finalGeneListFileName = sys.argv[2]
    DNaseFileName = sys.argv[3]
    genesToDNaseFileName = sys.argv[4]
    hypersensFileName = sys.argv[5]
    conservFileName = sys.argv[6]
    genesToHypersensFileName = sys.argv[7]
    genesToConservFileName = sys.argv[8]
    genesToDistFileName = sys.argv[9]
    genesToUpFileName = sys.argv[10]
    distCutoff = int(sys.argv[11])
                        
    getDNaseFeaturesPlus(genesFileName, finalGeneListFileName, DNaseFileName, genesToDNaseFileName, hypersensFileName, conservFileName, genesToHypersensFileName, genesToConservFileName, genesToDistFileName, genesToUpFileName, distCutoff)
