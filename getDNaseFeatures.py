def getDNaseFeatures(genesFileName, DNaseFileName, genesToDNaseFileName, hypersensFileName, conservFileName, genesToHypersensFileName, genesToConservFileName, genesToDistFileName, genesToUpFileName):
    # Find the DNase features for each DNase site that is associated witih each gene
    genesFile = open(genesFileName)
    genesFile.readline()
    startListArray = []
    lastGeneName = ""
    geneIndex = -1
    for line in genesFile:
	# Iterate through genes and make a list of each gene's TSS's
	geneElements = line.split("\t")
	geneName = geneElements[12].strip()[1:-1]
	if geneName != lastGeneName:
		# At a new gene
		startListArray.append([])
		geneIndex = geneIndex + 1
		lastGeneName = geneName
	start = int(geneElements[4].strip())
        startListArray[geneIndex].append(start)
    genesFile.close()

    DNaseFile = open(DNaseFileName)
    DNaseMiddleArray = []
    for line in DNaseFile:
	# Iterate through DNase hypersensitivity regions and make a list of the middle of each
	DNaseElements = line.split("\t")
	DNaseStart = int(DNaseElements[1])
    	DNaseEnd = int(DNaseElements[2])
    	DNaseMiddle = DNaseStart + (float(DNaseEnd - DNaseStart) / float(2))
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
		genesToHypersensFile.write(str(hypersensArray[DNaseIndex]))
		genesToHypersensFile.write("\t")
		genesToConservFile.write(str(conservArray[DNaseIndex]))
		genesToConservFile.write("\t")
		DNaseMiddle = DNaseMiddleArray[DNaseIndex]
		distArray = []
		for start in startListArray[geneCount]:
			# Find the distances from each TSS to the current DNase
			dist = abs(start - DNaseMiddle)
			distArray.append(dist)
		minDist = min(distArray)
		minDistFeat = float(minDist)/float(5000)
		genesToDistFile.write(str(minDistFeat))
		genesToDistFile.write("\t")
		minDistIndex = distArray.index(minDist)
		trueDist = startListArray[geneCount][minDistIndex] - DNaseMiddle
		isUp = 1
		if trueDist < 0:
			# DNase hypersensitivity site is downstream
			isUp = 0
		genesToUpFile.write(str(isUp))
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
    genesFileName = sys.argv[1]
    DNaseFileName = sys.argv[2]
    genesToDNaseFileName = sys.argv[3]
    hypersensFileName = sys.argv[4]
    conservFileName = sys.argv[5]
    genesToHypersensFileName = sys.argv[6]
    genesToConservFileName = sys.argv[7]
    genesToDistFileName = sys.argv[8]
    genesToUpFileName = sys.argv[9]
                        
    getDNaseFeatures(genesFileName, DNaseFileName, genesToDNaseFileName, hypersensFileName, conservFileName, genesToHypersensFileName, genesToConservFileName, genesToDistFileName, genesToUpFileName)
