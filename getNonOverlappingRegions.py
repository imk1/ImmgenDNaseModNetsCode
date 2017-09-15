def getNonOverlappingRegions(regionListFileNameOne, regionListFileNameTwo, outputFileNameOne, outputFileNameTwo):
	# For 2 lists of regions, make a list of the regions in each list that do not overlap with regions in other list
	regionListFileOne = open(regionListFileNameOne)
	regionListFileTwo = open(regionListFileNameTwo)
	outputFileOne = open(outputFileNameOne, 'w+')
	outputFileTwo = open(outputFileNameTwo, 'w+')
	regionListOne = regionListFileOne.readlines()
	regionListTwo = regionListFileTwo.readlines()
	indexTwo = 0
	lineTwo = regionListTwo[0]
	lineElementsTwo = lineTwo.split("\t")
	chrTwo = lineElementsTwo[0]
	startTwo = int(lineElementsTwo[1])
	endTwo = int(lineElementsTwo[2])
	for lineOne in regionListOne:
		# Iterate through regions in the first list and find those in each list that do not overlap
		lineElementsOne = lineOne.split("\t")
		chrOne = lineElementsOne[0]
		startOne = int(lineElementsOne[1])
		endOne = int(lineElementsOne[2])
		if indexTwo >= len(regionListTwo):
			# All regions in the 2nd list have been considered, so record the region in list 1
			outputFileOne.write(lineOne)
			continue
		while chrOne > chrTwo:
			# Region list 2 has more elements on a chromosome than region list one, so record those elements
			outputFileTwo.write(lineTwo)
			indexTwo = indexTwo + 1
			if indexTwo >= len(regionListTwo):
				# All regions in the second file have been considered, so stop
				break
			lineTwo = regionListTwo[indexTwo]
			lineElementsTwo = lineTwo.split("\t")
			chrTwo = lineElementsTwo[0]
			startTwo = int(lineElementsTwo[1])
			endTwo = int(lineElementsTwo[2])
		while startOne > endTwo:
			# Region list 2 has elements on the chromosome that do not overlap with region list 1, so record
			outputFileTwo.write(lineTwo)
			indexTwo = indexTwo + 1
			if indexTwo >= len(regionListTwo):
				# All regions in the second file have been considered, so stop
				break
			lineTwo = regionListTwo[indexTwo]
			lineElementsTwo = lineTwo.split("\t")
			chrTwo = lineElementsTwo[0]
			startTwo = int(lineElementsTwo[1])
			endTwo = int(lineElementsTwo[2])
		if chrTwo > chrOne:
			# Region list 1 has more elements on a chromosome than region list 2, so record those elements
			outputFileOne.write(lineOne)
			continue
		if startTwo > endOne:
			# Region list 1 has elements on the chromosome that do not overlap with region list 2, so record
			outputFileOne.write(lineOne)
			continue
		indexTwo = indexTwo + 1
		if indexTwo >= len(regionListTwo):
			# All regions in the second file have been considered, so stop
			break
		lineTwo = regionListTwo[indexTwo]
		lineElementsTwo = lineTwo.split("\t")
		chrTwo = lineElementsTwo[0]
		startTwo = int(lineElementsTwo[1])
		endTwo = int(lineElementsTwo[2])
	# DOES NOT HANDLE CASE WHERE A REGION IN 1 LIST OVERLAPS 2 REGIONS IN THE OTHER LIST

if __name__=="__main__":
   import sys
   regionListFileNameOne = sys.argv[1]
   regionListFileNameTwo = sys.argv[2]
   outputFileNameOne = sys.argv[3]
   outputFileNameTwo = sys.argv[4];
                            
   getNonOverlappingRegions(regionListFileNameOne, regionListFileNameTwo, outputFileNameOne, outputFileNameTwo)
