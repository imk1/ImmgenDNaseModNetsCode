def intersect(listOne, listTwo):
	# Find all of the members of listOne that are not in listTwo
	intersectList = []
	for element in listOne:
		# Iterate through the elements of listOne and add those that are not in listTwo to diffList
		if element in listTwo:
			# Element is not in listTwo, so add it to diffList
			intersectList.append(element)
	return intersectList

def makeList(fileName):
	# Makes a list from the file, where each line in the file is an element in the list
	# ASSUMES THAT EACH LINE IN THE FILE HAS EXACTLY 1 WORD/NUMBER/etc.
	listFile = open(fileName)
	listFromFile = []
	for line in listFile:
		# Iterate through the lines of the file and add each line to the list
		listFromFile.append(string.upper(line.strip()))
	listFile.close()
	return listFromFile

def writeList(listForFile, fileName):
	# Writes a list to a file
	# ASSUMES THAT ALL ELEMENTS IN THE LIST ARE STRINGS
	listFile = open(fileName, 'w+')
	for element in listForFile:
		# Iterate through elements for a list and write each element to the file
		listFile.write(element)
		listFile.write("\n")
	listFile.close()

if __name__=="__main__":
    import sys
    import string
    listOneFileName = sys.argv[1]
    listTwoFileName = sys.argv[2]
    diffListFileName = sys.argv[3]
    listOne = makeList(listOneFileName)
    listTwo = makeList(listTwoFileName)
    intersectList = intersect(listOne, listTwo)
    writeList(intersectList, diffListFileName)
