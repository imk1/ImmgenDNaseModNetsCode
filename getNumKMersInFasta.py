import math
import sys as sys
import string


def updateKMerCounts(kMerList,refseq,kMerCounts,merLength):
   #for a given sequence in refseq, scans the sequence and counts the number of each k-mer
   l1 = merLength
   l2 = len(refseq)
   marg = l2 - l1 + 1
   
   if l1>l2:
      return kMerCounts
   for i in range(0,marg):
      currentMer = string.upper(refseq[i:i+l1])
      if currentMer not in kMerList:
	# The current mer is not in the list of k-mers (should have an "N"), so do not consider it
	continue
      kMerCounts[currentMer] = kMerCounts[currentMer] + 1
   return kMerCounts


def getNumKMersInFasta(input_file, kMerListFileName):
   # Get the counts of each k-mer in a fasta file
   kMerListFile = open(kMerListFileName)
   kMerList = []
   for line in kMerListFile:
   	# Make the list of k-mers
   	kMerList.append(line.strip())
   kMerCounts = {}
   for kMer in kMerList:
   	# Initialize the counts of each k-mer to 0
   	kMerCounts[kMer] = 0
   merLength = len(kMerList[0])
   print merLength
   #for each pair of lines of input_file (assuming to be associated with a sequence)
   #scans the sequence for k-mers
   o = open(input_file)
   lines = o.readlines()
   o.close()
   lineIndex = 0
   while lineIndex < len(lines):
      if (lineIndex % 500 == 0) or (lineIndex % 729 == 0):
         print lineIndex
      prefix = lines[lineIndex][1:].strip()
      s = ""
      while lines[lineIndex+1][0] != ">":
         # Iterate through the base sequence and add it to the string
         s = s + lines[lineIndex+1].strip()
         lineIndex = lineIndex + 1
         if lineIndex >= len(lines)-1:
            # At end of file, so stop
            break
      kMerCounts = updateKMerCounts(kMerList,s,kMerCounts,merLength) #sequence should be the last field
      lineIndex = lineIndex + 1
   return [kMerList, kMerCounts]


def writeNumsAndFreqs(kMerList, kMerCounts, numOutputFileName, freqOutputFileName):
	# Write the numbers of occurrences and frequencies of each k-mer to a file
	numOutputFile = open(numOutputFileName, 'w+')
	freqOutputFile = open(freqOutputFileName, 'w+')
	kMerCountsNums = []
	for i in range(len(kMerList)):
		# Iterate through k-mers and make an array of the counts of each
		kMerCountsNums.append(kMerCounts[kMerList[i]])
	sumCounts = sum(kMerCountsNums)
	for count in kMerCountsNums:
		# Iterate through the counts of each k-mer and write the count and the k-mer's frequency to the files
		numOutputFile.write(str(count))
		numOutputFile.write("\n")
		freq = float(count)/float(sumCounts)
		freqOutputFile.write(str(freq))
		freqOutputFile.write("\n")
	numOutputFile.close()
	freqOutputFile.close()


if __name__=="__main__":
   input_file = sys.argv[1]
   kMerListFileName = sys.argv[2]
   numOutputFileName = sys.argv[3]
   freqOutputFileName = sys.argv[4];
                            
   [kMerList, kMerCounts] = getNumKMersInFasta(input_file, kMerListFileName)
   writeNumsAndFreqs(kMerList, kMerCounts, numOutputFileName, freqOutputFileName)
