def combinePWMScoreFiles(scoreFileNamePrefix, scoreFileNameListFile):
    # For each enhancer, put the score for each TF in a 2-D array
    scoreFileNameList = open(scoreFileNameListFile)
    combinedScoreArray = []
    firstFile = True
    for scoreFileName in scoreFileNameList:
        # Iterate through files with scores and put them all in an array
        scoreFile = open(scoreFileNamePrefix + scoreFileName[:len(scoreFileName)-1])
        scoreCount = 0
        for line in scoreFile:
            # Iterate through scores and put each score into the array
            if firstFile == True:
                # Initialize the combined score array
                combinedScoreArray.append([])
            lineElements = line.split("\t")
            score = float(lineElements[1])
            combinedScoreArray[scoreCount].append(score)
            scoreCount = scoreCount + 1
        if firstFile == True:
            # Finished with the first file
            firstFile = False
        scoreFile.close()
    scoreFileNameList.close()
    return combinedScoreArray

def writeCombinedScoreArray(combinedScoreArray, combinedScoreFileName):
    # Write the 2-D array of scores to a file
    combinedScoreFile = open(combinedScoreFileName, 'w+')
    for scoreList in combinedScoreArray:
        # Iterate through scores for each enhancer
        for score in scoreList:
            # Iterate through an enhancers scores
            combinedScoreFile.write(str(score))
            combinedScoreFile.write("\t")
        combinedScoreFile.write("\n")
    combinedScoreFile.close()

if __name__=="__main__":
    import sys
    scoreFileNamePrefix = sys.argv[1]
    scoreFileNameListFile = sys.argv[2]
    combinedScoreFileName = sys.argv[3]
                            
    combinedScoreArray = combinePWMScoreFiles(scoreFileNamePrefix, scoreFileNameListFile)
    writeCombinedScoreArray(combinedScoreArray, combinedScoreFileName)
