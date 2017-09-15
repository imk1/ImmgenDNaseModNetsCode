def findCommonProbes(expressionProbesFileName, annotationProbesFileName, annotationProbesOutFileName):
    # Remove information on annotated probes whose expression data is not present
    expressionProbesFile = open(expressionProbesFileName)
    annotationProbesFile = open(annotationProbesFileName)
    annotationProbesOutFile = open(annotationProbesOutFileName, 'w+')
    annotationLines = annotationProbesFile.readlines()
    annotationLinesCount = 0
    currentLineElements = annotationLines[annotationLinesCount].split("\t")
    currentProbe = int(currentLineElements[0])
    for line in expressionProbesFile:
        # Iterate through probes with expression data and find corresponding annotation data
        lineElements = line.split("\t")
        probe = int(lineElements[0])
        while currentProbe < probe:
            # Iterate through annotations until the correct probe has been found
            annotationLinesCount = annotationLinesCount + 1
            if annotationLinesCount >= len(annotationLines):
                print "Too many probes!"
                print str(probe)
                break
            currentLineElements = annotationLines[annotationLinesCount].split("\t")
            currentProbe = int(currentLineElements[0])
        if probe == currentProbe:
            # The correct probe has been found, so write its annotation to a file
            annotationProbesOutFile.write(annotationLines[annotationLinesCount])
        else:
            print "Problem!"
            print str(probe)
    expressionProbesFile.close()
    annotationProbesFile.close()
    annotationProbesOutFile.close()

if __name__=="__main__":
    import sys
    expressionProbesFileName = sys.argv[1]
    annotationProbesFileName = sys.argv[2]
    annotationProbesOutFileName = sys.argv[3]
                            
    findCommonProbes(expressionProbesFileName, annotationProbesFileName, annotationProbesOutFileName)
