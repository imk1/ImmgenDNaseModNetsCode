def makeTFList(TFFileName):
    # Make a list of TFs in which all of the letters are capitalized
    TFFile = open(TFFileName)
    TFList = []
    for line in TFFile:
        # Iterate through TF file and make a list of all of the TFs in upper case
        TFCap = string.upper(line.strip())
        TFList.append(TFCap)
    return TFList

def findTFRows(TFList, expressionFileName, TFRowFileName):
    # Find the row with the expression data for each TF
    TFRowFile = open(TFRowFileName, 'w+')
    for TF in TFList:
        # Iterate through TFs and find the row with the expression for each
        TFLocation = -1
        expressionFile = open(expressionFileName)
        expressionCount = 0
        for line in expressionFile:
            # Iterate through expressions and find the row with the expression for the TF
            # ASSUMES THAT GENES ARE NOT REPEATED
            lineElements = line.split("\t")
            geneName = string.upper(lineElements[1])
            if geneName == TF:
                # The row with the expression for the TF has been found
                TFLocation = expressionCount
                break
            expressionCount = expressionCount + 1
        expressionFile.close()
        TFRowFile.write(str(TFLocation))
        TFRowFile.write("\n")
    TFRowFile.close()
    
if __name__=="__main__":
    import sys
    import string
    TFFileName = sys.argv[1] 
    expressionFileName = sys.argv[2]
    TFRowFileName = sys.argv[3]
                            
    TFList = makeTFList(TFFileName)
    findTFRows(TFList, expressionFileName, TFRowFileName)
