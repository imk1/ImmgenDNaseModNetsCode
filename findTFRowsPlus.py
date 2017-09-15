def makeTFList(TFFileName):
    # Make a list of TFs in which all of the letters are capitalized
    TFFile = open(TFFileName)
    TFList = []
    for line in TFFile:
        # Iterate through TF file and make a list of all of the TFs in upper case
        TFCap = string.upper(line.strip())
        TFList.append(TFCap)
    return TFList

def findTFRows(TFList, geneFileName, TFRowFileName):
    # Find the row with the expression data for each TF by finding the row with the gene name in the final list of genes, where the gene order in this list is the same as the gene order in the expression file
    TFRowFile = open(TFRowFileName, 'w+')
    for TF in TFList:
        # Iterate through TFs and find the row with the expression for each
        TFLocation = -1
        geneFile = open(geneFileName)
        count = 0
        for line in geneFile:
            # Iterate through expressions and find the row with the expression for the TF
            # ASSUMES THAT GENES ARE NOT REPEATED
            geneName = string.upper(line.strip())
            if geneName == TF:
                # The row with the expression for the TF has been found, add 1 so that rows are not 0-indexed
                TFLocation = count + 1
                break
            count = count + 1
        geneFile.close()
        TFRowFile.write(str(TFLocation))
        TFRowFile.write("\n")
    TFRowFile.close()
    
if __name__=="__main__":
    import sys
    import string
    TFFileName = sys.argv[1] 
    geneFileName = sys.argv[2]
    TFRowFileName = sys.argv[3]
                            
    TFList = makeTFList(TFFileName)
    findTFRows(TFList, geneFileName, TFRowFileName)
