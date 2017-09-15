def filterProbeInfo(probeFileName, probeOutputFileName):
    # Iterate through probe information and remove probes that correspond to unknown places in the genome
    probeFile = open(probeFileName)
    probeOutputFile = open(probeOutputFileName, 'w+')
    for line in probeFile:
        # Iterate through probe information and write the appropriate information to the file
        lineElements = line.split(",")
        if (lineElements[2][1:4] != "chr"): #or ("random" in lineElements[2]):
            # Chromosome is not known
            continue
        if ((lineElements[2][1:4] == "---") or (lineElements[3][1:4] == "---")) or ((lineElements[4][1:4] == "---") or (lineElements[5][1:4] == "---")):
            # Essential information about chromosome is not known
            continue
        probeOutputFile.write(lineElements[1][1:-1])
        probeOutputFile.write("\t")
        probeOutputFile.write(lineElements[2][1:-1])
        probeOutputFile.write("\t")
        probeOutputFile.write(lineElements[3][1:-1])
        probeOutputFile.write("\t")
        probeOutputFile.write(lineElements[4][1:-1])
        probeOutputFile.write("\t")
        probeOutputFile.write(lineElements[5][1:-1])
        probeOutputFile.write("\t")
        probeOutputFile.write(lineElements[7][1:-1])
        probeOutputFile.write("\t")
        probeOutputFile.write(lineElements[8][1:-1])
        probeOutputFile.write("\t")
        probeOutputFile.write(lineElements[9][1:-1])
        probeOutputFile.write("\t")
        probeOutputFile.write(lineElements[10][1:-1])
        probeOutputFile.write("\t")
        probeOutputFile.write(lineElements[11][1:-1])
        probeOutputFile.write("\t")
        probeOutputFile.write(lineElements[12][1:-1])
        probeOutputFile.write("\t")
        probeOutputFile.write(lineElements[13][1:-1])
        probeOutputFile.write("\t")
        probeOutputFile.write(lineElements[14][1:-1])
        probeOutputFile.write("\t")
        probeOutputFile.write(lineElements[15][1:-1])
        probeOutputFile.write("\n")
    probeFile.close()
    probeOutputFile.close()

if __name__=="__main__":
    import sys
    probeFileName = sys.argv[1] 
    probeOutputFileName = sys.argv[2]
                            
    filterProbeInfo(probeFileName, probeOutputFileName)
