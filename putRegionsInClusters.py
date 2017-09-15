def putRegionsInClusters(regionsFileName, clustersFileName, clustersOutFileNamePrefix, numClusters):
    # Iterate through clusters and write the regions for each cluster to a file
    clustFileArray = []
    clustSizesArray = []
    for i in range(numClusters):
        # Iterate through clusters and create a file for each cluster
        clustFileName = clustersOutFileNamePrefix + str(i+1) + '.bed'
        clustFile = open(clustFileName, 'w+')
        clustFileArray.append(clustFile)
        clustSizesArray.append(0)
    regionsFile = open(regionsFileName)
    clustersFile = open(clustersFileName)
    for line in regionsFile:
        # Iterate through regions and write each region to the file with the appropriate cluster
        clust = int(clustersFile.readline())
        clustFileArray[clust-1].write(line)
        clustSizesArray[clust-1] = clustSizesArray[clust-1] + 1
    for i in range(numClusters):
        # Iterate through clusters and close the file for each cluster
        clustFileArray[i].close()
    regionsFile.close()
    clustersFile.close()
    return clustSizesArray

def writeClustSizes(clustSizesArray, clustSizesFileName):
    # Write the size of each cluster to a file
    clustSizesFile = open(clustSizesFileName, 'w+')
    for clustSize in clustSizesArray:
        # Write the size of each cluster to a file
        clustSizesFile.write(str(clustSize))
        clustSizesFile.write('\n')

if __name__=="__main__":
    import sys
    regionsFileName = sys.argv[1]
    clustersFileName = sys.argv[2]
    clustersOutFileNamePrefix = sys.argv[3]
    numClusters = int(sys.argv[4])
    clustSizesFileName = sys.argv[5]

    clustSizesArray = putRegionsInClusters(regionsFileName, clustersFileName, clustersOutFileNamePrefix, numClusters)
    writeClustSizes(clustSizesArray, clustSizesFileName)
