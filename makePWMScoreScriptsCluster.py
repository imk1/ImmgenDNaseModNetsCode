def makePWMScoreScriptsCluster(PWMFileNamePrefix, PWMFileNameListFile, DNaseFileNameListFile, scoreFileNamePrefix, scriptFileNamePrefix, cutoff):
    # Make a script that will run each of the PWM score commands
    DNaseFileNameList = open(DNaseFileNameListFile)
    PWMFileNameList = open(PWMFileNameListFile)
    PWMFileNames = PWMFileNameList.readlines()
    for line in DNaseFileNameList:
        # Iterate through DNase file names and make entries in the scrpt for each
        lineElements = line.split("/")
        DNaseCellLine = lineElements[len(lineElements)-1][:len(lineElements[len(lineElements)-1])-8]
        for PWMFN in PWMFileNames:
            # Iterate through names of PWM files and make an entry in the script for each
            scriptFile = open(scriptFileNamePrefix+DNaseCellLine + "_" + PWMFN[:len(PWMFN)-5] + ".sh", 'w+')
            fullFN = PWMFileNamePrefix + PWMFN[:len(PWMFN)-1]
            scriptFile.write("python /afs/cs.stanford.edu/u/imk1/ImmgenProject/src/get_pwm_score_plus.py ")
            scriptFile.write(fullFN + " " + line[:len(line)-2] + " ")
            scriptFile.write(scoreFileNamePrefix + PWMFN[:len(PWMFN)-5] + "_" + DNaseCellLine + "_Scores.txt " + str(cutoff))
            scriptFile.close()
    DNaseFileNameList.close()
    PWMFileNameList.close()
    
if __name__=="__main__":
   import sys as sys
   PWMFileNamePrefix = "/afs/cs.stanford.edu/u/imk1/ImmgenProject/PWMsJASPAR/"
   PWMFileNameListFile = "C:/Users/Irene/Documents/ImmgenProject/PWMsJASPAR/PWMsJasparFileNames"
   DNaseFileNameListFile = "C:/Users/Irene/Documents/ImmgenProject/MouseDNase/mouseDNaseFastaFileNamesPlus.txt"
   scoreFileNamePrefix = "/afs/cs.stanford.edu/u/imk1/ImmgenProject/PWMExpCutNormScores/"
   scriptFileNamePrefix = "C:/Users/Irene/Documents/ImmgenProject/src/run_PwmScore_"
   cutoff = 1
                            
   makePWMScoreScriptsCluster(PWMFileNamePrefix, PWMFileNameListFile, DNaseFileNameListFile, scoreFileNamePrefix, scriptFileNamePrefix, cutoff)
