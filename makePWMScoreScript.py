def makePWMScoreScript(PWMFileNamePrefix, PWMFileNameListFile, DNaseFileNameListFile, scoreFileNamePrefix, scriptFileName, cutoff, pseudocount):
    # Make a script that will run all of the PWM score commands
    DNaseFileNameList = open(DNaseFileNameListFile)
    PWMFileNameList = open(PWMFileNameListFile)
    PWMFileNames = PWMFileNameList.readlines()
    scriptFile = open(scriptFileName, 'w+')
    for line in DNaseFileNameList:
        # Iterate through DNase file names and make entries in the scrpt for each
        lineElements = line.split("/")
        #DNaseCellLine = lineElements[len(lineElements)-1][:len(lineElements[len(lineElements)-1])-8]
        for PWMFN in PWMFileNames:
            # Iterate through names of PWM files and make an entry in the script for each
            fullFN = PWMFileNamePrefix + PWMFN.strip()
            scriptFile.write("python get_pwm_score_cisbp.py ")
            scriptFile.write(fullFN + " " + line.strip() + " ")
            #scriptFile.write(scoreFileNamePrefix + PWMFN[:len(PWMFN)-5] + "_" + DNaseCellLine + "_Scores.txt " + str(cutoff))
	    scriptFile.write(scoreFileNamePrefix + PWMFN.strip() + "_Scores.txt " + str(cutoff) + " " + str(pseudocount))
            scriptFile.write("\n")
    DNaseFileNameList.close()
    PWMFileNameList.close()
    scriptFile.close()
    
if __name__=="__main__":
   import sys as sys
   PWMFileNamePrefix = "/scr/ImmgenProject/CisBP_2012_05_07/"
   PWMFileNameListFile = "/scr/ImmgenProject/CisBP_2012_05_07/fileNamesPWMsFiltered"
   DNaseFileNameListFile = "/scr/ImmgenProject/MouseDNase/mouseDNaseFastaFileNamesPlus.txt"
   scoreFileNamePrefix = "/scr/ImmgenProject/PWMExpCutNormScores/"
   scriptFileName = "/scr/ImmgenProject/src/runGetPwmScoreCisBP.sh"
   cutoff = -10
   pseudocount = .0001
                            
   makePWMScoreScript(PWMFileNamePrefix, PWMFileNameListFile, DNaseFileNameListFile, scoreFileNamePrefix, scriptFileName, cutoff, pseudocount)
