
import math
import sys as sys
import string

def make_pwm_cisbp(txt,pseudo_count):
   #takes an input PWM file (in cis-bp format) and outputs a PWM array ...
   #txt: is the name of the PWM file (7 descriptor rows, 1 column with numbers and 4 columns for ACGT respectively, the number additional rows corresponds ...
   #to the length of the motfif
   #psuedo_count can be set to a small value (e.g. 0.001) to avoid over-fitting
       

   o = open(txt)
   lines = o.readlines()[7:] # Remove descriptor lines
   o.close()

   seq = 'ACGT'
   mat = {}
   total_scores = {}
    
   for j in range(0,len(lines)):
       # Iterate through positions and compute the total score for each position
       total_scores[j] = 0
       for i in range(0,4):
	   # Add the score for each base to the total score
           total_scores[j] +=  float(lines[j].strip().split()[i+1])

   for i in range(0,4):
      # Compute the score for each base at each position by dividing the score for the base by the total score and incorporating pseudocounts
      mat[i] = []     
      for j in range(0,len(lines)):
	 # Iterate through positions and compute the score for each position
	 s = lines[j].strip().split()
         mat[i].append((float(pseudo_count + float(s[i+1])))/float(4*pseudo_count + total_scores[j]))
   print len(mat[0])
   return mat


def get_max_pwm(mat):
   #takes a PWM array and gets the maximum possible score

   maxScores = []
   print "Getting max!"
   for i in range(0,4):
      # Iterate through bases for each position
      for j in range(0, len(mat[0])):
         # Iterate through positions and bases
         if i == 0:
            # Initialize score
            maxScores.append(-10000)
         if mat[i][j] >= maxScores[j]:
            # Replace the score with the current largest score
            maxScores[j] = mat[i][j]
   total_score = 0
   for j in range(0, len(mat[0])):
      # Compute the total maximum score
      total_score += math.log(maxScores[j],2)
   print total_score
   return total_score


def get_max_pwm_base(mat, baseFreqs):
   #takes a PWM array and gets the maximum possible score
   #accounts for base frequencies

   seq = 'ACGT'
   maxScores = []
   print "Getting max!"
   for i in range(0,4):
      # Iterate through bases for each position
      for j in range(0, len(mat[0])):
         # Iterate through positions and bases
         if i == 0:
            # Initialize score
            maxScores.append(-10000)
	 currentScore = float(mat[i][j])/float(baseFreqs[seq[i]])
         if currentScore >= maxScores[j]:
            # Replace the score with the current largest score
            maxScores[j] = currentScore
   total_score = 0
   for j in range(0, len(mat[0])):
      # Compute the total maximum score
      total_score += math.log(maxScores[j],2)
   print total_score
   return total_score


def get_pwm_score_exp(mat,refseq,cutoff):
   #for a given PWM matrix in mat, and a given sequence in refseq, scans mat across refseq
        
   seq = 'ACGT'
   idx = {}
   for i in range(0,4):
      idx[seq[i]] = i

   l1 = len(mat[0])
   l2 = len(refseq)
   marg = l2 - l1 + 1
   
   if l1>l2:
      return {}
   total_score = {}
   zeroScoreList = []
   for i in range(0,marg):
      total_score[i] = 0
      for j in range(0,l1):
         try:
            total_score[i] += math.log( mat[idx[string.upper(refseq[i+j])]][j],2)
         except:
	    total_score[i] -= 2 # Subract the average of the possible scores
            continue
      if total_score[i] < cutoff:
	# Score is NOT sufficiently high to consider it
        total_score[i] = -10000
	zeroScoreList.append(i)
   for i in range(0,marg):
      # Iterate through scores and take the exp of each
      if i not in zeroScoreList:
         # Score is sufficiently high to consider it
         total_score[i] = 2 ** total_score[i]
	 if total_score[i] > 1:
	    print str(total_score[i])
      else:
         total_score[i] = 0
   return total_score


def get_pwm_score_exp_norm(mat,refseq,cutoff,max_total_score):
   #for a given PWM matrix in mat, and a given sequence in refseq, scans mat across refseq
        
   seq = 'ACGT'
   idx = {}
   for i in range(0,4):
      idx[seq[i]] = i

   l1 = len(mat[0])
   l2 = len(refseq)
   marg = l2 - l1 + 1
   
   if l1>l2:
      return {}
   total_score = {}
   zeroScoreList = []
   for i in range(0,marg):
      total_score[i] = 0
      for j in range(0,l1):
         try:
            total_score[i] += math.log( mat[idx[string.upper(refseq[i+j])]][j],2)
         except:
	    if (((string.upper(refseq[i+j]) == 'A') or (string.upper(refseq[i+j]) == 'C')) or ((string.upper(refseq[i+j]) == 'G') or (string.upper(refseq[i+j]) == 'T'))) and mat[idx[string.upper(refseq[i+j])]][j] != 0:
	    	print "Problem!"
		print mat[idx[string.upper(refseq[i+j])]][j]
	    total_score[i] -= 2 # Subract the average of the possible scores
	    #total_score[i] -= 10000 # Subtract large quantity so that score does not make cutoff and region is therefore not considered
            continue
      if total_score[i] >= cutoff:
	# Score is sufficiently high to consider it
      	total_score[i] -= max_total_score
      else:
        total_score[i] = -10000
	zeroScoreList.append(i)
   for i in range(0,marg):
      # Iterate through scores and take the exp of each
      if i not in zeroScoreList:
         # Score is sufficiently high to consider it
         total_score[i] = 2 ** total_score[i]
	 if total_score[i] > 1:
	    print str(total_score[i])
      else:
         total_score[i] = 0
   return total_score


def get_pwm_score_exp_base(mat,refseq,cutoff,max_total_score,baseFreqs):
   #for a given PWM matrix in mat, and a given sequence in refseq, scans mat across refseq
   #divides PWM probabilities by base frequencies
        
   seq = 'ACGT'
   idx = {}
   for i in range(0,4):
      idx[seq[i]] = i

   l1 = len(mat[0])
   l2 = len(refseq)
   marg = l2 - l1 + 1
   
   if l1>l2:
      return {}
   total_score = {}
   zeroScoreList = []
   for i in range(0,marg):
      total_score[i] = 0
      for j in range(0,l1):
         try:
            currentBase = string.upper(refseq[i+j])
            total_score[i] += (math.log( mat[idx[currentBase]][j],2) - math.log( baseFreqs[currentBase],2))
         except:
	    if (((string.upper(refseq[i+j]) == 'A') or (string.upper(refseq[i+j]) == 'C')) or ((string.upper(refseq[i+j]) == 'G') or (string.upper(refseq[i+j]) == 'T'))) and mat[idx[string.upper(refseq[i+j])]][j] != 0:
	    	print "Problem!"
		print mat[idx[string.upper(refseq[i+j])]][j]
	    total_score[i] -= 2 # Subract the average of the possible scores
            continue
      #if total_score[i] >= cutoff:
	# Score is sufficiently high to consider it
      	#total_score[i] -= max_total_score
      #else:
      if total_score[i] < cutoff:
        total_score[i] = -10000
	zeroScoreList.append(i)
   for i in range(0,marg):
      # Iterate through scores and take the exp of each
      if i not in zeroScoreList:
         # Score is sufficiently high to consider it
         total_score[i] = 2 ** total_score[i]
	 #if total_score[i] > 1.0001:
	 #   print str(total_score[i])
      else:
         total_score[i] = 0
   return total_score


def score_each_dnase_window_sum(input_file,mat,output_file,cutoff):
   #for each pair of lines of input_file (assuming to be a DNase position associated with a sequence)
   #scans the motif in mat across, output_file is the name of the output file
   o = open(input_file)
   lines = o.readlines()
   o.close()
   scores = {}
   o = open(output_file,'w')
   k = 0
   lineIndex = 0
   while lineIndex < len(lines):
      prefix = lines[lineIndex][1:].strip()
      s = ""
      while lines[lineIndex+1][0] != ">":
         # Iterate through the base sequence and add it to the string
         s = s + lines[lineIndex+1].strip()
         lineIndex = lineIndex + 1
         if lineIndex >= len(lines)-1:
            # At end of file, so stop
            break
      r = get_pwm_score_exp(mat,s,cutoff) #sequence should be the last field
      #print s[0]+'\t'+str(max(r.values()))
      sumScore = sum(r.values())
      if sumScore < 0:
         # Set all scores that are less than 0 to 0
         sumScore = 0
      o.write(prefix+'\t'+str(sumScore)+'\n')
      lineIndex = lineIndex + 1
   o.close()


def score_each_dnase_window_norm(input_file,mat,output_file,cutoff):
   #for each pair of lines of input_file (assuming to be a DNase position associated with a sequence)
   #scans the motif in mat across, output_file is the name of the output file
   o = open(input_file)
   lines = o.readlines()
   o.close()
   scores = {}
   o = open(output_file,'w')
   k = 0
   lineIndex = 0
   max_total_score = get_max_pwm(mat)
   while lineIndex < len(lines):
      prefix = lines[lineIndex][1:].strip()
      s = ""
      while lines[lineIndex+1][0] != ">":
         # Iterate through the base sequence and add it to the string
         s = s + lines[lineIndex+1].strip()
         lineIndex = lineIndex + 1
         if lineIndex >= len(lines)-1:
            # At end of file, so stop
            break
      r = get_pwm_score_exp_norm(mat,s,cutoff,max_total_score) #sequence should be the last field
      #print s[0]+'\t'+str(max(r.values()))
      sumScore = sum(r.values())
      if sumScore < 0:
         # Set all scores that are less than 0 to 0
         sumScore = 0
      o.write(prefix+'\t'+str(sumScore)+'\n')
      lineIndex = lineIndex + 1
   o.close()


def score_each_dnase_window_base(input_file,mat,output_file,cutoff,baseFreqFileName):
   #for each pair of lines of input_file (assuming to be a DNase position associated with a sequence)
   #scans the motif in mat across, output_file is the name of the output file
   #divdes each score for each base by the base frequency and normalizes each score for a sequence by the maximum score
   baseFreqFile = open(baseFreqFileName)
   baseFreqs = {}
   seq = 'ACGT'
   lineCount = 0
   for line in baseFreqFile:
      # Iterate through frequencies of bases and add each frequency to the array of base frequencies
      baseFreqs[seq[lineCount]] = float(line.strip())
      lineCount = lineCount + 1
   o = open(input_file)
   lines = o.readlines()
   o.close()
   scores = {}
   o = open(output_file,'w')
   k = 0
   lineIndex = 0
   max_total_score = get_max_pwm_base(mat, baseFreqs)
   while lineIndex < len(lines):
      prefix = lines[lineIndex][1:].strip()
      s = ""
      while lines[lineIndex+1][0] != ">":
         # Iterate through the base sequence and add it to the string
         s = s + lines[lineIndex+1].strip()
         lineIndex = lineIndex + 1
         if lineIndex >= len(lines)-1:
            # At end of file, so stop
            break
      r = get_pwm_score_exp_base(mat,s,cutoff,max_total_score,baseFreqs) #sequence should be the second field
      sumScore = sum(r.values())
      if sumScore < 0:
         # Set all scores that are less than 0 to 0
	 print "Problem! -- sumScore < 0"
         sumScore = 0
      o.write(prefix+'\t'+str(sumScore)+'\n')
      lineIndex = lineIndex + 1
   o.close()


if __name__=="__main__":
   motif = sys.argv[1] 
   input_file = sys.argv[2]
   output_file = sys.argv[3]
   cutoff = int(sys.argv[4])
   pseudo_count = float(sys.argv[5])
   if len(sys.argv) > 6:
      # There is an additional input argument because base frequencies will be accounted for
      baseFreqFileName = sys.argv[6]
                            
   mat = make_pwm_cisbp(motif,pseudo_count)
   # SWITCH SCORING METHOD TO NORM FOR NO BASES AND TO SUM FOR NO BASES AND NO NORMALIZATION!
   score_each_dnase_window_base(input_file,mat,output_file,cutoff,baseFreqFileName) # Use 0 for cutoff (1, but log space)
   #score_each_dnase_window_norm(input_file,mat,output_file,cutoff)
   #score_each_dnase_window_sum(input_file,mat,output_file,cutoff)
   

