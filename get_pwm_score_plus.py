# FIX CASE FOR ACGT!!!!!!!!!!!!!!

import math
import sys as sys

def make_pwm(txt,pseudo_count):
   #takes an input PFM file (in JASPAR format) and outputs a PWM array ...
   #txt: is the name of the PFM file (4 rows, ACGT respectively, the number of columns corresponds ...
   #to the length of the motfif
   #psuedo_count can be set to a small value (e.g. 0.001) to avoid over-fitting
       

   o = open(txt)
   lines = o.readlines()
   o.close()

   seq = 'ACGT'
   mat = {}
   total_scores = {}
    
   for j in range(0,len(lines[0].strip().split())):
       total_scores[j] = 0
       for i in range(0,4):
           total_scores[j] +=  int(lines[i].strip().split()[j])

   for i in range(0,4):
      mat[i] = []
      s = lines[i].strip().split()
      for j in range(0,len(s)):
         mat[i].append((pseudo_count + float(s[j]))/float(4*pseudo_count + total_scores[j]))
 
   return mat

def get_max_pwm(mat):
   #takes a PWM array and gets the maximum possible score

   maxScores = []
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
   for i in range(0, len(mat[0])):
      # Compute the total maximum score
      total_score += math.log(maxScores[i],2) - math.log(0.25,2)
   return total_score


def get_pwm_score(mat,refseq):
   #for a given PWM matrix in mat, and a given sequence in refseq, scans mat across refseq
        
   seq = 'ACGT'
   idx = {}
   for i in range(0,4):
      idx[seq[i]] = i


   l1 = len(mat[0])
   l2 = len(refseq)
   marg = l2 - l1
   
   if l1>l2:
      return {}

   total_score = {}
   for i in range(0,marg):
      total_score[i] = 0
      for j in range(0,l1):
         try:
            total_score[i] += math.log( mat[idx[refseq[i+j]]][j],2) - math.log(0.25,2)
         except:
            continue
   return total_score

def get_pwm_score_exp(mat,refseq,cutoff):
   #for a given PWM matrix in mat, and a given sequence in refseq, scans mat across refseq
        
   seq = 'ACGT'
   idx = {}
   for i in range(0,4):
      idx[seq[i]] = i


   l1 = len(mat[0])
   l2 = len(refseq)
   marg = l2 - l1
   
   if l1>l2:
      return {}

   total_score = {}
   for i in range(0,marg):
      total_score[i] = 0
      for j in range(0,l1):
         try:
            total_score[i] += math.log( mat[idx[refseq[i+j]]][j],2) - math.log(0.25,2)
         except:
            continue
   for i in range(0,marg):
      # Iterate through scores and take the exp of each
      if total_score[i] >= cutoff:
         # Score is sufficiently high to consider it
         total_score[i] = 2 ** total_score[i]
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
   marg = l2 - l1
   
   if l1>l2:
      return {}
   total_score = {}
   zeroScoreList = []
   for i in range(0,marg):
      total_score[i] = 0
      for j in range(0,l1):
         try:
            total_score[i] += math.log( mat[idx[refseq[i+j]]][j],2) - math.log(0.25,2)
         except:
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

def get_pwm_score_dynamic(mat,refseq):
   #same as above, except the window to scan PWM is dynamic based on the length of the motif:
   #the window ONLY includes all bases that intersect 1 position with the center of refseq

   seq = 'ACGT'
   idx = {}
   for i in range(0,4):
      idx[seq[i]] = i


   l1 = len(mat[0])
   l2 = len(refseq)
   snppos = int((len(refseq)-1)/2.0)
   refseqD = refseq[snppos-l1+1:snppos+l1]
   marg = len(refseqD) - l1
   if l1>l2:
      return {}

   total_score = {}
   for i in range(0,l1):
      total_score[i] = 0
      for j in range(0,l1):
         try:
            total_score[i] += math.log( mat[idx[refseqD[i+j]]][j],2) - math.log(0.25,2)
         except:
            continue

   return total_score



def reverse_seq(seq):
   map = {}
   map['A'] = 'T'
   map['T'] = 'A'
   map['C'] = 'G'
   map['G'] = 'C'
   map['N'] = 'N'
   map['-'] = 'N'
   revseq = ''
   for i in seq:
      revseq += map[i]
   return revseq


def score_each_dnase_window(input_file,mat,output_file):
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
      prefix = lines[lineIndex][5:].strip()
      s = ""
      while lines[lineIndex+1][0] != ">":
         # Iterate through the base sequence and add it to the string
         s = s + lines[lineIndex+1].strip()
         lineIndex = lineIndex + 1
         if lineIndex >= len(lines)-1:
            # At end of file, so stop
            break
      r = get_pwm_score(mat,s) #sequence should be the last field
      #print s[0]+'\t'+str(max(r.values()))
      maxScore = max(r.values())
      if maxScore < 0:
         # Set all scores that are less than 0 to 0
         maxScore = 0
      o.write(prefix+'\t'+str(maxScore)+'\n')
      lineIndex = lineIndex + 1
   o.close()

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
      prefix = lines[lineIndex][5:].strip()
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

def score_each_snp_window_mat(input_file,mat):
    #same as above but returs matrix
    o = open(input_file)
    lines = o.readlines()
    o.close()
    scores = {}
    
    k = 0
    for l in lines:
        s = l.strip().split()
        r = get_pwm_score(mat,s[-1]) #sequence should be the last field
        scores[s[0]] = (s[-1],max(r.values()))
    return scores

if __name__=="__main__":

# FIX CASE FOR ACGT!!!!!!!!!!!!!!
# FIX CASE FOR TRY/EXCEPT!!!!!!!!

   motif = sys.argv[1] 
   input_file = sys.argv[2]
   output_file = sys.argv[3]
   cutoff = int(sys.argv[4])
   pseudo_count = 0;
                            
   mat = make_pwm(motif,pseudo_count)
   #score_each_dnase_window(input_file,mat,output_file)
   score_each_dnase_window_norm(input_file,mat,output_file,cutoff)
   

