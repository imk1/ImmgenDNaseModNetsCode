import math
import sys as sys

def get_sequence_lengths(input_file, output_file):
   #for each pair of lines of input_file (assuming to be a DNase position associated with a sequence)
   #gets the length of the sequence
   o = open(input_file)
   lines = o.readlines()
   o.close()
   o = open(output_file,'w')
   k = 0
   lineIndex = 0
   while lineIndex < len(lines):
      prefix = lines[lineIndex][5:].strip()
      # Getting length of sequence
      s = 0
      while lines[lineIndex+1][0] != ">":
         # Iterate through the base sequence and add it to the string
         s = s + len(lines[lineIndex+1].strip())
         lineIndex = lineIndex + 1
         if lineIndex >= len(lines)-1:
            # At end of file, so stop
            break
      o.write(prefix+'\t'+str(s)+'\n')
      lineIndex = lineIndex + 1
   o.close()
   
if __name__=="__main__":
   input_file = sys.argv[1]
   output_file = sys.argv[2]
   
   get_sequence_lengths(input_file,output_file)
   

