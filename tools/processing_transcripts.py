import numpy as np
import pandas as pd
import sys

'''
simple python function to process biomart outputs of transcripts
'''
#file = "mart_export_mouse.txt"
file = input("enter file name > ")
fasta = open(file,'r')
print(sys.argv)
lines = fasta.readlines()

count = 0
for line in lines:   
    if ">" in line:
        count += 1

fileOut = np.empty((count,2),dtype=np.object)

count = -1
for line in lines:   
    if ">" in line:
        print(line[1:])
        fileOut[count,0] = line[1:].rstrip()
        count +=1
    else:
        if fileOut[count-1,1] is None:
            fileOut[count-1,1] = line.rstrip()
        else:
            fileOut[count-1,1] = fileOut[count-1,1] + line.rstrip()

pd.DataFrame(fileOut).to_csv("transcripts.tsv")        
