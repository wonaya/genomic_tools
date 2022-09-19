import os,sys
import numpy as np

def rm_scaffold(input_fa) :
    list = []
    chr_list = range(1,11)
    chr_list_str = map(str, chr_list)
    sequence = []
    for a in open(input_fa) : 
        if a[1:2] in chr_list_str : 
            if len(sequence) == 0 :
                sequence.append(a)
            else :
                ## write output and reset
                chr = str(sequence[0].split(" ")[0][1:])
                outfile = open("chr"+chr+".txt", 'w') 
                for sequence_lines in sequence :
                    outfile.write(sequence_lines)
                outfile.close()
                sequence = []
                sequence.append(a)
        else: 
            sequence.append(a)
            if a[:1] == ">" :
                chr = str(sequence[0].split(" ")[0][1:])
                outfile = open("chr"+chr+".txt", 'w') 
                for sequence_lines in sequence :
                    outfile.write(sequence_lines)
                outfile.close()
		break          
                
rm_scaffold(sys.argv[1])

chr_list = range(1,11)
chr_list_str = map(str, chr_list)
outfile = open(sys.argv[2], 'w')
for chr in chr_list_str: 
    print chr
    for a in open("chr"+chr+".txt", 'r') :
        outfile.write(a)
outfile.close()

