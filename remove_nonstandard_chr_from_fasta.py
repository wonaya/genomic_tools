import os,sys

outfile = open(sys.argv[1].split(".fa")[0]+"_rmscaffold.fa",'w')
seq = []
for a in open(sys.argv[1],'r') :
    if a[0] == ">" :
        if len(seq) == 0 :
            seq.append(a)
        else :
            if seq[0][:2] == ">C" or seq[0][:2] == ">M": 
                #write
                print seq[0].split(" ")
                for b in seq :
                    outfile.write(b)
            seq = []
            seq.append(a)
    else :
        seq.append(a) 

if "dna_sm:chromosome" in seq[0].split(" ") :
    for b in seq :
        outfile.write(b)
outfile.close()
