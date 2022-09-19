import os,sys 
outfile = open("MaizeCyc_2.0.2_annot.txt", 'w')
for a1 in open(sys.argv[1], 'r') :
    for a2 in a1.split("\t")[1:272] :
        if len(a2) > 0 :
            outfile.write(str(a2.split("_")[0])+"\t"+str(a1.split("\t")[0])+"\n") 
outfile.close()