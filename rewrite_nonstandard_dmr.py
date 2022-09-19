import os,sys

infile = open(sys.argv[1], 'r')
inlines = infile.readlines()
infile.close()
outfile = open(sys.argv[2], 'w')
outfile.write("chr\tstart\tend\twidth\tstrand\tmean.meth.diff\tnum.CpG\tnum.DMC\tDMR.pvalue\tDMR.qvalue\n")
chr = sys.argv[3]
for line in inlines :
    if line.split("\t")[0] != '"seqnames"' and line.split("\t")[1] == '"'+str(chr)+'"' and float(line.split("\t")[10].strip("\n")) <= 0.05 :
        for a in line.split("\t")[1:] :
            if line.split("\t").index(a) == 1 :
                outfile.write(str(chr)+"\t")
            else :
                outfile.write(a.strip('"\n')+"\t")
        outfile.write("\n")
outfile.close()
#os.system("mv "+str(file1ID)+"_"+str(file2ID)+"_"+str(type)+"_dmr_chr"+str(chr)+"_temp.txt "+str(file1ID)+"_"+str(file2ID)+"_"+str(type)+"_dmr_chr"+str(chr)+".txt")
    
