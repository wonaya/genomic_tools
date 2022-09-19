import os,sys

# strandedness save 
pos = []
neg = []
dict = {}
## save strand status
for a in open(sys.argv[2] , 'r') : # example: CpG_context_B73_ACAGTG_merged.sorted_chr10.out bsseq/methylKit/swDMR
    dict[a.split("\t")[3]] = a.split("\t")[1]

# bedgraph 
outlist = []

if sys.argv[1] == "bsseq" :
    for b in open(sys.argv[2].split(".out")[0]+".bedGraph", 'r') :
        outlist.append(b.split("\t")[0].strip("chr")+"\t"+b.split("\t")[1]+"\t"+dict[str(int(b.split("\t")[1])+1)]+"\t"+"CG"+"\t"+b.split("\t")[4]+"\t"+str(int(b.split("\t")[4])+int(b.split("\t")[5].strip("\n"))))
    outfile = open(sys.argv[2].split(".out")[0]+".bsseq", 'w')
    for out in outlist :
        outfile.write(out)
        outfile.write("\n")
    outfile.close()

elif sys.argv[1] == "methylKit" :
    outfile = open(sys.argv[2].split(".out")[0]+".methylKit", 'w')
    for b in open(sys.argv[2].split(".out")[0]+".bedGraph", 'r') :
        if  dict[str(int(b.split("\t")[1])+1)] == "+" :
            strand = "F"
        elif  dict[str(int(b.split("\t")[1])+1)] == "-" :
            strand = "R"
        #print str(b.split("\t")[0])+"."+str(b.split("\t")[1]), str(b.split("\t")[0]), str(b.split("\t")[1]), strand
        total = float(int(b.split("\t")[4])+int(b.split("\t")[5].strip("\n")))
        freqc = float(int(b.split("\t")[4]))
        freqt = float(int(b.split("\t")[5].strip("\n")))
        percc = freqc/total*100
        perct = freqt/total*100
        #print str(b.split("\t")[0])+"."+str(b.split("\t")[1]), str(b.split("\t")[0]), str(b.split("\t")[1]), strand, int(total), percc, perct
        outfile.write(str(b.split("\t")[0])+"."+str(b.split("\t")[1]))
        outfile.write("\t")
        outfile.write(str(b.split("\t")[0])+"\t"+str(b.split("\t")[1]))
        outfile.write("\t")
        outfile.write(str(strand)+"\t"+str(int(total))+"\t"+str(percc)+"\t"+str(perct)+"\n")
    outfile.close()

elif sys.argv[1] == "swDMR" :
    outfile = open(sys.argv[2].split(".")[0]+".swDMR", 'w')
    for b in open(sys.argv[2].split(".")[0]+".bedGraph", 'r') :
        chr = str(b.split("\t")[0])
        pos = str(b.split("\t")[1])
        strand = str(dict[str(int(b.split("\t")[1])+1)])
        context = "CG"
        methno = str(b.split("\t")[4])
        unmethno = str(b.split("\t")[5].strip("\n"))
        outfile.write(chr+"\t"+pos+"\t"+strand+"\t"+context+"\t"+methno+"\t"+unmethno+"\n")
    outfile.close()