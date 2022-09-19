import os,sys

if str(sys.argv[1]).split(".")[-1] != "gff3" :
    print "place .gff file"
    sys.exit()

type = []
for a in open(sys.argv[1], 'r') :
    if a[0] != "#" :
        type.append(a.split("\t")[2])

type = set(type)

outfile = open(str(sys.argv[1]).split(".")[0]+".bed", 'w')
for a in open(sys.argv[1], 'r') :
    if a[0] != "#" and a.split("\t")[2] == "gene" :
        outfile.write(a.split("\t")[0])
        outfile.write("\t")
        outfile.write(a.split("\t")[3])
        outfile.write("\t")
        outfile.write(a.split("\t")[4])
        outfile.write("\t")
        ### change this as needed
        #outfile.write(a.split("\t")[-1])
        outfile.write(a.split(";")[0].split("=")[1])
        outfile.write("\n")
outfile.close()
sys.exit()
## split ratio_segmentation ## comment out if needed
timing = ['ES', 'ESMS', 'ESLS', 'ESMSLS','MS','MSLS','ESMSLS']

for time in timing :
    print time
    outfile = open(time+".bed", 'w')
    for a in open(str(sys.argv[1]).split(".")[0]+".bed", 'r') :
        if a.split("\t")[-1].strip("\n") == time :
            for x in range(int(a.split("\t")[1][:-3]), int(a.split("\t")[2][:-3])+1) : 
                print x
            #print a.split("\t")[1][:-3], a.split("\t")[2][:-3]
            sys.exit()
            #outfile.write(a)
    outfile.close()

