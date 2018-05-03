import os,sys
#python 1xnorm.py file.txt reffile.txt prefix

file = []
reffile = []
for a in open(sys.argv[1], 'r'):
    file.append(a.strip("\n"))
for a in open(sys.argv[2], 'r') :
    reffile.append(a.strip("\n"))

readlen = []
for files in file :
    name=files.split(".bam")[0]
    print files
    os.system("export LC_ALL=C; bedtools bamtobed -i %s | cut -f 1-3 | sort -S 1G -k1,1 -k2,2n > %s"%(files, name+".bed"))
    os.system("bedtools intersect -a %s -b %s -c -sorted > %s"%("maize_4.33_w1000.bed", name+".bed", name+".bedgraph"))
    outfile = open(name+"_1xnorm.bedgraph", 'w')
    from collections import defaultdict 
    from operator import itemgetter
    BGF = open(name+".bedgraph", 'r')
    noComments = filter(lambda x: x[0] != '#', BGF)
    genome = defaultdict(list)
    # Read Data
    for line in noComments:
        tmp = line.rstrip('\n').split('\t')
	s, e, v = int(tmp[1]), int(tmp[2]), float(tmp[3])
	genome[tmp[0]].append((s,e,v))
    # Process each chrom
    sChroms = sorted(genome.keys())
    for chrom in sChroms:
	total = sum(map(itemgetter(2), genome[chrom]))
	nBases = sum(map(lambda x: x[1]-x[0], genome[chrom]))
	scale = nBases/float(total)
	for record in genome[chrom]:
	    s, e, v = record
	    oV = v*scale
	    if int(oV) == oV:
		outfile.write('%s\t%i\t%i\t%i\n'%(chrom, s, e, oV))
	    else:
		outfile.write('%s\t%i\t%i\t%.3f\n'%(chrom, s, e, oV))
    outfile.close()

for files in reffile :
    name=files.split(".bam")[0]
    #os.system("/opt/apps/samtools/1.3.1/bin/samtools index "+files)
    #os.system("/home1/02114/wonaya/.local/bin/bamCoverage --bam "+files+" --binSize 1000 --normalizeTo1x 2300000000 -of bedgraph --numberOfProcessors 16 -o "+name+"_1xnorm.bedgraph")
    print files
    os.system("export LC_ALL=C; bedtools bamtobed -i %s | cut -f 1-3 | sort -S 1G -k1,1 -k2,2n > %s"%(files, name+".bed"))
    os.system("bedtools intersect -a %s -b %s -c -sorted > %s"%("maize_4.33_w1000.bed", name+".bed", name+".bedgraph"))
    outfile = open(name+"_1xnorm.bedgraph", 'w')
    from collections import defaultdict
    from operator import itemgetter
    BGF = open(name+".bedgraph", 'r')
    noComments = filter(lambda x: x[0] != '#', BGF)
    genome = defaultdict(list)
    # Read Data
    for line in noComments:
        tmp = line.rstrip('\n').split('\t')
        s, e, v = int(tmp[1]), int(tmp[2]), float(tmp[3])
        genome[tmp[0]].append((s,e,v))
    # Process each chrom
    sChroms = sorted(genome.keys())
    for chrom in sChroms:
        total = sum(map(itemgetter(2), genome[chrom]))
        nBases = sum(map(lambda x: x[1]-x[0], genome[chrom]))
        scale = nBases/float(total)
        for record in genome[chrom]:
            s, e, v = record
            oV = v*scale
            if int(oV) == oV:
                outfile.write('%s\t%i\t%i\t%i\n'%(chrom, s, e, oV))
            else:
                outfile.write('%s\t%i\t%i\t%.3f\n'%(chrom, s, e, oV))
    outfile.close()
   
print "### intersectBed DONE###"
os.system("rm -Rf tmp*bed")

filR = open("maize_4.33_w1000.bed", 'r')
linR = filR.readlines()
fil1 = open(file[0].split(".bam")[0]+"_1xnorm.bedgraph", 'r')
lin1 = fil1.readlines()
fil2 = open(file[1].split(".bam")[0]+"_1xnorm.bedgraph", 'r')
lin2 = fil2.readlines()
fil3 = open(file[2].split(".bam")[0]+"_1xnorm.bedgraph", 'r')
lin3 = fil3.readlines()
filR1 = open(reffile[0].split(".bam")[0]+"_1xnorm.bedgraph", 'r')
ref1 = filR1.readlines()
filR2 = open(reffile[1].split(".bam")[0]+"_1xnorm.bedgraph", 'r')
ref2 = filR2.readlines()
filR3 = open(reffile[2].split(".bam")[0]+"_1xnorm.bedgraph", 'r')
ref3 = filR3.readlines()

outfile = open(str(sys.argv[3])+"_merged.bedgraph", 'w')
index = 0
for lineR in linR : 
    sum = float(lin1[index].split("\t")[-1].strip("\n"))+float(lin2[index].split("\t")[-1].strip("\n"))+float(lin3[index].split("\t")[-1].strip("\n"))
    mean = float(sum)/float(len(file))
    outfile.write("\t".join(lineR.split("\t")[:3]).strip("\n")+"\t")
    outfile.write(str(mean)+"\n")
    index += 1
outfile.close()

outfile = open(str(sys.argv[3])+"_mergeref.bedgraph", 'w')
index = 0
for lineR in linR :
    sum = float(ref1[index].split("\t")[-1].strip("\n"))+float(ref2[index].split("\t")[-1].strip("\n"))+float(ref3[index].split("\t")[-1].strip("\n"))
    mean = float(sum)/float(len(reffile))
    outfile.write("\t".join(lineR.split("\t")[:3]).strip("\n")+"\t")
    outfile.write(str(mean)+"\n")
    index += 1
outfile.close()

### divide it
outfile = open("ratio_"+str(sys.argv[3])+".bedgraph", 'w')

mergef = open(str(sys.argv[3])+"_merged.bedgraph",'r')
merger = open(str(sys.argv[3])+"_mergeref.bedgraph", 'r')
mergel=mergef.readlines()
mergerl=merger.readlines()
index = 0
for a in mergel :
    outfile.write("\t".join(a.split("\t")[:3]).strip("\n")+"\t")
    if float(mergerl[index].split("\t")[-1].strip("\n")) > 0 :
        outfile.write(str(float(a.split("\t")[-1].strip("\n"))/float(mergerl[index].split("\t")[-1].strip("\n")))+"\n")
    else : 
        outfile.write("0\n")
outfile.close()
os.system("rm -Rf "+str(sys.argv[3])+"_merged.bed "+str(sys.argv[3])+"_mergeref.bed")

