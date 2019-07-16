import os,sys
#python 1xnorm.py file.txt reffile.txt windowfile prefix readlength fai

readlength=int(sys.argv[5])
windowfile=sys.argv[3]
prefix=sys.argv[4]

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
    os.system("sortBed -i "+name+".bed -faidx "+sys.argv[6]+" > tmp.bed")
    os.system("mv tmp.bed "+name+".bed")

    ###make into windowsize tiles 
    os.system("bedtools intersect -a %s -b %s -c -sorted > %s"%(windowfile, name+".bed", name+".bedgraph"))
    
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
    ## Calculate scaling factor
    total_reads = 0
    for a in open(name+".bedgraph", 'r') :
        total_reads += float(a.split("\t")[-1].strip("\n"))
    scale = 1/(total_reads*readlength/2.3e9)
    print scale
    # Process each chrom to 1X
    sChroms = sorted(genome.keys())
    for chrom in sChroms:
        for record in genome[chrom]:
	    s, e, v = record
	    oV = v*scale
	    if int(oV) == oV:
		outfile.write('%s\t%i\t%i\t%i\n'%(chrom, s, e, oV))
	    else:
		outfile.write('%s\t%i\t%i\t%.3f\n'%(chrom, s, e, oV))
    outfile.close()
  
## do the same for input file
for files in reffile :
    name=files.split(".bam")[0]
    print files
    os.system("export LC_ALL=C; bedtools bamtobed -i %s | cut -f 1-3 | sort -S 1G -k1,1 -k2,2n > %s"%(files, name+".bed"))
    os.system("sortBed -i "+name+".bed -faidx "+sys.argv[6]+" > tmp.bed")
    os.system("mv tmp.bed "+name+".bed")
    os.system("bedtools intersect -a %s -b %s -c -sorted > %s"%(windowfile, name+".bed", name+".bedgraph"))
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
    ## Calculate scaling factor and normalize
    total_reads = 0
    for a in open(name+".bedgraph", 'r') :
        total_reads += float(a.split("\t")[-1].strip("\n"))
    #print total_reads
    scale = 1/(total_reads*readlength/2.3e9)
    print scale
    # Process each chrom
    sChroms = sorted(genome.keys())
    for chrom in sChroms:
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

filR = open(windowfile, 'r')
linR = filR.readlines()
collect = []
for fil in file :
    file_content = open(fil.split(".bam")[0]+"_1xnorm.bedgraph", 'r')
    lines = file_content.readlines()
    collect.append(lines)

ref_collect = []
for fil in reffile :
    file_content = open(fil.split(".bam")[0]+"_1xnorm.bedgraph", 'r')
    lines = file_content.readlines()
    ref_collect.append(lines)

outfile = open(str(prefix+"_merged.bedgraph", 'w')
index = 0
for lineR in linR : 
    sum = 0
    for file in collect :
        sum += float(file[index].split("\t")[-1].strip("\n"))
    mean = float(sum)/float(len(collect))
    outfile.write("\t".join(lineR.split("\t")[:3]).strip("\n")+"\t")
    outfile.write(str(mean)+"\n")
    index += 1
outfile.close()

outfile = open(str(prefix)+"_mergeref.bedgraph", 'w')
index = 0
for lineR in linR :
    sum = 0
    for file in ref_collect :
        sum += float(file[index].split("\t")[-1].strip("\n"))
    mean = float(sum)/float(len(reffile))
    outfile.write("\t".join(lineR.split("\t")[:3]).strip("\n")+"\t")
    outfile.write(str(mean)+"\n")
    index += 1
outfile.close()

### divide it
outfile = open("ratio_1x_"+str(prefix)+".bedgraph", 'w')

mergef = open(str(prefix)+"_merged.bedgraph",'r')
merger = open(str(prefix)+"_mergeref.bedgraph", 'r')
mergel=mergef.readlines()
mergerl=merger.readlines()
index = 0
for a in mergel :
    outfile.write("\t".join(a.split("\t")[:3]).strip("\n")+"\t")
    if float(mergerl[index].split("\t")[-1].strip("\n")) > 0 :
        outfile.write(str(float(a.split("\t")[-1].strip("\n"))/float(mergerl[index].split("\t")[-1].strip("\n")))+"\n")
    else : 
        outfile.write("0\n")
    index += 1
outfile.close()
### log ratio
#import math
#outfile = open("logratio_1x_"+str(sys.argv[3])+".bedgraph", 'w')
#index = 0
#for a in mergel :
#    outfile.write("\t".join(a.split("\t")[:3]).strip("\n")+"\t")
#    if float(mergerl[index].split("\t")[-1].strip("\n")) > 0 :
#        print math.log(float(a.split("\t")[-1].strip("\n"))/float(mergerl[index].split("\t")[-1].strip("\n")))
#        outfile.write(str(math.log(float(a.split("\t")[-1].strip("\n"))/float(mergerl[index].split("\t")[-1].strip("\n"))))+"\n")
#    else :
#        outfile.write("0\n")
#outfile.close()
os.system("rm -Rf "+str(prefix)+"_merged.bedgraph "+str(prefix)+"_mergeref.bedgraph")
