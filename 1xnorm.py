#Updated 6/5/19
#USAGE python 1xnorm.py file.txt reffile.txt fai windowsize prefix
#file.txt contains list of test .bam files e.g. ES1.bam\nES2.bam\nES3.bam
#reffile.txt contains list of reference .bam files e.g. control1.bam\ncontrol2.bam

#import individual file names
import os,sys
file = []
reffile = []
for a in open(sys.argv[1], 'r'):
    file.append(a.strip("\n"))
for a in open(sys.argv[2], 'r') :
    reffile.append(a.strip("\n"))

### make genome windows
window = open("tmp.window", 'w')
for a in open(sys.argv[3], 'r') :
    start = 0
    for x in range(0, int(float(a.split("\t")[1])/float(sys.argv[4]))) :
        window.write(a.split("\t")[0]+"\t")
        window.write(str(x*float(sys.argv[4]))+"\t")
        window.write(str((x+1)*float(sys.argv[4]))+"\n")
    window.write(a.split("\t")[0]+"\t"+str((x+1)*float(sys.argv[4]))+"\t"+a.split("\t")[1]+"\n")
window.close()

### extract name of each sample from file.txt file names
readlen = []
for files in file :
    name=files.split(".bam")[0]
    print files
    # run bed tools to convert bam file to bed file and extract columns and sort
    os.system("export LC_ALL=C; bedtools bamtobed -i %s | cut -f 1-3 | sort -S 1G -k1,1 -k2,2n > %s"%(files, name+".bed"))
    # overlap with maize genome interval to output to file chr1\n1\1000\nES etc. 
    os.system("bedtools intersect -a %s -b %s -c -sorted > %s"%("tmp.window", name+".bed", name+".bedgraph"))
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
    ## Calculate scaling factor to normalize
    total_reads = 0
    for a in open(name+".bedgraph", 'r') :
        total_reads += float(a.split("\t")[-1].strip("\n"))
    scale = 1/(total_reads*130/2.3e9) ### given read length=130, genome size 2.3e9
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
for files in reffile :
    name=files.split(".bam")[0]
    print files
    #os.system("export LC_ALL=C; bedtools bamtobed -i %s | cut -f 1-3 | sort -S 1G -k1,1 -k2,2n > %s"%(files, name+".bed"))
    #os.system("bedtools intersect -a %s -b %s -c -sorted > %s"%("maize_4.33_w1000.bed", name+".bed", name+".bedgraph"))
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
    scale = 1/(total_reads*130/2.3e9)
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

filR = open("maize_4.33_w1000.bed", 'r')
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

outfile = open(str(sys.argv[3])+"_merged.bedgraph", 'w')
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

outfile = open(str(sys.argv[3])+"_mergeref.bedgraph", 'w')
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
outfile = open("ratio_1x_"+str(sys.argv[3])+".bedgraph", 'w')

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
    index += 1
outfile.close()
os.system("rm -Rf "+str(sys.argv[3])+"_merged.bedgraph "+str(sys.argv[3])+"_mergeref.bedgraph")
