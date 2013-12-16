### using multiple bam files, counts reads into intervals specified and calculate differences between input ChIP
### 12/06/2013 by Jawon Song

import os,sys
import math
import numbers
import readline, glob
def complete(text, state):
    return (glob.glob(text+'*')+[None])[state]
readline.set_completer_delims(' \t\n;')
readline.parse_and_bind("tab: complete")
readline.set_completer(complete)
 
####### Setup the run #######
 
specie = raw_input("Type specie of interest (eg. maize or arabidopsis) : ")
specie_val = None
while specie_val == None : 
    if specie == "maize" or specie == "arabidopsis" :
        specie_val = specie
    else :
        specie = raw_input("Type specie of interest (maize or arabidopsis) : ")
specie = specie_val

window = raw_input("Type window size (bp) in integer (eg. 1000) : ")
window_val = None
while window_val == None :
    try: 
        int(window)
        window_val = window
    except ValueError: 
        window = raw_input("Type window size (bp) in integer (eg. 1000) : ")
window = window_val

noofbam = raw_input("Type in number of bam files to compute (including control IP) : ")
noofbam_val = None
while noofbam_val == None :
    try :
        int(noofbam)
        noofbam_val = noofbam
    except ValueError: 
        noofbam = raw_input("Type in number of bam files to compute (including control IP) : ")
noofbam = noofbam_val

list_of_bams = []
for x in range(1,int(noofbam)) :
    list_of_bams.append(raw_input("Type in name of file "+str(x)+" (eg. 1_S1_L001_sorted.bam) : "))
    x += 1
list_of_bams.append(raw_input("Type in name of the control bam file (eg. 5_S5_L001_sorted.bam) : "))

#outname = raw_input("Name of MultiBamCov output file : ")
location_of_mbc = raw_input("Full path of multiBamCov (eg. /usr/bin/bedtools/bin/multiBamCov) : ")
location_of_samtools = raw_input("Full path of Samtools (eg. /usr/bin/samtools) : ")

####### Begin running #######
print "\n### Begin running ###\n"

outfile = open(str(specie)+"_"+str(window)+"bp.gff", 'w')
outfile.write("##gff-version 3\n")

print "###generating GFF file for "+str(specie)+" in window size "+str(window)+", please wait ..."
if specie == "arabidopsis" :
    chromosome_size = {1:30600,2:19800,3:23600,4:18700,5:27100}
    max_chr = 6
elif specie == "maize" :
    chromosome_size = {1:301600,2:238000,3:232200,4:241500,5:217000,6:169200,7:176800,8:175800,9:156740,10:150180}
    max_chr = 11
    
for a in range(1, max_chr) : 
    for b in range(0, chromosome_size[a]) :
        diff = 10 - len(str(b*1000+1))
        ID = '0'*diff
        outfile.write(str(a)+"\t.\tWindow\t"+str(b*1000+1)+"\t"+str((b+1)*1000)+"\t.\t.\t.\tID="+str(a)+str(ID+str(b*1000+1))+"\n")
outfile.close()

scr = str(location_of_mbc)+" -bams "+str(" ".join(list_of_bams))+" -bed "+str(specie)+"_"+str(window)+"bp.gff > tmp.out"
print "###running "+str(scr)+" ..., (takes a while)"
os.system(scr)

list_of_reads_mapped = []
print "###running samtools to get number of reads mapped ..."
for x in range(1,int(noofbam)+1) :
    flagscr = str(location_of_samtools)+" view -c -F 4 "+list_of_bams[x-1]+" > tmp.flagstat"
    os.system(flagscr)
    for flag in open("tmp.flagstat", 'r') :
        list_of_reads_mapped.append(int(flag.strip("\n")))
print list_of_reads_mapped
os.system("rm tmp.flagstat")

factor = {}
for count in list_of_reads_mapped[:-1] :
    factor[list_of_reads_mapped.index(count)+1] = float(list_of_reads_mapped[-1])/float(count)
#factor = {1:float(float(6914218)/float(6322434)), 2:float(float(6914218)/float(5918587)), 3:float(float(6914218)/float(6145466)), 4:float(float(6914218)/float(5399519))}

print "###writing bedGraphs ..."
for bams in list_of_bams[:-1] :
    front = bams.split(".")[0]
    back = str(list_of_bams[-1]).split(".")[0]
    outfile = open(str(front)+"_"+str(back)+".bedGraph", 'w')
    for a in open("tmp.out", 'r') :
        ip = float(a.split("\t")[9+int(list_of_bams.index(bams))])
        control = float(a.split("\t")[9+int(len(list_of_bams[:-1]))].strip("\n"))
        k = factor[int(list_of_bams.index(bams))+1]
        if ip > 0 :
            outfile.write("Chr"+str(a.split("\t")[0])+"\t")
            outfile.write(str(a.split("\t")[3])+"\t")
            outfile.write(str(a.split("\t")[4])+"\t")
            outfile.write(str(math.log(((ip*k)/(control+1)),2))+"\n")    
    outfile.close()
os.system("rm -Rf tmp.out")
print "###All calculations complete!"