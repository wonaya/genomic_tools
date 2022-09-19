## gc composition calculation
import multiprocessing
import os,sys

def get_genome():
    genome_list = []
    whole_list = []
    chr_list = []
    print "reading in genome"
    for genome in open("/work/02114/wonaya/genome/Zea_mays.AGPv2.14.fa", 'r') :
        if genome[0] == ">" :
            if len(chr_list) == 0 :
                chr_list.append(genome.split(" dna")[0][1:])
            else :
                whole_list.append("".join(genome_list))
                chr_list.append(genome.split(" dna")[0][1:])
                genome_list = []
        else :
            genome_list.append(genome.strip("\n"))
    whole_list.append("".join(genome_list))
    
    return whole_list
    
def gc_calc(chr, window) :
    outfile = open("gc_composition_chr"+str(chr)+".txt", 'w')
    list_window = []
    count = 0
    for x in get_genome()[int(chr)-1] : 
        if len(list_window) == window :
            #print int(list_window.count('G'))+int(list_window.count('C')), int(list_window.count('A'))+int(list_window.count('T'))
            #print float(int(list_window.count('G'))+int(list_window.count('C')))/float(window)
            outfile.write(str(chr))
            outfile.write("\t")
            outfile.write(str(count*window+1))
            outfile.write("\t")
            outfile.write(str(count*window+window))
            outfile.write("\t")
            outfile.write(str(float(int(list_window.count('G'))+int(list_window.count('C')))/float(window)))
            outfile.write("\n")
            list_window = []
            count += 1
        else :
            list_window.append(x)
    outfile.write(str(chr))
    outfile.write("\t")
    outfile.write(str(count*window+1))
    outfile.write("\t")
    outfile.write(str(count*window+len(list_window)))
    outfile.write("\t")
    outfile.write(str(float(int(list_window.count('G'))+int(list_window.count('C')))/len(list_window)))
    outfile.write("\n")
    outfile.close()        
    #print len(get_genome()[int(chr)-1])

jobs = []
for chr in range(1,11) :
    window = 10000
    print chr
    s = multiprocessing.Process(target=gc_calc, args=(str(chr), window ))
    jobs.append(s)
    s.start()
[x.join() for x in jobs]

outfile = open("maize_gc_composition.bedGraph",'w')
for chr in range(1,11):
    for a in open("gc_composition_chr"+str(chr)+".txt", 'r') :
        outfile.write(a)
outfile.close()
