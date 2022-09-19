import os,sys
import numpy as np

chrom_list = []
seq_list = []
genome_list = []
for a in open(sys.argv[1],'r') :
    if a[0] == ">" :
        #save new chrom name 
        chrom_list.append(a.split(" dna")[0][1:])
        if len(chrom_list) > 1 : 
            #merge sequences
            genome_list.append("".join(seq_list))
            #reset seq_list
            seq_list = []
    else :
        seq_list.append(a.strip("\n"))

#deal with last chrom sequence
genome_list.append("".join(seq_list))
outfile = open("AGPv4_cg_sites_test.bed", 'w') 
outfile.write("chr\tstart\tend\t#cg sites\t#chg sites\t#chh sites\n")

for chrom in chrom_list :
    cg_loc_list = []
    loc = 0
    #print genome_list[0][800:900], len(genome_list[0][300:400])
    #sys.exit()
    for base in genome_list[chrom_list.index(chrom)] :
        loc += 1
        if base == "C" :
            if genome_list[chrom_list.index(chrom)][loc-1:loc+1] == "CG" :
                cg_loc_list.append(loc)
        elif base == "G" :
            if genome_list[chrom_list.index(chrom)][loc-2:loc] == "CG" :
                cg_loc_list.append(loc)
    #print cg_loc_list
    #sys.exit()
    loc = 0
    chg_loc_list = []
    chg_list = ["CAG", "CTG", "CCG"]
    chg_rev_list = ["CTG", "CAG", "CGG"]
    for base in genome_list[chrom_list.index(chrom)] :
        loc += 1
        if base == "C" :
            if genome_list[chrom_list.index(chrom)][loc-1:loc+2] in chg_list  :
                chg_loc_list.append(loc)
        elif base == "G" :
            if genome_list[chrom_list.index(chrom)][loc-3:loc] in chg_rev_list :
                chg_loc_list.append(loc)

    loc = 0
    chh_loc_list = []
    chh_list = ["CAA", "CAC", "CAT", "CCA", "CCC", "CCT", "CTA", "CTC", "CTT"]
    chh_rev_list = ["TTG", "GTG", "ATG", "TGG", "GGG", "AGG", "TAG", "GAG", "AAG"]
    for base in genome_list[chrom_list.index(chrom)] :
        loc += 1
        if base == "C" :
            if genome_list[chrom_list.index(chrom)][loc-1:loc+2] in chh_list  :
                chh_loc_list.append(loc)
        elif base == "G" :
            if genome_list[chrom_list.index(chrom)][loc-3:loc] in chh_rev_list :
                chh_loc_list.append(loc)

    count_list = np.zeros(((len(genome_list[chrom_list.index(chrom)])/100)+2,3))
    for cg_loc in cg_loc_list :
        count_list[(int(cg_loc)-1)/100][0] += 1
    for chg_loc in chg_loc_list :
        count_list[(int(chg_loc)-1)/100][1] += 1
    for chh_loc in chh_loc_list :
        count_list[(int(chh_loc)-1)/100][2] += 1

    for x in range(0, len(count_list)) :
        outfile.write(str(chrom)+"\t")
        outfile.write(str(x*100+1)+"\t")
        outfile.write(str((x+1)*100)+"\t")
        listval = count_list[x].tolist()
        outfile.write(str(int(listval[0]))+"\t")
        outfile.write(str(int(listval[1]))+"\t")
        outfile.write(str(int(listval[2]))+"\n")
        #if x == 10 :  
        #    sys.exit()
outfile.close()
    
    
