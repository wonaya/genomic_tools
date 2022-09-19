import os,sys

#86795	87511
#128, 140 and 150 for Oh43, Tx303
def check_genome(chr):
    print "###calculating total number of possible methylation sites"
    outfile = open("temp_"+str(chr)+".allC", 'w')
    genome_list = []
    whole_list = []
    chr_list = []
    print "reading in genome"
    for genome in open("/work/02114/wonaya/genome/Zea_mays.AGPv2.14/Zea_mays.AGPv2.14.fa", 'r') :
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
    
    chh_list = ["CAA", "CAC", "CAT", "CCA", "CCC", "CCT", "CTA", "CTC", "CTT"]
    chh_rev_list = ["TTG", "GTG", "ATG", "TGG", "GGG", "AGG", "TAG", "GAG", "AAG"]
    chg_list = ["CAG", "CTG", "CCG"]
    chg_rev_list = ["CTG", "CAG", "CGG"]
    cg_file = open("temp_cpg_"+str(chr)+".txt", 'w')
    chg_file = open("temp_chg_"+str(chr)+".txt", 'w')
    chh_file = open("temp_chh_"+str(chr)+".txt", 'w')
    cg_loc_list = []
    chg_loc_list = []
    chh_loc_list = []
    for a in open("temp_"+str(chr)+".bed", 'r') :
        loc = 0
        for base in whole_list[int(chr)-1][int(a.split("\t")[1])-1:int(a.split("\t")[2])] : 
            if base == "C" :
                if whole_list[int(chr)-1][int(a.split("\t")[1])+loc-1:int(a.split("\t")[1])+loc+1] == "CG" :
                    #cg_file.write(str(int(a.split("\t")[1])+loc)+"\n")
                    if int(a.split("\t")[1])+loc not in cg_loc_list :
                        cg_loc_list.append(int(a.split("\t")[1])+loc)
                    #else :
                    #    print "CG", int(a.split("\t")[1])+loc
                    #    print "something not right" ; sys.exit()
                if whole_list[int(chr)-1][int(a.split("\t")[1])+loc-1:int(a.split("\t")[1])+loc+2] in chg_list :
                    #chg_file.write(str(int(a.split("\t")[1])+loc)+"\n")
                    if int(a.split("\t")[1])+loc not in chg_loc_list :
                        chg_loc_list.append(int(a.split("\t")[1])+loc)
                    #else :
                    #    print "CHG", int(a.split("\t")[1])+loc
                    #    print "something not right" ; sys.exit()
                if whole_list[int(chr)-1][int(a.split("\t")[1])+loc-1:int(a.split("\t")[1])+loc+2] in chh_list :
                    if int(a.split("\t")[1])+loc not in chh_loc_list :
                        chh_loc_list.append(int(a.split("\t")[1])+loc)
                    #else :
                    #    print "CHH", int(a.split("\t")[1])+loc
                    #    print "something not right" ; sys.exit()
            elif base == "G" :
                #print whole_list[int(chr)-1][int(a.split("\t")[1])+loc-2:int(a.split("\t")[1])+loc]
                #sys.exit()
                if whole_list[int(chr)-1][int(a.split("\t")[1])+loc-2:int(a.split("\t")[1])+loc] == "CG" :
                    #print int(a.split("\t")[1])+loc
                    #sys.exit()
                    if int(a.split("\t")[1])+loc not in cg_loc_list :
                        cg_loc_list.append(int(a.split("\t")[1])+loc) 
                    #else :
                    #    print "CGr", int(a.split("\t")[1])+loc
                    #    print "something not right" ; sys.exit()
                if whole_list[int(chr)-1][int(a.split("\t")[1])+loc-3:int(a.split("\t")[1])+loc] in chg_rev_list :
                    if int(a.split("\t")[1])+loc not in chg_loc_list :
                        chg_loc_list.append(int(a.split("\t")[1])+loc) 
                    #else :
                    #    print "CHGr", int(a.split("\t")[1])+loc
                    #    print "something not right" ; sys.exit()
                if whole_list[int(chr)-1][int(a.split("\t")[1])+loc-3:int(a.split("\t")[1])+loc] in chh_rev_list :
                    if int(a.split("\t")[1])+loc not in chh_loc_list :
                        chh_loc_list.append(int(a.split("\t")[1])+loc) 
                    #else :
                    #    print "CHHr", int(a.split("\t")[1])+loc
                    #    print "something not right" ; sys.exit()
            loc += 1
        count_cpg = 0
        count_chh = 0
        count_chg = 0
        for chh in chh_list :
            count_chh += whole_list[int(chr)-1][int(a.split("\t")[1])-1:int(a.split("\t")[2])].count(chh)
        for chg in chg_list :
            count_chg += whole_list[int(chr)-1][int(a.split("\t")[1])-1:int(a.split("\t")[2])].count(chg)
        count_cpg += whole_list[int(chr)-1][int(a.split("\t")[1])-1:int(a.split("\t")[2])].count("CG")
        for chh_rev in chh_rev_list :
            count_chh += whole_list[int(chr)-1][int(a.split("\t")[1])-1:int(a.split("\t")[2])].count(chh_rev)
        for chg_rev in chg_rev_list :
            count_chg += whole_list[int(chr)-1][int(a.split("\t")[1])-1:int(a.split("\t")[2])].count(chg_rev)
        count_cpg += whole_list[int(chr)-1][int(a.split("\t")[1])-1:int(a.split("\t")[2])].count("CG")
        outfile.write("\t".join(a.split("\t")[:3])+"\t"+str(count_cpg)+"\t"+str(count_chg)+"\t"+str(count_chh)+"\n")
        
    for cpgs in cg_loc_list :
        cg_file.write(str(cpgs)+"\n")
    for chgs in chg_loc_list :
        chg_file.write(str(chgs)+"\n")
    for chhs in chh_loc_list :
        chh_file.write(str(chhs)+"\n")
    outfile.close()


def check_cov():
    dmr = range(86795,87512)
    print "Oh43"
    cpg_marked = []
    for a in open("Oh43_all3/CpG_context_Oh43_all3_bt202_chr1.bismark.cov", 'r') :
        if int(a.split("\t")[1]) in dmr : 
            cpg_marked.append(int(a.split("\t")[1]))
            print "CpG", a
        if int(a.split("\t")[1]) >= 90000 : 
            break
    chg_marked = []
    for a in open("Oh43_all3/CHG_context_Oh43_all3_bt202_chr1.bismark.cov", 'r') :
        if int(a.split("\t")[1]) in dmr : 
            chg_marked.append(int(a.split("\t")[1]))
            print "CHG", a
        if int(a.split("\t")[1]) >= 90000 : 
            break
    chh_marked = []
    for a in open("Oh43_all3/CHH_context_Oh43_all3_bt202_chr1.bismark.cov", 'r') :
        if int(a.split("\t")[1]) in dmr : 
            chh_marked.append(int(a.split("\t")[1]))
            print "CHH", a
        if int(a.split("\t")[1]) >= 90000 : 
            break
    
    print len(cpg_marked),len(chg_marked),len(chh_marked)
    print set(cpg_marked)&set(chg_marked)
    print set(cpg_marked)&set(chh_marked)
check_genome(1)
sys.exit()  
check_cov()