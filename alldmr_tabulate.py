import os,sys

context = sys.argv[1]
chr_dict = {1:302000000,2:240000000,3:233000000,4:243000000,5:219000000,6:171000000,7:178000000,8:178000000,9:158000000,10:152000000}


for chr in range(1,11) :
    print chr, context
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
    #chg_list = ["CAG", "CTG", "CCG"]
    #chg_rev_list = ["CTG", "CAG", "CGG"]
    chh_list = ["CAA", "CAC", "CAT", "CCA", "CCC", "CCT", "CTA", "CTC", "CTT"]
    chh_rev_list = ["TTG", "GTG", "ATG", "TGG", "GGG", "AGG", "TAG", "GAG", "AAG"]
    
    large_list = []
    for x in range(0,chr_dict[chr]):
        large_list.append([0]*15)
    print "read B73"
    for a in open(str(context)+"_B73_all3_R1_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
    print "read Mo17"
    for a in open(str(context)+"_Mo17_all3_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][3] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][4] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][5] += 1
    print "read Oh43"
    for a in open(str(context)+"_Oh43_all3_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][6] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][7] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][8] += 1
    print "read CML322"
    for a in open(str(context)+"_CML322_all3_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][9] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][10] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][11] += 1
    print "read Tx303"
    for a in open(str(context)+"_Tx303_all3_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][12] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][13] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][14] += 1
    
    outfile = open("temp_alldmr_withcov_"+str(chr)+"_"+str(context)+".txt",'w')
    for a in open("/corral-tacc/tacc/iplant/vaughn/springer_vaughn/eichten/jawon/alldmr/temp_alldmr_"+str(chr)+"_"+str(context)+".txt",'r') :
        if context == "CpG" :
            valid_cg = []
            loc = 0
            for base in whole_list[int(chr)-1][int(a.split("\t")[2])-1:int(a.split("\t")[3].strip("\n"))+1] : 
                if base == "C" :
                    if whole_list[int(chr)-1][int(a.split("\t")[2])+loc-1:int(a.split("\t")[2])+loc+1] == "CG" :
                        if int(a.split("\t")[2])+loc not in valid_cg :
                            valid_cg.append(int(a.split("\t")[2])+loc)
                elif base == "G" :
                    if whole_list[int(chr)-1][int(a.split("\t")[2])+loc-2:int(a.split("\t")[2])+loc] == "CG" :
                        if int(a.split("\t")[2])+loc not in valid_cg :
                            valid_cg.append(int(a.split("\t")[2])+loc) 
                loc += 1
        
        elif context == "CHG" :
            valid_cg = []
            loc = 0
            for base in whole_list[int(chr)-1][int(a.split("\t")[2])-1:int(a.split("\t")[3].strip("\n"))+1] : 
                if base == "C" :
                    if whole_list[int(chr)-1][int(a.split("\t")[2])+loc-1:int(a.split("\t")[2])+loc+2] in chg_list :
                        if int(a.split("\t")[2])+loc not in valid_cg :
                            valid_cg.append(int(a.split("\t")[2])+loc)
                elif base == "G" :
                    if whole_list[int(chr)-1][int(a.split("\t")[2])+loc-3:int(a.split("\t")[2])+loc] in chg_rev_list :
                        if int(a.split("\t")[2])+loc not in valid_cg :
                            valid_cg.append(int(a.split("\t")[2])+loc) 
                loc += 1

        elif context == "CHH" :
            valid_cg = []
            loc = 0
            for base in whole_list[int(chr)-1][int(a.split("\t")[2])-1:int(a.split("\t")[3].strip("\n"))+1] : 
                if base == "C" :
                    if whole_list[int(chr)-1][int(a.split("\t")[2])+loc-1:int(a.split("\t")[2])+loc+2] in chh_list :
                        if int(a.split("\t")[2])+loc not in valid_cg :
                            valid_cg.append(int(a.split("\t")[2])+loc)
                elif base == "G" :
                    if whole_list[int(chr)-1][int(a.split("\t")[2])+loc-3:int(a.split("\t")[2])+loc] in chh_rev_list :
                        if int(a.split("\t")[2])+loc not in valid_cg :
                            valid_cg.append(int(a.split("\t")[2])+loc) 
                loc += 1
        
        countmet_b73 = 0
        countunmet_b73 = 0
        uniquec_b73 = 0
        if a.split("\t")[1] == "Mo17" :
            countmet_mo17 = 0
            countunmet_mo17 = 0
            uniquec_mo17 = 0
            for x in set(range(int(a.split("\t")[2]),(int(a.split("\t")[3].strip("\t"))+1)))&set(valid_cg) :
                countmet_b73 += large_list[x][0]
                countunmet_b73 += large_list[x][1]
                uniquec_b73 += large_list[x][2]
                countmet_mo17 += large_list[x][3]
                countunmet_mo17 += large_list[x][4]
                uniquec_mo17 += large_list[x][5]
            outfile.write(a.strip("\n"))
            outfile.write("\t")
            if float(countmet_b73)+float(countunmet_b73) == 0 :
                outfile.write("NA")
            else :
                outfile.write(str(float(countmet_b73)/(float(countmet_b73)+float(countunmet_b73))*100))
            outfile.write("\t")
            outfile.write(str(uniquec_b73))
            outfile.write("\t")
            if float(countmet_mo17)+float(countunmet_mo17) == 0 :
                outfile.write("NA")
            else :
                outfile.write(str(float(countmet_mo17)/(float(countmet_mo17)+float(countunmet_mo17))*100))
            outfile.write("\t")
            outfile.write(str(uniquec_mo17))
            outfile.write("\n")
        elif a.split("\t")[1] == "Oh43" :
            countmet_oh43 = 0
            countunmet_oh43 = 0
            uniquec_oh43 = 0
            for x in set(range(int(a.split("\t")[2]),(int(a.split("\t")[3].strip("\t"))+1)))&set(valid_cg) :
                countmet_b73 += large_list[x][0]
                countunmet_b73 += large_list[x][1]
                uniquec_b73 += large_list[x][2]
                countmet_oh43 += large_list[x][6]
                countunmet_oh43 += large_list[x][7]
                uniquec_oh43 += large_list[x][8]
            outfile.write(a.strip("\n"))
            outfile.write("\t")
            if float(countmet_b73)+float(countunmet_b73) == 0 :
                outfile.write("NA")
            else :
                outfile.write(str(float(countmet_b73)/(float(countmet_b73)+float(countunmet_b73))*100))
            outfile.write("\t")
            outfile.write(str(uniquec_b73))
            outfile.write("\t")
            if float(countmet_oh43)+float(countunmet_oh43) == 0 :
                outfile.write("NA")
            else :
                outfile.write(str(float(countmet_oh43)/(float(countmet_oh43)+float(countunmet_oh43))*100))
            outfile.write("\t")
            outfile.write(str(uniquec_oh43))
            outfile.write("\n")
        elif a.split("\t")[1] == "CML322" :
            countmet_cml322 = 0
            countunmet_cml322 = 0
            uniquec_cml322 = 0
            for x in set(range(int(a.split("\t")[2]),(int(a.split("\t")[3].strip("\t"))+1)))&set(valid_cg) :
                countmet_b73 += large_list[x][0]
                countunmet_b73 += large_list[x][1]
                uniquec_b73 += large_list[x][2]
                countmet_cml322 += large_list[x][9]
                countunmet_cml322 += large_list[x][10]
                uniquec_cml322 += large_list[x][11]
            outfile.write(a.strip("\n"))
            outfile.write("\t")
            if float(countmet_b73)+float(countunmet_b73) == 0 :
                outfile.write("NA")
            else :
                outfile.write(str(float(countmet_b73)/(float(countmet_b73)+float(countunmet_b73))*100))
            outfile.write("\t")
            outfile.write(str(uniquec_b73))
            outfile.write("\t")
            if float(countmet_cml322)+float(countunmet_cml322) == 0 :
                outfile.write("NA")
            else :
                outfile.write(str(float(countmet_cml322)/(float(countmet_cml322)+float(countunmet_cml322))*100))
            outfile.write("\t")
            outfile.write(str(uniquec_cml322))
            outfile.write("\n")
        elif a.split("\t")[1] == "Tx303" :
            countmet_tx303 = 0
            countunmet_tx303 = 0
            uniquec_tx303 = 0
            for x in set(range(int(a.split("\t")[2]),(int(a.split("\t")[3].strip("\t"))+1)))&set(valid_cg) :
                countmet_b73 += large_list[x][0]
                countunmet_b73 += large_list[x][1]
                uniquec_b73 += large_list[x][2]
                countmet_tx303 += large_list[x][12]
                countunmet_tx303 += large_list[x][13]
                uniquec_tx303 += large_list[x][14]
                outfile.write(a.strip("\n"))
            outfile.write("\t")
            if float(countmet_b73)+float(countunmet_b73) == 0 :
                outfile.write("NA")
            else :
                outfile.write(str(float(countmet_b73)/(float(countmet_b73)+float(countunmet_b73))*100))
            outfile.write("\t")
            outfile.write(str(uniquec_b73))
            outfile.write("\t")
            if float(countmet_tx303)+float(countunmet_tx303) == 0 :
                outfile.write("NA")
            else :
                outfile.write(str(float(countmet_tx303)/(float(countmet_tx303)+float(countunmet_tx303))*100))
            outfile.write("\t")
            outfile.write(str(uniquec_tx303))
            outfile.write("\n")
    outfile.close()
    del large_list


outfile = open("temp_alldmr_withcov_"+str(context)+".txt",'w')
for chr in range(1,11) :
    print chr
    for a in open("temp_alldmr_withcov_"+str(chr)+"_"+str(context)+".txt",'r') :
        outfile.write(str(chr)+"\t")
        outfile.write(a)
outfile.close()    