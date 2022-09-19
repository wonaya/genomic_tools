import os,sys
"""
## merge between CpG, CHG
genos = ['Mo17','Oh43', 'CML322', 'Tx303']
for geno in genos :
    print geno
    for chr in range(1,11) :
        outfile = open("test_dmr_cpg_chg_merge_"+str(geno)+"_"+str(chr)+".txt", 'w')
        print chr
        cpg_loc = []
        for a in open("temp_alldmr_"+str(chr)+"_CpG.txt", 'r') :
            if a.split("\t")[1] == geno :
                cpg_loc.append(range(int(a.split("\t")[2]), int(a.split("\t")[3].strip("\n"))+1))
        chg_loc = []
        for a in open("temp_alldmr_"+str(chr)+"_CHG.txt", 'r') :
            if a.split("\t")[1] == geno :
                chg_loc.append(range(int(a.split("\t")[2]), int(a.split("\t")[3].strip("\n"))+1))
        overlap_loc = []
        for chg in chg_loc :
            for cpg in cpg_loc :
                if len(set(chg)&set(cpg)) > 0 :
                    #print min(min(chg), max(chg), min(cpg), max(cpg)), max(min(chg), max(chg), min(cpg), max(cpg))
                    overlap_loc.append(range(min(min(chg), max(chg), min(cpg), max(cpg)),max(min(chg), max(chg), min(cpg), max(cpg))+1))
        #print len(cpg_loc), len(chg_loc), len(overlap_loc)
        
        for cpg in cpg_loc :
            count = 0
            for overlap in overlap_loc :
                if len(set(cpg)&set(overlap)) > 0 :
                    count += 1
            #print min(min(cpg), max(cpg)), max(min(cpg), max(cpg)), count
            if count > 0 :
                outfile.write("B73\t"+str(geno)+"\t"+str(chr)+"\t"+str(min(overlap))+"\t"+str(max(overlap))+"\tCpG CHG\n")
            else : 
                outfile.write("B73\t"+str(geno)+"\t"+str(chr)+"\t"+str(min(cpg))+"\t"+str(max(cpg))+"\tCpG\n")
        for chg in chg_loc :
            count = 0
            for overlap in overlap_loc :
                if len(set(chg)&set(overlap)) > 0 :
                    count += 1
            #print min(min(chg), max(chg)), max(min(chg), max(chg)), count
            if count == 0 :
                outfile.write("B73\t"+str(geno)+"\t"+str(chr)+"\t"+str(min(chg))+"\t"+str(max(chg))+"\tCHG\n")
        outfile.close()
        os.system("sort -nk4 test_dmr_cpg_chg_merge_"+str(geno)+"_"+str(chr)+".txt > test_dmr_cpg_chg_merge_"+str(geno)+"_"+str(chr)+"_sorted.txt")
        os.system("rm -Rf test_dmr_cpg_chg_merge_"+str(geno)+"_"+str(chr)+".txt")
        os.system("cat test_dmr_cpg_chg_merge_"+str(geno)+"_"+str(chr)+"_sorted.txt >> test_dmr_cpg_chg_merge_"+str(geno)+"_sorted.txt")
        os.system("rm -Rf test_dmr_cpg_chg_merge_"+str(geno)+"_"+str(chr)+"_sorted.txt")

## merge between genos
sys.exit()
"""


genos = ['Mo17','Oh43', 'CML322', 'Tx303']
for chr in range(1,11) :
    print chr
    geno1 = {}
    for a in open("test_dmr_cpg_chg_merge_Mo17_sorted.txt", 'r') :
        if int(a.split("\t")[2]) == chr :
            geno1[str(a.split("\t")[3])+"-"+str(a.split("\t")[4])] = ['Mo17']
    
    for a in open("test_dmr_cpg_chg_merge_Oh43_sorted.txt", 'r') : 
        if int(a.split("\t")[2]) == chr :
            count = 0
            overlap_coord = []
            for gen1 in geno1.keys(): 
                if len(set(range(int(gen1.split("-")[0]),int(gen1.split("-")[1])+1))&set(range(int(a.split("\t")[3]), int(a.split("\t")[4])+1))) > 0 : 
                    count += 1
                    overlap_coord.append(str(gen1.split("-")[0])+"-"+str(gen1.split("-")[1]))
            coords = []
            if count > 0 :
                coords = []
                coords.append(int(a.split("\t")[3]))
                coords.append(int(a.split("\t")[4]))
                for overlap in overlap_coord :
                    coords.append(int(overlap.split("-")[0]))
                    coords.append(int(overlap.split("-")[1]))
                    del geno1[overlap]
                geno1[str(min(coords))+"-"+str(max(coords))] = ['Mo17', 'Oh43']
            else :
                geno1[a.split("\t")[3]+"-"+a.split("\t")[4]] = ['Oh43']
        
    for a in open("test_dmr_cpg_chg_merge_CML322_sorted.txt", 'r') : 
        if int(a.split("\t")[2]) == chr :
            count = 0
            overlap_coord = []
            overlap_geno = []
            for gen1 in geno1.keys(): 
                if len(set(range(int(gen1.split("-")[0]),int(gen1.split("-")[1])+1))&set(range(int(a.split("\t")[3]), int(a.split("\t")[4])+1))) > 0 : 
                    count += 1
                    overlap_coord.append(str(gen1.split("-")[0])+"-"+str(gen1.split("-")[1]))
                    for genos1 in geno1[gen1] :
                        overlap_geno.append(genos1)
            coords = []
            if count > 0 :
                coords = []
                coords.append(int(a.split("\t")[3]))
                coords.append(int(a.split("\t")[4]))
                for overlap in overlap_coord :
                    coords.append(int(overlap.split("-")[0]))
                    coords.append(int(overlap.split("-")[1]))
                    del geno1[overlap]
                
                overlap_geno.append('CML322')
                geno1[str(min(coords))+"-"+str(max(coords))]  = overlap_geno
            else :
                geno1[a.split("\t")[3]+"-"+a.split("\t")[4]] = ['CML322']
        
    for a in open("test_dmr_cpg_chg_merge_Tx303_sorted.txt", 'r') : 
        if int(a.split("\t")[2]) == chr :
            count = 0
            overlap_coord = []
            overlap_geno = []
            for gen1 in geno1.keys(): 
                 if len(set(range(int(gen1.split("-")[0]),int(gen1.split("-")[1])+1))&set(range(int(a.split("\t")[3]), int(a.split("\t")[4])+1))) > 0 : 
                    count += 1
                    overlap_coord.append(str(gen1.split("-")[0])+"-"+str(gen1.split("-")[1]))
                    for genos1 in geno1[gen1] :
                        overlap_geno.append(genos1)
            coords = []
            if count > 0 :
                coords = []
                coords.append(int(a.split("\t")[3]))
                coords.append(int(a.split("\t")[4]))
                for overlap in overlap_coord :
                    coords.append(int(overlap.split("-")[0]))
                    coords.append(int(overlap.split("-")[1]))
                    del geno1[overlap]
                
                overlap_geno.append('Tx303')
                geno1[str(min(coords))+"-"+str(max(coords))]  = overlap_geno
                    
            else :
                geno1[a.split("\t")[3]+"-"+a.split("\t")[4]] = ['Tx303']
        
    print len(geno1)
    dmr_total = {}
    outfile = open("test_dmr_cpg_chg_merge_genos_merge_"+str(chr)+".txt", 'w')
    for dmr in geno1.keys() :
        if len(geno1[dmr]) == 1 :
            for a in open("test_dmr_cpg_chg_merge_"+str(geno1[dmr][0])+"_sorted.txt", 'r') :
                if int(a.split("\t")[2]) == chr and a.split("\t")[3] == dmr.split("-")[0] :
                    dmr_total[dmr] = [[geno1[dmr]],[a.split("\t")[5].strip("\n")]]
        elif len(geno1[dmr]) > 1 :
            context =[]
            for gen1 in geno1[dmr] :
                for a in open("test_dmr_cpg_chg_merge_"+str(gen1)+"_sorted.txt", 'r') :
                    if int(a.split("\t")[2]) == chr and len(set(range(int(dmr.split("-")[0]),int(dmr.split("-")[1])+1))&set(range(int(a.split("\t")[3]),int(a.split("\t")[4])+1))) > 0 : 
                        context.append(a.split("\t")[5].strip("\n"))
            dmr_total[dmr] = [[geno1[dmr]],context]
        #print chr, dmr, geno1[dmr], ",".join(set(dmr_total[dmr][0][0])), ",".join(set(dmr_total[dmr][1]))
        #outfile.write(str(chr)+"\t"+str(dmr.split("-")[0])+"\t"+str(dmr.split("-")[1])+"\t"+",".join(set(dmr_total[dmr][0][0]))+"\t"+",".join(set(dmr_total[dmr][1]))+"\n")
        outfile.write(str(chr)+"\t"+str(dmr.split("-")[0])+"\t"+str(dmr.split("-")[1])+"\t"+"1\n")
    outfile.close()
    os.system("sort -nk2 test_dmr_cpg_chg_merge_genos_merge_"+str(chr)+".txt > test_dmr_cpg_chg_merge_genos_merge_"+str(chr)+"_sorted.txt")
    os.system("rm -Rf test_dmr_cpg_chg_merge_genos_merge_"+str(chr)+".txt")
    os.system("cat test_dmr_cpg_chg_merge_genos_merge_"+str(chr)+"_sorted.txt >> test_dmr_cpg_chg_merge_genos_merge.txt")
    os.system("rm -Rf test_dmr_cpg_chg_merge_genos_merge_"+str(chr)+"_sorted.txt")
    print "done"