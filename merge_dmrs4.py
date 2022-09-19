### use largemem
### use single command line instead of multiple
### check DMR files before hand

import os,sys
import numpy
import multiprocessing
from datetime import datetime
from types import *

def check_dmr_file(geno, chr, context) :
    dmr_file = open("B73_"+str(geno)+"_"+str(context)+"_dmr_chr"+str(chr)+".txt")
    dmr_lines = dmr_file.readlines()
    for dmr_line in dmr_lines[1:] :
        try:
            int(dmr_line.split("\t")[1])
            return "NA"
        except ValueError:
            return "B73_"+str(geno)+"_"+str(context)+"_dmr_chr"+str(chr)+".txt"
    
def merge_dmr(chr, context) :
    outfile = open("temp_alldmr_"+str(chr)+"_"+str(context)+".txt" , 'w')
    dmr_loc1 = []
    if os.path.isfile("B73_Mo17_"+str(context)+"_dmr_chr"+str(chr)+"_50bp_4DMC.txt") :
        for line1 in open("B73_Mo17_"+str(context)+"_dmr_chr"+str(chr)+"_50bp_4DMC.txt", 'r') :
            if line1.split("\t")[0] != 'chr' and float(line1.split("\t")[9].strip("\n")) <= 0.05:
                dmr_loc1.append(str(line1.split("\t")[1])+"-"+str(int(line1.split("\t")[2])+1))
                outfile.write("B73\tMo17\t"+str(line1.split("\t")[1])+"\t"+str(int(line1.split("\t")[2]))+"\n")
    dmr_loc2 = []                    
    if os.path.isfile("B73_Oh43_"+str(context)+"_dmr_chr"+str(chr)+"_50bp_4DMC.txt") :
        for line2 in open("B73_Oh43_"+str(context)+"_dmr_chr"+str(chr)+"_50bp_4DMC.txt", 'r') :
            if line2.split("\t")[0] != 'chr' and float(line2.split("\t")[9].strip("\n")) <= 0.05:
                dmr_loc2.append(str(line2.split("\t")[1])+"-"+str(int(line2.split("\t")[2])+1))
                outfile.write("B73\tOh43\t"+str(line2.split("\t")[1])+"\t"+str(int(line2.split("\t")[2]))+"\n")
    new_dmr = []
    no_overlap1 = []    
    no_overlap2 = []  
    with_overlap1 = []
    with_overlap2 = []  
    print "chr", chr, "before merging dmrs B73_Mo17,Oh43_"+str(context), len(dmr_loc1), len(dmr_loc2)   
    for dmr1 in dmr_loc1 :
        for dmr2 in dmr_loc2 :
            if dmr1 != dmr2 :
                if len(set(range(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]))) & set(range(int(dmr2.split("-")[0]), int(dmr2.split("-")[1])))) > 0 :
                    if str(min(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]), int(dmr2.split("-")[0]), int(dmr2.split("-")[1])))+"-"+str(max(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]), int(dmr2.split("-")[0]), int(dmr2.split("-")[1]))) not in new_dmr :
                        new_dmr.append(str(min(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]), int(dmr2.split("-")[0]), int(dmr2.split("-")[1])))+"-"+str(max(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]), int(dmr2.split("-")[0]), int(dmr2.split("-")[1]))))
                        with_overlap1.append(dmr1)
                        with_overlap2.append(dmr2)
                else :
                    if dmr1 not in no_overlap1 :
                        no_overlap1.append(dmr1)
    for dmr2 in dmr_loc2 :
        for dmr1 in dmr_loc1 :
            if dmr1 != dmr2 :
                if len(set(range(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]))) & set(range(int(dmr2.split("-")[0]), int(dmr2.split("-")[1])))) == 0 :
                    if dmr2 not in no_overlap2 :
                        no_overlap2.append(dmr2)
    for dmr1 in dmr_loc1 :
        if dmr1 not in with_overlap1 :
            new_dmr.append(dmr1)
    for dmr2 in dmr_loc2 :
        if dmr2 not in with_overlap2 :
            new_dmr.append(dmr2)
    print "chr", chr, "individual dmrs B73_Mo17,Oh43_"+str(context), len(dmr_loc1), len(dmr_loc2)   
    dmr_loc1 = new_dmr
    dmr_loc2 = []
    if os.path.isfile("B73_CML322_"+str(context)+"_dmr_chr"+str(chr)+"_50bp_4DMC.txt") :
        for line3 in open("B73_CML322_"+str(context)+"_dmr_chr"+str(chr)+"_50bp_4DMC.txt", 'r') :
            if line3.split("\t")[0] != 'chr' and float(line3.split("\t")[9].strip("\n")) <= 0.05:
                dmr_loc2.append(str(line3.split("\t")[1])+"-"+str(int(line3.split("\t")[2])+1))
                outfile.write("B73\tCML322\t"+str(line3.split("\t")[1])+"\t"+str(int(line3.split("\t")[2]))+"\n")
    print "chr", chr, "individual dmrs B73_Mo17,Oh43,CML322_"+str(context), len(dmr_loc1), len(dmr_loc2)
    
    new_dmr = []
    no_overlap1 = []    
    no_overlap2 = []  
    with_overlap1 = []
    with_overlap2 = []  
    for dmr1 in dmr_loc1 :
        for dmr2 in dmr_loc2 :
            if dmr1 != dmr2 :
                if len(set(range(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]))) & set(range(int(dmr2.split("-")[0]), int(dmr2.split("-")[1])))) > 0 :
                    if str(min(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]), int(dmr2.split("-")[0]), int(dmr2.split("-")[1])))+"-"+str(max(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]), int(dmr2.split("-")[0]), int(dmr2.split("-")[1]))) not in new_dmr :
                        new_dmr.append(str(min(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]), int(dmr2.split("-")[0]), int(dmr2.split("-")[1])))+"-"+str(max(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]), int(dmr2.split("-")[0]), int(dmr2.split("-")[1]))))
                        with_overlap1.append(dmr1)
                        with_overlap2.append(dmr2)
                else :
                    if dmr1 not in no_overlap1 :
                        no_overlap1.append(dmr1)
    for dmr2 in dmr_loc2 :
        for dmr1 in dmr_loc1 :
            if dmr1 != dmr2 :
                if len(set(range(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]))) & set(range(int(dmr2.split("-")[0]), int(dmr2.split("-")[1])))) == 0 :
                    if dmr2 not in no_overlap2 :
                        no_overlap2.append(dmr2)
    for dmr1 in dmr_loc1 :
        if dmr1 not in with_overlap1 :
            new_dmr.append(dmr1)
    for dmr2 in dmr_loc2 :
        if dmr2 not in with_overlap2 :
            new_dmr.append(dmr2)
    dmr_loc1 = new_dmr
    
    dmr_loc2 = []
    if os.path.isfile("B73_Tx303_"+str(context)+"_dmr_chr"+str(chr)+"_50bp_4DMC.txt") :
        for line1 in open("B73_Tx303_"+str(context)+"_dmr_chr"+str(chr)+"_50bp_4DMC.txt", 'r') :
            if line1.split("\t")[0] != 'chr' and float(line1.split("\t")[9].strip("\n")) <= 0.05:
                dmr_loc2.append(str(line1.split("\t")[1])+"-"+str(int(line1.split("\t")[2])+1))
                outfile.write("B73\tTx303\t"+str(line1.split("\t")[1])+"\t"+str(int(line1.split("\t")[2]))+"\n")
    
    print "chr", chr, "individual B73_Mo17,Oh43,CML322,Tx303_"+str(context), len(dmr_loc1), len(dmr_loc2)
    
    new_dmr = []
    no_overlap1 = []    
    no_overlap2 = []  
    with_overlap1 = []
    with_overlap2 = []  
    for dmr1 in dmr_loc1 :
        for dmr2 in dmr_loc2 :
            if dmr1 != dmr2 :
                if len(set(range(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]))) & set(range(int(dmr2.split("-")[0]), int(dmr2.split("-")[1])))) > 0 :
                    if str(min(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]), int(dmr2.split("-")[0]), int(dmr2.split("-")[1])))+"-"+str(max(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]), int(dmr2.split("-")[0]), int(dmr2.split("-")[1]))) not in new_dmr :
                        new_dmr.append(str(min(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]), int(dmr2.split("-")[0]), int(dmr2.split("-")[1])))+"-"+str(max(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]), int(dmr2.split("-")[0]), int(dmr2.split("-")[1]))))
                        with_overlap1.append(dmr1)
                        with_overlap2.append(dmr2)
                else :
                    if dmr1 not in no_overlap1 :
                        no_overlap1.append(dmr1)
    for dmr2 in dmr_loc2 :
        for dmr1 in dmr_loc1 :
            if dmr1 != dmr2 :
                if len(set(range(int(dmr1.split("-")[0]), int(dmr1.split("-")[1]))) & set(range(int(dmr2.split("-")[0]), int(dmr2.split("-")[1])))) == 0 :
                    if dmr2 not in no_overlap2 :
                        no_overlap2.append(dmr2)
    for dmr1 in dmr_loc1 :
        if dmr1 not in with_overlap1 :
            new_dmr.append(dmr1)
    for dmr2 in dmr_loc2 :
        if dmr2 not in with_overlap2 :
            new_dmr.append(dmr2)
    outfile.close()
    
    final_dmr = new_dmr
    
    start =[]
    for dmr in final_dmr :
        start.append(int(dmr.split("-")[0]))
    start.sort()
    
    sorted_dmr = []
    dmr_dict = {}
    for point in start :
        for dmr in final_dmr :
            if int(dmr.split("-")[0]) == point :
                if dmr not in sorted_dmr :
                    sorted_dmr.append(dmr)
                    dmr_dict[int(dmr.split("-")[0])] = int(dmr.split("-")[1])
    
    print "final no. of DMRs B73_Mo17_Oh43_CML322_Tx303", str(context), chr, len(sorted_dmr), len(dmr_dict)
    
    print "###preparing bed file for calculating reads mapped for DMR regions"
    ### get coverage
    outfile = open("temp_"+str(chr)+"_"+str(context)+".bed", 'w')
    for dmr in sorted_dmr :
        outfile.write("chr"+str(chr)+"\t"+str(dmr.split("-")[0])+"\t"+str(dmr.split("-")[1])+"\t0\n")
    outfile.close()

#def check_pos(chr, context) :
def check_pos(chr) :
    print "###calculating total number of possible methylation sites"
    #outfile = open("temp_"+str(chr)+"_"+str(context)+".allC", 'w')
    outfile = open("temp_"+str(chr)+"_merged.allC", 'w')
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
    chh_list = ["CAA", "CAC", "CAT", "CCA", "CCC", "CCT", "CTA", "CTC", "CTT"]
    chh_rev_list = ["TTG", "GTG", "ATG", "TGG", "GGG", "AGG", "TAG", "GAG", "AAG"]
    chg_list = ["CAG", "CTG", "CCG"]
    chg_rev_list = ["CTG", "CAG", "CGG"]
    cg_file = open("temp_cpg_"+str(chr)+"_merged.txt", 'w')
    chg_file = open("temp_chg_"+str(chr)+"_merged.txt", 'w')
    chh_file = open("temp_chh_"+str(chr)+"_merged.txt", 'w')
    #cg_file = open("temp_cpg_"+str(chr)+"_"+str(context)+".txt", 'w')
    #chg_file = open("temp_chg_"+str(chr)+"_"+str(context)+".txt", 'w')
    #chh_file = open("temp_chh_"+str(chr)+"_"+str(context)+".txt", 'w')
    #for a in open("temp_"+str(chr)+"_"+str(context)+".bed", 'r') :
    #for a in open("test_dmr_cpg_chg_merge_genos_merge.bedGraph", 'r'):
    a_s = []
    afile = open("test_dmr_cpg_chg_merge_genos_merge.bedGraph", 'r')
    alines = afile.readlines()
    for aline in alines :
        if aline.split("\t")[0] == str(chr) :
            a_s.append(aline)
    for a in a_s :
        cg_loc_list = []
        chg_loc_list = []
        chh_loc_list = []
        loc = 0
        for base in whole_list[int(chr)-1][int(a.split("\t")[1])-1:int(a.split("\t")[2])+1] : 
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
        cg_loc_list.sort()
        chg_loc_list.sort()
        chh_loc_list.sort()
        cg_file.write("\t".join(a.split("\t")[:3]))
        cg_file.write("\t")
        cg_file.write(",".join(str(e) for e in cg_loc_list)) 
        cg_file.write("\n")
        chg_file.write("\t".join(a.split("\t")[:3]))
        chg_file.write("\t")
        chg_file.write(",".join(str(e) for e in chg_loc_list)) 
        chg_file.write("\n")
        chh_file.write("\t".join(a.split("\t")[:3]))
        chh_file.write("\t")
        chh_file.write(",".join(str(e) for e in chh_loc_list)) 
        chh_file.write("\n")
        outfile.write("\t".join(a.split("\t")[:3])+"\t"+str(len(cg_loc_list))+"\t"+str(len(chg_loc_list))+"\t"+str(len(chh_loc_list))+"\n")
    cg_file.close()
    chg_file.close()
    chh_file.close()
    outfile.close()

def coverage(chr):
#def coverage(chr, context):
    ### read coverage calculation
    outfile = open("test_dmr_cpg_chg_merge_genos_merge_"+str(chr)+".bed", 'w')
    for a in open("test_dmr_cpg_chg_merge_genos_merge.bedGraph",'r'):
        if int(a.split("\t")[0]) == chr :
            outfile.write("chr"+str(a.split("\t")[0]))
            outfile.write("\t")
            outfile.write("\t".join(a.split("\t")[1:3]))
            outfile.write("\t0\n")
    outfile.close()
    #os.system("/opt/apps/bedtools/2.19.0/bin/multiBamCov -bams B73_all3/B73_all3_R1_bt202_sorted.bam Mo17_all3/Mo17_all3_bt202_sorted.bam Oh43_all3/Oh43_all3_bt202_sorted.bam CML322_all3/CML322_all3_bt202_sorted.bam Tx303_all3/Tx303_all3_bt202_sorted.bam -bed temp_"+str(chr)+"_"+str(context)+".bed > temp_bam_count_"+str(chr)+"_"+str(context)+".txt")
    os.system("/opt/apps/bedtools/2.19.0/bin/multiBamCov -bams ../10-01-13_5genos/B73_all3/B73_all3_R1_bt202_sorted.bam ../10-01-13_5genos/Mo17_all3/Mo17_all3_bt202_sorted.bam ../10-01-13_5genos/Oh43_all3/Oh43_all3_bt202_sorted.bam ../10-01-13_5genos/CML322_all3/CML322_all3_bt202_sorted.bam ../10-01-13_5genos/Tx303_all3/Tx303_all3_bt202_sorted.bam -bed test_dmr_cpg_chg_merge_genos_merge_"+str(chr)+".bed > temp_bam_count_"+str(chr)+"_merged.txt")
    print "###calculating length of DMR"
    ### calculate length
    outfile = open("temp_"+str(chr)+"_merged.length", 'w')
    for a in open("test_dmr_cpg_chg_merge_genos_merge_"+str(chr)+".bed", 'r') :
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(int(a.split("\t")[2])-int(a.split("\t")[1]))+"\n")
    outfile.close()

def cpg_metrate(chr):
    ### calculate non-covered C
    ### B73 % methylation
    print "###getting B73 methylation rates for individual DMRs", datetime.now()
    cg_loc = open("temp_cpg_"+str(chr)+"_merged.txt" ,'r') 
    #cg_loc = open("temp_cpg_"+str(chr)+"_"+str(context)+".txt" ,'r') 
    cg_loc_lines = cg_loc.readlines()
    valid_cg = []
    for cg_loc_line in cg_loc_lines :
        if len(cg_loc_line.split("\t")[3].strip("\n").split(",")) > 0 and len(cg_loc_line.split("\t")[3].strip("\n").split(",")[0]) != 0:
            for x in cg_loc_line.split("\t")[3].strip("\n").split(",") :
                valid_cg.append(int(x))
    valid_cg.sort()
    print "read B73", datetime.now()
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    for a in open("../07-16-14_bismar12.3/CpG_B73_all3_R1_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
    print "done", datetime.now()
    
    outfile = open("temp_cpg_b73_"+str(chr)+"_merged.cov2", 'w')
    for a in cg_loc_lines :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for x in set(range(int(a.split("\t")[1]),(int(a.split("\t")[2])+1)))&set(valid_cg) :
            countmet += large_list[x][0]
            countunmet += large_list[x][1]
            uniquec += large_list[x][2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "Mo17", datetime.now()
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    for a in open("../07-16-14_bismar12.3/CpG_Mo17_all3_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
    print "done", datetime.now()
    outfile = open("temp_cpg_mo17_"+str(chr)+"_merged.cov2", 'w')
    for a in cg_loc_lines :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for x in set(range(int(a.split("\t")[1]),(int(a.split("\t")[2])+1)))&set(valid_cg) :
            countmet += large_list[x][0]
            countunmet += large_list[x][1]
            uniquec += large_list[x][2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "Oh43", datetime.now()
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    for a in open("../07-16-14_bismar12.3/CpG_Oh43_all3_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
    print "done", datetime.now()
    outfile = open("temp_cpg_oh43_"+str(chr)+"_merged.cov2", 'w')
    for a in cg_loc_lines :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for x in set(range(int(a.split("\t")[1]),(int(a.split("\t")[2])+1)))&set(valid_cg) :
            countmet += large_list[x][0]
            countunmet += large_list[x][1]
            uniquec += large_list[x][2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "CML322", datetime.now()
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    for a in open("../07-16-14_bismar12.3/CpG_CML322_all3_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
    print "done", datetime.now()
    outfile = open("temp_cpg_cml322_"+str(chr)+"_merged.cov2", 'w')
    for a in cg_loc_lines :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for x in set(range(int(a.split("\t")[1]),(int(a.split("\t")[2])+1)))&set(valid_cg) :
            countmet += large_list[x][0]
            countunmet += large_list[x][1]
            uniquec += large_list[x][2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "Tx303", datetime.now()
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    for a in open("../07-16-14_bismar12.3/CpG_Tx303_all3_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
    print "done", datetime.now()
    outfile = open("temp_cpg_tx303_"+str(chr)+"_merged.cov2", 'w')
    for a in cg_loc_lines :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for x in set(range(int(a.split("\t")[1]),(int(a.split("\t")[2])+1)))&set(valid_cg) :
            countmet += large_list[x][0]
            countunmet += large_list[x][1]
            uniquec += large_list[x][2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
def chg_metrate(chr):    
    print "###getting B73 methylation rates for individual DMRs", datetime.now()
    chg_loc = open("temp_chg_"+str(chr)+"_merged.txt" ,'r') 
    chg_loc_lines = chg_loc.readlines()
    valid_chg = []
    for chg_loc_line in chg_loc_lines :
        for x in chg_loc_line.split("\t")[3].strip("\n").split(",") :
            if len(x) > 1 :
                valid_chg.append(int(x))
    valid_chg.sort()
    print len(valid_chg)
    
    ### B73 CHG % methylation
    print "### calculate CHG methylation rate", datetime.now()
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    for a in open("../07-16-14_bismar12.3/CHG_B73_all3_R1_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
    print "done", datetime.now()
    outfile = open("temp_chg_b73_"+str(chr)+"_merged.cov2", 'w')
    for a in chg_loc_lines :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for x in set(range(int(a.split("\t")[1]),(int(a.split("\t")[2])+1)))&set(valid_chg) :
            countmet += large_list[x][0]
            countunmet += large_list[x][1]
            uniquec += large_list[x][2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "Mo17"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    for a in open("../07-16-14_bismar12.3/CHG_Mo17_all3_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
    print "done", datetime.now()
    outfile = open("temp_chg_mo17_"+str(chr)+"_merged.cov2", 'w')
    for a in chg_loc_lines :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for x in set(range(int(a.split("\t")[1]),(int(a.split("\t")[2])+1)))&set(valid_chg) :
            countmet += large_list[x][0]
            countunmet += large_list[x][1]
            uniquec += large_list[x][2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "Oh43"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    for a in open("../07-16-14_bismar12.3/CHG_Oh43_all3_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
    print "done", datetime.now()
    outfile = open("temp_chg_oh43_"+str(chr)+"_merged.cov2", 'w')
    for a in chg_loc_lines :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for x in set(range(int(a.split("\t")[1]),(int(a.split("\t")[2])+1)))&set(valid_chg) :
            countmet += large_list[x][0]
            countunmet += large_list[x][1]
            uniquec += large_list[x][2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "CML322"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    for a in open("../07-16-14_bismar12.3/CHG_CML322_all3_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
    print "done", datetime.now()
    outfile = open("temp_chg_cml322_"+str(chr)+"_merged.cov2", 'w')
    for a in chg_loc_lines :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for x in set(range(int(a.split("\t")[1]),(int(a.split("\t")[2])+1)))&set(valid_chg) :
            countmet += large_list[x][0]
            countunmet += large_list[x][1]
            uniquec += large_list[x][2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list    

    print "Tx303"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    for a in open("../07-16-14_bismar12.3/CHG_Tx303_all3_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
    print "done", datetime.now()
    outfile = open("temp_chg_tx303_"+str(chr)+"_merged.cov2", 'w')
    for a in chg_loc_lines :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for x in set(range(int(a.split("\t")[1]),(int(a.split("\t")[2])+1)))&set(valid_chg) :
            countmet += large_list[x][0]
            countunmet += large_list[x][1]
            uniquec += large_list[x][2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list

def chh_metrate(chr):     
    print "###getting B73 methylation rates for individual DMRs", datetime.now()
    chh_loc = open("temp_chh_"+str(chr)+"_merged.txt" ,'r') 
    chh_loc_lines = chh_loc.readlines()
    valid_chh = []
    for chh_loc_line in chh_loc_lines :
        for x in chh_loc_line.split("\t")[3].strip("\n").split(",") :
            if len(x) > 1 :
                valid_chh.append(int(x))
    valid_chh.sort()
    print len(valid_chh)
    
    
    ### B73 CHH % methylation
    print "### calculate CHH methylation rate"
    print "B73"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    for a in open("../07-16-14_bismar12.3/CHH_B73_all3_R1_bt202_sortedn_"+str(chr)+".cov", 'r') :
        if a.split("\t")[0] == "chr"+str(chr) : 
            large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
            large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
            large_list[int(a.split("\t")[1])][2] += 1
    print "done", datetime.now()
    outfile = open("temp_chh_b73_"+str(chr)+"_merged.cov2", 'w')
    for a in chh_loc_lines :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for x in set(range(int(a.split("\t")[1]),(int(a.split("\t")[2])+1)))&set(valid_chh) :
            countmet += large_list[x][0]
            countunmet += large_list[x][1]
            uniquec += large_list[x][2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    """
    print "Mo17"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    print "making list done"
    for a in open("../07-16-14_bismar12.3/CHH_Mo17_all3_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
    print "done", datetime.now()
    outfile = open("temp_chh_mo17_"+str(chr)+"_merged.cov2", 'w')
    for a in chh_loc_lines :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for x in set(range(int(a.split("\t")[1]),(int(a.split("\t")[2])+1)))&set(valid_chh) :
            countmet += large_list[x][0]
            countunmet += large_list[x][1]
            uniquec += large_list[x][2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "Oh43"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    for a in open("../07-16-14_bismar12.3/CHH_Oh43_all3_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
    print "done", datetime.now()
    outfile = open("temp_chh_oh43_"+str(chr)+"_merged.cov2", 'w')
    for a in chh_loc_lines :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for x in set(range(int(a.split("\t")[1]),(int(a.split("\t")[2])+1)))&set(valid_chh) :
            countmet += large_list[x][0]
            countunmet += large_list[x][1]
            uniquec += large_list[x][2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "CML322"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    for a in open("../07-16-14_bismar12.3/CHH_CML322_all3_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
    print "done", datetime.now()
    outfile = open("temp_chh_cml322_"+str(chr)+"_merged.cov2", 'w')
    for a in chh_loc_lines :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for x in set(range(int(a.split("\t")[1]),(int(a.split("\t")[2])+1)))&set(valid_chh) :
            countmet += large_list[x][0]
            countunmet += large_list[x][1]
            uniquec += large_list[x][2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "Tx303"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    for a in open("../07-16-14_bismar12.3/CHH_Tx303_all3_bt202_sortedn_"+str(chr)+".cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
    print "done", datetime.now()
    outfile = open("temp_chh_tx303_"+str(chr)+"_merged.cov2", 'w')
    for a in chh_loc_lines :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for x in set(range(int(a.split("\t")[1]),(int(a.split("\t")[2])+1)))&set(valid_chh) :
            countmet += large_list[x][0]
            countunmet += large_list[x][1]
            uniquec += large_list[x][2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    """
def medip(chr):
    print "###start meDIP analysis"
    ### meDIP
    outfile = open("temp_meDIP_"+str(chr)+"_merged.txt", 'w')
    for a in open("test_dmr_cpg_chg_merge_genos_merge_"+str(chr)+".bed", 'r') :
        probe_name = []
        b73 = []
        mo17 = []
        oh43 = []
        cml322 = []
        tx303 = []
        for b in open("5geno_fromDiDIP_meDIP_normalized_values_30-1-14.txt" ,'r' ) :
            if b.split("\t")[0] != "chromosome" and b.split("\t")[0] == "chr"+str(chr) :
                if int(b.split("\t")[1]) in range(int(a.split("\t")[1])-300,int(a.split("\t")[2])+301) or int(b.split("\t")[2]) in range(int(a.split("\t")[1])-300,int(a.split("\t")[2])+301) :
                    probe_name.append(b.split("\t")[3])
                    b73.append(float(b.split("\t")[4]))
                    mo17.append(float(b.split("\t")[5]))
                    oh43.append(float(b.split("\t")[6]))
                    cml322.append(float(b.split("\t")[7]))
                    tx303.append(float(b.split("\t")[8].strip("\n")))
        if len(probe_name) == 0 :
            outfile.write("\t".join(a.split("\t")[:3])+"\tNA\tNA\tNA\tNA\tNA\tNA\n")
        else :
            outfile.write("\t".join(a.split("\t")[:3])+"\t"+",".join(probe_name)+"\t"+str(numpy.mean(b73))+"\t"+str(numpy.mean(mo17))+"\t"+str(numpy.mean(oh43))+"\t"+str(numpy.mean(cml322))+"\t"+str(numpy.mean(tx303))+"\n")
    outfile.close()

def medip_tmp(chr, context, tmp):
    if not os.path.isfile("temp_"+str(chr)+"_"+str(context)+".bed_"+str(tmp)) :
        sys.exit()
    print "###start meDIP analysis", chr, context, tmp
    ### meDIP
    outfile = open("temp_meDIP_"+str(chr)+"_"+str(context)+".txt_"+str(tmp), 'w')
    for a in open("temp_"+str(chr)+"_"+str(context)+".bed_"+str(tmp), 'r') :
        probe_name = []
        b73 = []
        mo17 = []
        oh43 = []
        cml322 = []
        tx303 = []
        for b in open("5geno_fromDiDIP_meDIP_normalized_values_30-1-14.txt" ,'r' ) :
            if b.split("\t")[0] != "chromosome" and b.split("\t")[0] == "chr"+str(chr) :
                if int(b.split("\t")[1]) in range(int(a.split("\t")[1])-300,int(a.split("\t")[2])+301) or int(b.split("\t")[2]) in range(int(a.split("\t")[1])-300,int(a.split("\t")[2])+301) :
                    probe_name.append(b.split("\t")[3])
                    b73.append(float(b.split("\t")[4]))
                    mo17.append(float(b.split("\t")[5]))
                    oh43.append(float(b.split("\t")[6]))
                    cml322.append(float(b.split("\t")[7]))
                    tx303.append(float(b.split("\t")[8].strip("\n")))
        if len(probe_name) == 0 :
            outfile.write("\t".join(a.split("\t")[:3])+"\tNA\tNA\tNA\tNA\tNA\tNA\n")
        else :
            outfile.write("\t".join(a.split("\t")[:3])+"\t"+",".join(probe_name)+"\t"+str(numpy.mean(b73))+"\t"+str(numpy.mean(mo17))+"\t"+str(numpy.mean(oh43))+"\t"+str(numpy.mean(cml322))+"\t"+str(numpy.mean(tx303))+"\n")
    outfile.close()
   
    
def combine(chr) :    
    outfile = open("Merged_DMR_5genos_chr"+str(chr)+"_08-08-14.txt", 'w')
    a = open("test_dmr_cpg_chg_merge_genos_merge_"+str(chr)+".bed", 'r')
    b1 = open("temp_cpg_b73_"+str(chr)+"_merged.cov2", 'r')
    b2 = open("temp_cpg_mo17_"+str(chr)+"_merged.cov2", 'r')
    b3 = open("temp_cpg_oh43_"+str(chr)+"_merged.cov2", 'r')
    b4 = open("temp_cpg_cml322_"+str(chr)+"_merged.cov2", 'r')
    b5 = open("temp_cpg_tx303_"+str(chr)+"_merged.cov2", 'r')
    c = open("temp_"+str(chr)+"_merged.allC", 'r')
    d = open("temp_bam_count_"+str(chr)+"_merged.txt", 'r')
    e = open("temp_meDIP_"+str(chr)+"_merged.txt", 'r')
    f1 = open("temp_chg_b73_"+str(chr)+"_merged.cov2", 'r')
    f2 = open("temp_chg_mo17_"+str(chr)+"_merged.cov2", 'r')
    f3 = open("temp_chg_oh43_"+str(chr)+"_merged.cov2", 'r')
    f4 = open("temp_chg_cml322_"+str(chr)+"_merged.cov2", 'r')
    f5 = open("temp_chg_tx303_"+str(chr)+"_merged.cov2", 'r')
    g1 = open("temp_chh_b73_"+str(chr)+"_merged.cov2", 'r')
    g2 = open("temp_chh_mo17_"+str(chr)+"_merged.cov2", 'r')
    g3 = open("temp_chh_oh43_"+str(chr)+"_merged.cov2", 'r')
    g4 = open("temp_chh_cml322_"+str(chr)+"_merged.cov2", 'r')
    g5 = open("temp_chh_tx303_"+str(chr)+"_merged.cov2", 'r')
    #h = open("temp_alldmr_"+str(chr)+"_merged.txt", 'r')
    h1 = open("test_dmr_cpg_chg_merge_Mo17_sorted.txt", 'r')
    h2 = open("test_dmr_cpg_chg_merge_Oh43_sorted.txt", 'r')
    h3 = open("test_dmr_cpg_chg_merge_CML322_sorted.txt", 'r')
    h4 = open("test_dmr_cpg_chg_merge_Tx303_sorted.txt", 'r')
    
    alines = a.readlines()
    b1lines_all = b1.readlines()
    b1lines = []
    for b1line in b1lines_all : 
        if b1line.split("\t")[0] == str(chr) :
            b1lines.append(b1line)
    b2lines_all = b2.readlines()
    b2lines = []
    for b2line in b2lines_all : 
        if b2line.split("\t")[0] == str(chr) :
            b2lines.append(b2line)
    b3lines_all = b3.readlines()
    b3lines = []
    for b3line in b3lines_all : 
        if b3line.split("\t")[0] == str(chr) :
            b3lines.append(b3line)
    b4lines_all = b4.readlines()
    b4lines = []
    for b4line in b4lines_all : 
        if b4line.split("\t")[0] == str(chr) :
            b4lines.append(b4line)
    b5lines_all = b5.readlines()
    b5lines = []
    for b5line in b5lines_all : 
        if b5line.split("\t")[0] == str(chr) :
            b5lines.append(b5line)
    clines_all = c.readlines()
    clines = []
    for cline in clines_all : 
        if cline.split("\t")[0] == str(chr) :
            clines.append(cline)
    dlines = d.readlines()
    elines = e.readlines()
    f1lines_all = f1.readlines()
    f1lines = []
    for f1line in f1lines_all : 
        if f1line.split("\t")[0] == str(chr) :
            f1lines.append(f1line)
    f2lines_all = f2.readlines()
    f2lines = []
    for f2line in f2lines_all : 
        if f2line.split("\t")[0] == str(chr) :
            f2lines.append(f2line)
    f3lines_all = f3.readlines()
    f3lines = []
    for f3line in f3lines_all : 
        if f3line.split("\t")[0] == str(chr) :
            f3lines.append(f3line)
    f4lines_all = f4.readlines()
    f4lines = []
    for f4line in f4lines_all : 
        if f4line.split("\t")[0] == str(chr) :
            f4lines.append(f4line)
    f5lines_all = f5.readlines()
    f5lines = []
    for f5line in f5lines_all : 
        if f5line.split("\t")[0] == str(chr) :
            f5lines.append(f5line)
    
    g1lines_all = g1.readlines()
    g1lines = []
    for g1line in g1lines_all : 
        if g1line.split("\t")[0] == str(chr) :
            g1lines.append(g1line)
    g2lines_all = g2.readlines()
    g2lines = []
    for g2line in g2lines_all : 
        if g2line.split("\t")[0] == str(chr) :
            g2lines.append(g2line)
    g3lines_all = g3.readlines()
    g3lines = []
    for g3line in g3lines_all : 
        if g3line.split("\t")[0] == str(chr) :
            g3lines.append(g3line)
    g4lines_all = g4.readlines()
    g4lines = []
    for g4line in g4lines_all : 
        if g4line.split("\t")[0] == str(chr) :
            g4lines.append(g4line)
    g5lines_all = g5.readlines()
    g5lines = []
    for g5line in g5lines_all : 
        if g5line.split("\t")[0] == str(chr) :
            g5lines.append(g5line)
    hlines = []
    h1lines = h1.readlines()
    for h1line in h1lines :
        if h1line.split("\t")[2] == str(chr) :
            hlines.append(h1line)
    h2lines = h2.readlines()
    for h2line in h2lines :
        if h2line.split("\t")[2] == str(chr) :
            hlines.append(h2line)
    h3lines = h3.readlines()
    for h3line in h3lines :
        if h3line.split("\t")[2] == str(chr) :
            hlines.append(h3line)
    h4lines = h4.readlines()
    for h4line in h4lines :
        if h4line.split("\t")[2] == str(chr) :
            hlines.append(h4line)
            
    
    print chr, len(alines), len(b1lines), len(b2lines), len(b3lines), len(b4lines), len(b5lines), len(clines), len(dlines), len(elines), len(f1lines), len(f2lines), len(f3lines), len(f4lines), len(f5lines), len(g1lines), len(g2lines), len(g3lines), len(g4lines), len(g5lines), len(hlines)
    
    outfile.write("chr\tstart\tend\tlength\tContributing genotypes\tContributing contexts\tB73 cov\tMo17 cov\tOh43 cov\tCML322 cov\tTx303 cov\tB73 CpG meth%\tMo17 CpG meth%\tOh43 CpG meth%\tCML322 CpG meth%\tTx303 CpG meth%\tB73 CHG meth%\tMo17 CHG meth%\tOh43 CHG meth%\tCML322 CHG meth%\tTx303 CHG meth%\tB73 CHH meth%\tMo17 CHH meth%\tOh43 CHH meth%\tCML322 CHH meth%\tTx303 CHH meth%\tB73 CG mapped\tMo17 CG mapped\tOh43 CG mapped\tCML322 CG mapped\tTx303 CG mapped\tB73 CHG mapped\tMo17 CHG mapped\tOh43 CHG mapped\tCML322 CHG mapped\tTx303 CHG mapped\tB73 CHH mapped\tMo17 CHH mapped\tOh43 CHH mapped\tCML322 CHH mapped\tTx303 CHH mapped\tNo. of probe +/-300bp\tB73 mean probe val\tMo17 mean probe val\tOh43 mean probe val\tCML322 mean probe val\tTx303 mean probe val\tTotal CpG\tTotal CHG\tTotal CHH\n")
    for aline in alines :
        coordinates = aline.split("\t")[:3]
        length = str(int(aline.split("\t")[2])-int(aline.split("\t")[1]))
        ## get dmr genos
        contrib_geno_list = []
        contrib_type_list = []
        for hline in hlines :
            if hline.split("\t")[2] == str(chr) :
                #list1 = range(int(hline.split("\t")[2]),int(hline.split("\t")[3].strip("\n"))+1)
                list1 = range(int(hline.split("\t")[3]),int(hline.split("\t")[4].strip("\n"))+1)
                list2 = range(int(aline.split("\t")[1]),int(aline.split("\t")[2])+1)
                if len(set(list1)&set(list2)) > 0 : 
                    if hline.split("\t")[1] not in contrib_geno_list :
                        contrib_geno_list.append(hline.split("\t")[1])
                    if hline.split("\t")[5].strip("\n") not in contrib_type_list :
                        contrib_type_list.append(hline.split("\t")[5].strip("\n"))
                        
        if len(contrib_geno_list) == 0 :
            contrib_geno = "None" ; sys.exit()
        else :
            contrib_geno = str(",".join(contrib_geno_list))
        if len(contrib_type_list) == 0 :
            contrib_type = "None" ; sys.exit()
        else :
            contrib_type = str(",".join(contrib_type_list))
        ## b73_1
        b73_uniquenoofC = str(b1lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(b1lines[alines.index(aline)].split("\t")[2])+float(b1lines[alines.index(aline)].split("\t")[3])) == 0 :
            b73_metrate = "NA"
        else : 
            b73_metrate  = str((float(b1lines[alines.index(aline)].split("\t")[2]))/(float(b1lines[alines.index(aline)].split("\t")[2])+float(b1lines[alines.index(aline)].split("\t")[3])))
        ## mo17_1
        mo17_uniquenoofC = str(b2lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(b2lines[alines.index(aline)].split("\t")[2])+float(b2lines[alines.index(aline)].split("\t")[3])) == 0 :
            mo17_metrate = "NA"
        else : 
            mo17_metrate  = str((float(b2lines[alines.index(aline)].split("\t")[2]))/(float(b2lines[alines.index(aline)].split("\t")[2])+float(b2lines[alines.index(aline)].split("\t")[3])))
        ## oh43_1
        oh43_uniquenoofC = str(b3lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(b3lines[alines.index(aline)].split("\t")[2])+float(b3lines[alines.index(aline)].split("\t")[3])) == 0 :
            oh43_metrate = "NA"
        else : 
            oh43_metrate  = str((float(b3lines[alines.index(aline)].split("\t")[2]))/(float(b3lines[alines.index(aline)].split("\t")[2])+float(b3lines[alines.index(aline)].split("\t")[3])))
        ## cml322_1
        cml322_uniquenoofC = str(b4lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(b4lines[alines.index(aline)].split("\t")[2])+float(b4lines[alines.index(aline)].split("\t")[3])) == 0 :
            cml322_metrate = "NA"
        else : 
            cml322_metrate  = str((float(b4lines[alines.index(aline)].split("\t")[2]))/(float(b4lines[alines.index(aline)].split("\t")[2])+float(b4lines[alines.index(aline)].split("\t")[3])))
        ## tx303_1
        tx303_uniquenoofC = str(b5lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(b5lines[alines.index(aline)].split("\t")[2])+float(b5lines[alines.index(aline)].split("\t")[3])) == 0 :
            tx303_metrate = "NA"
        else : 
            tx303_metrate  = str((float(b5lines[alines.index(aline)].split("\t")[2]))/(float(b5lines[alines.index(aline)].split("\t")[2])+float(b5lines[alines.index(aline)].split("\t")[3])))
        
        ## CHG
        ## b73_1
        b73_chg_uniquenoofC = str(f1lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(f1lines[alines.index(aline)].split("\t")[2])+float(f1lines[alines.index(aline)].split("\t")[3])) == 0 :
            b73_chg_metrate = "NA"
        else : 
            b73_chg_metrate  = str((float(f1lines[alines.index(aline)].split("\t")[2]))/(float(f1lines[alines.index(aline)].split("\t")[2])+float(f1lines[alines.index(aline)].split("\t")[3])))
        ## mo17_1
        mo17_chg_uniquenoofC = str(f2lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(f2lines[alines.index(aline)].split("\t")[2])+float(f2lines[alines.index(aline)].split("\t")[3])) == 0 :
            mo17_chg_metrate = "NA"
        else : 
            mo17_chg_metrate  = str((float(f2lines[alines.index(aline)].split("\t")[2]))/(float(f2lines[alines.index(aline)].split("\t")[2])+float(f2lines[alines.index(aline)].split("\t")[3])))
        ## oh43_1
        oh43_chg_uniquenoofC = str(f3lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(f3lines[alines.index(aline)].split("\t")[2])+float(f3lines[alines.index(aline)].split("\t")[3])) == 0 :
            oh43_chg_metrate = "NA"
        else : 
            oh43_chg_metrate  = str((float(f3lines[alines.index(aline)].split("\t")[2]))/(float(f3lines[alines.index(aline)].split("\t")[2])+float(f3lines[alines.index(aline)].split("\t")[3])))
        ## cml322_1
        cml322_chg_uniquenoofC = str(f4lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(f4lines[alines.index(aline)].split("\t")[2])+float(f4lines[alines.index(aline)].split("\t")[3])) == 0 :
            cml322_chg_metrate = "NA"
        else : 
            cml322_chg_metrate  = str((float(f4lines[alines.index(aline)].split("\t")[2]))/(float(f4lines[alines.index(aline)].split("\t")[2])+float(f4lines[alines.index(aline)].split("\t")[3])))
        ## tx303_1
        tx303_chg_uniquenoofC = str(f5lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(f5lines[alines.index(aline)].split("\t")[2])+float(f5lines[alines.index(aline)].split("\t")[3])) == 0 :
            tx303_chg_metrate = "NA"
        else : 
            tx303_chg_metrate  = str((float(f5lines[alines.index(aline)].split("\t")[2]))/(float(f5lines[alines.index(aline)].split("\t")[2])+float(f5lines[alines.index(aline)].split("\t")[3])))
        
        ## CHH
        ## b73_1
        b73_chh_uniquenoofC = str(g1lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(g1lines[alines.index(aline)].split("\t")[2])+float(g1lines[alines.index(aline)].split("\t")[3])) == 0 :
            b73_chh_metrate = "NA"
        else : 
            b73_chh_metrate  = str((float(g1lines[alines.index(aline)].split("\t")[2]))/(float(g1lines[alines.index(aline)].split("\t")[2])+float(g1lines[alines.index(aline)].split("\t")[3])))
        ## mo17_1
        mo17_chh_uniquenoofC = str(g2lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(g2lines[alines.index(aline)].split("\t")[2])+float(g2lines[alines.index(aline)].split("\t")[3])) == 0 :
            mo17_chh_metrate = "NA"
        else : 
            mo17_chh_metrate  = str((float(g2lines[alines.index(aline)].split("\t")[2]))/(float(g2lines[alines.index(aline)].split("\t")[2])+float(g2lines[alines.index(aline)].split("\t")[3])))
        ## oh43_1
        oh43_chh_uniquenoofC = str(g3lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(g3lines[alines.index(aline)].split("\t")[2])+float(g3lines[alines.index(aline)].split("\t")[3])) == 0 :
            oh43_chh_metrate = "NA"
        else : 
            oh43_chh_metrate  = str((float(g3lines[alines.index(aline)].split("\t")[2]))/(float(g3lines[alines.index(aline)].split("\t")[2])+float(g3lines[alines.index(aline)].split("\t")[3])))
        ## cml322_1
        cml322_chh_uniquenoofC = str(g4lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(g4lines[alines.index(aline)].split("\t")[2])+float(g4lines[alines.index(aline)].split("\t")[3])) == 0 :
            cml322_chh_metrate = "NA"
        else : 
            cml322_chh_metrate  = str((float(g4lines[alines.index(aline)].split("\t")[2]))/(float(g4lines[alines.index(aline)].split("\t")[2])+float(g4lines[alines.index(aline)].split("\t")[3])))
        ## tx303_1
        tx303_chh_uniquenoofC = str(g5lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(g5lines[alines.index(aline)].split("\t")[2])+float(g5lines[alines.index(aline)].split("\t")[3])) == 0 :
            tx303_chh_metrate = "NA"
        else : 
            tx303_chh_metrate  = str((float(g5lines[alines.index(aline)].split("\t")[2]))/(float(g5lines[alines.index(aline)].split("\t")[2])+float(g5lines[alines.index(aline)].split("\t")[3])))
        
        totalCpG = str(clines[alines.index(aline)].split("\t")[3])
        totalCHG = str(clines[alines.index(aline)].split("\t")[4])
        totalCHH = str(clines[alines.index(aline)].split("\t")[5].strip("\n"))
        b73_cov = str(dlines[alines.index(aline)].split("\t")[4])
        mo17_cov = str(dlines[alines.index(aline)].split("\t")[5])
        oh43_cov = str(dlines[alines.index(aline)].split("\t")[6])
        cml322_cov = str(dlines[alines.index(aline)].split("\t")[7])
        tx303_cov = str(dlines[alines.index(aline)].split("\t")[8].strip("\n"))
        if elines[alines.index(aline)].split("\t")[3] == "NA" :
            probe_count = "0"
        else : 
            probe_count = str(len(elines[alines.index(aline)].split("\t")[3].split(",")))
        b73_probe_val = str(elines[alines.index(aline)].split("\t")[4])
        mo17_probe_val = str(elines[alines.index(aline)].split("\t")[5])
        oh43_probe_val = str(elines[alines.index(aline)].split("\t")[6])
        cml322_probe_val = str(elines[alines.index(aline)].split("\t")[7])
        tx303_probe_val = str(elines[alines.index(aline)].split("\t")[8].strip("\n"))
        outfile.write("\t".join(coordinates)+"\t"+length+"\t"+str(contrib_geno)+"\t"+str(contrib_type)+"\t"+b73_cov+"\t"+mo17_cov+"\t"+oh43_cov+"\t"+cml322_cov+"\t"+tx303_cov+"\t"+b73_metrate+"\t"+mo17_metrate+"\t"+oh43_metrate+"\t"+cml322_metrate+"\t"+tx303_metrate+"\t"+b73_chg_metrate+"\t"+mo17_chg_metrate+"\t"+oh43_chg_metrate+"\t"+cml322_chg_metrate+"\t"+tx303_chg_metrate+"\t"+b73_chh_metrate+"\t"+mo17_chh_metrate+"\t"+oh43_chh_metrate+"\t"+cml322_chh_metrate+"\t"+tx303_chh_metrate+"\t"+b73_uniquenoofC+"\t"+mo17_uniquenoofC+"\t"+oh43_uniquenoofC+"\t"+cml322_uniquenoofC+"\t"+tx303_uniquenoofC+"\t"+b73_chg_uniquenoofC+"\t"+mo17_chg_uniquenoofC+"\t"+oh43_chg_uniquenoofC+"\t"+cml322_chg_uniquenoofC+"\t"+tx303_chg_uniquenoofC+"\t"+b73_chh_uniquenoofC+"\t"+mo17_chh_uniquenoofC+"\t"+oh43_chh_uniquenoofC+"\t"+cml322_chh_uniquenoofC+"\t"+tx303_chh_uniquenoofC+"\t"+probe_count+"\t"+b73_probe_val+"\t"+mo17_probe_val+"\t"+oh43_probe_val+"\t"+cml322_probe_val+"\t"+tx303_probe_val+"\t"+totalCpG+"\t"+totalCHG+"\t"+totalCHH+"\n")
        
    outfile.close()
    #os.system("rm -Rf temp*"+str(chr)+".*")

def bedgraph(chr) :
    a = open("CpG_DMR_5genos_chr1_2.txt", 'r')
    outfile1 = open("CpG_DMR_5genos_chr1_Tx303_CpG_meth_1.bedGraph", 'w')
    outfile2 = open("CpG_DMR_5genos_chr1_Tx303_CpG_meth_2.bedGraph", 'w')
    alines = a.readlines()
    for a in alines[1:]:
        if a.split("\t")[14] != "NA" and a.split("\t")[19] != "NA":
            outfile1.write("\t".join(a.split("\t")[:3]))
            outfile1.write("\t")
            outfile1.write(a.split("\t")[14])
            outfile1.write("\n")
            outfile2.write("\t".join(a.split("\t")[:3]))
            outfile2.write("\t")
            outfile2.write(a.split("\t")[19])
            outfile2.write("\n")
    outfile1.close()
    outfile2.close()

def all(chr) :
    #merge_dmr(chr,context)
    #coverage(chr,context)
    #check_pos(chr)
    #medip(chr,context)
    #medip_tmp(chr,context,tmp)
    #cpg_metrate(chr)
    #chg_metrate(chr)
    #chh_metrate(chr)
    combine(chr)

jobs = []
for chr in range(1,11):
    print chr
    s1 = multiprocessing.Process(target=all, args=(chr, ))
    jobs.append(s1)
    s1.start()
[x.join() for x in jobs]
#coverage(chr)
#cpg_metrate(chr)



outfile = open("Merged_DMR_5genos_08-08-14.txt", 'w')
for chr in range(1,11) :
    print chr
    afile = open("Merged_DMR_5genos_chr"+str(chr)+"_08-08-14.txt", 'r') 
    alines = afile.readlines()
    if chr == 1 :
        for aline in alines :
            outfile.write(aline)
    else :
        for aline in alines[1:] :
            outfile.write(aline)
outfile.close()
"""
for context in ['CHH'] :
    print context
    redo_list = []
    for chr in range(1,11) :
        for geno in genotypes :
            if check_dmr_file(geno, chr, context) not in redo_list and check_dmr_file(geno, chr, context) != "NA" :
                redo_list.append(check_dmr_file(geno, chr, context))
    print redo_list
    if len(redo_list) == 0 :
        jobs = []
        for chr in range(1,11):
            s1 = multiprocessing.Process(target=all, args=(chr, context, ))
            jobs.append(s1)
            s1.start()
        [x.join() for x in jobs]
        print "complete"; sys.exit()
    else :
        print "fix these before continuing", redo_list; sys.exit()


context = "CHH" 
outfile = open(str(context)+"_DMR_5genos_merged_new2.txt", 'w')
for chr in range(1,11) :
    print chr
    if chr == 1 :
        for a in open(str(context)+"_DMR_5genos_chr"+str(chr)+"_new2.txt", 'r') :
            outfile.write(a)
    else :
        a = open(str(context)+"_DMR_5genos_chr"+str(chr)+"_new2.txt", 'r')
        alines = a.readlines()
        for aline in alines[1:] :
            outfile.write(aline)
outfile.close()  
sys.exit()
"""     