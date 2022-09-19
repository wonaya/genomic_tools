import os,sys
import datetime
import random
import multiprocessing
import subprocess
import optparse

def conversion_rate(name) :
    print "extract chloroplast only" 
    print "OT"
    list_pt = []
    for a in open("CHH_OT_"+str(name)+"_pe.txt", 'r') :
        if len(a.split("\t")) > 1 and a.split("\t")[2] == "chloroplast" :
            list_pt.append(a)
    print "CTOT"
    if os.path.isfile("CHH_CTOT_"+str(name)+"_pe.txt") :
        for b in open("CHH_CTOT_"+str(name)+"_pe.txt", 'r') :
            if len(b.split("\t")) > 1 and b.split("\t")[2] == "chloroplast" :
                list_pt.append(b)
    print "OB"
    for c in open("CHH_OB_"+str(name)+"_pe.txt", 'r') :
        if len(c.split("\t")) > 1 and c.split("\t")[2] == "chloroplast" :
            list_pt.append(c)
    print "CTOB"
    if os.path.isfile("CHH_CTOB_"+str(name)+"_pe.txt") :
        for d in open("CHH_CTOB_"+str(name)+"_pe.txt", 'r') :
            if len(d.split("\t")) > 1 and d.split("\t")[2] == "chloroplast" :
                list_pt.append(d)
    large_list = []
    met = 0
    unmet = 0
    for list in list_pt :
        if list.split("\t")[1] == "+" :
            met += 1
        elif list.split("\t")[1] == "-" :
            unmet += 1
    print name, met, unmet, float(met)/float(int(met)+int(unmet))*100

def conversion_rate_5genos(name) :
    list_pt = []
    for a in open("CHH_context_"+str(name)+"_pe_sortedn.txt", 'r') :
        if len(a.split("\t")) > 1 and a.split("\t")[2] == "chloroplast" :
            list_pt.append(a)
    large_list = []
    met = 0
    unmet = 0
    for list in list_pt :
        if list.split("\t")[1] == "+" :
            met += 1
        elif list.split("\t")[1] == "-" :
            unmet += 1
    print name, met, unmet, float(met)/float(int(met)+int(unmet))*100
    
def tile(name, chr, window) :
    ## define genotype and chromosome
     
    ## make array of 8 columns and put 0
    ## CpG/TotC/CHG/TotC/CHH/TotC/all/TotC forward (first 8 columns) and reverse
    
    print "make array", "chr"+str(chr), datetime.datetime.now()
    large_list = []
    for x in range(0,500000000/int(window)):
        large_list.append([0]*16)

    ### CpG
    ### forward strand
    print "read CpG OT", "chr"+str(chr), datetime.datetime.now()
    for CpG_OT in open("CpG_OT_"+str(name)+".txt", 'r') :
        if CpG_OT[0:7] != "Bismark" :
            if CpG_OT.split("\t")[2] == "chr"+str(chr) :
                if CpG_OT.split("\t")[1] == "+" :
                    large_list[int(float(CpG_OT.split("\t")[3])/int(window))][0] += 1 
                    large_list[int(float(CpG_OT.split("\t")[3])/int(window))][1] += 1
                elif CpG_OT.split("\t")[1] == "-" :
                    large_list[int(float(CpG_OT.split("\t")[3])/int(window))][1] += 1
    del CpG_OT
    print "read CpG CTOT", "chr"+str(chr), datetime.datetime.now()
    for CpG_CTOT in open("CpG_CTOT_"+str(name)+".txt", 'r') :
        if CpG_CTOT[0:7] != "Bismark" :
            if CpG_CTOT.split("\t")[2] == "chr"+str(chr) :
                if CpG_CTOT.split("\t")[1] == "+" :
                    large_list[int(float(CpG_CTOT.split("\t")[3])/int(window))][0] += 1 
                    large_list[int(float(CpG_CTOT.split("\t")[3])/int(window))][1] += 1
                elif CpG_CTOT.split("\t")[1] == "-" :
                    large_list[int(float(CpG_CTOT.split("\t")[3])/int(window))][1] += 1
    del CpG_CTOT     
    ### reverse strand 
    print "read CpG OB","chr"+str(chr),datetime.datetime.now()
    for CpG_OB in open("CpG_OB_"+str(name)+".txt", 'r') :
        if CpG_OB[0:7] != "Bismark" :
            if CpG_OB.split("\t")[2] == "chr"+str(chr) :
                if CpG_OB.split("\t")[1] == "+" :
                    large_list[int(float(CpG_OB.split("\t")[3])/int(window))][8] += 1 
                    large_list[int(float(CpG_OB.split("\t")[3])/int(window))][9] += 1
                elif CpG_OB.split("\t")[1] == "-" :
                    large_list[int(float(CpG_OB.split("\t")[3])/int(window))][9] += 1
    del CpG_OB
    print "read CpG CTOB","chr"+str(chr),datetime.datetime.now()
    for CpG_CTOB in open("CpG_CTOB_"+str(name)+".txt", 'r') :
        if CpG_CTOB[0:7] != "Bismark" :
            if CpG_CTOB.split("\t")[2] == "chr"+str(chr) :
                if CpG_CTOB.split("\t")[1] == "+" :
                    large_list[int(float(CpG_CTOB.split("\t")[3])/int(window))][8] += 1 
                    large_list[int(float(CpG_CTOB.split("\t")[3])/int(window))][9] += 1
                elif CpG_CTOB.split("\t")[1] == "-" :
                    large_list[int(float(CpG_CTOB.split("\t")[3])/int(window))][9] += 1
    del CpG_CTOB
    ### CHG
    ### forward strand
    print "read CHG OT","chr"+str(chr),datetime.datetime.now()
    for CHG_OT in open("CHG_OT_"+str(name)+".txt", 'r') :
        if CHG_OT[0:7] != "Bismark" :
            if CHG_OT.split("\t")[2] == "chr"+str(chr) :
                if CHG_OT.split("\t")[1] == "+" :
                    large_list[int(float(CHG_OT.split("\t")[3])/int(window))][2] += 1 
                    large_list[int(float(CHG_OT.split("\t")[3])/int(window))][3] += 1
                elif CHG_OT.split("\t")[1] == "-" :
                    large_list[int(float(CHG_OT.split("\t")[3])/int(window))][3] += 1
    del CHG_OT
    print "read CHG CTOT","chr"+str(chr),datetime.datetime.now()
    for CHG_CTOT in open("CHG_CTOT_"+str(name)+".txt", 'r') :
        if CHG_CTOT[0:7] != "Bismark" :
            if CHG_CTOT.split("\t")[2] == "chr"+str(chr) :
                if CHG_CTOT.split("\t")[1] == "+" :
                    large_list[int(float(CHG_CTOT.split("\t")[3])/int(window))][2] += 1 
                    large_list[int(float(CHG_CTOT.split("\t")[3])/int(window))][3] += 1
                elif CHG_CTOT.split("\t")[1] == "-" :
                    large_list[int(float(CHG_CTOT.split("\t")[3])/int(window))][3] += 1
    del CHG_CTOT    
    ### reverse strand 
    print "read CHG OB","chr"+str(chr),datetime.datetime.now()
    for CHG_OB in open("CHG_OB_"+str(name)+".txt", 'r') :
        if CHG_OB[0:7] != "Bismark" :
            if CHG_OB.split("\t")[2] == "chr"+str(chr) :
                if CHG_OB.split("\t")[1] == "+" :
                    large_list[int(float(CHG_OB.split("\t")[3])/int(window))][10] += 1 
                    large_list[int(float(CHG_OB.split("\t")[3])/int(window))][11] += 1
                elif CHG_OB.split("\t")[1] == "-" :
                    large_list[int(float(CHG_OB.split("\t")[3])/int(window))][11] += 1
    del CHG_OB
    print "read CHG CTOB","chr"+str(chr),datetime.datetime.now()
    for CHG_CTOB in open("CHG_CTOB_"+str(name)+".txt", 'r') :
        if CHG_CTOB[0:7] != "Bismark" :
            if CHG_CTOB.split("\t")[2] == "chr"+str(chr) :
                if CHG_CTOB.split("\t")[1] == "+" :
                    large_list[int(float(CHG_CTOB.split("\t")[3])/int(window))][10] += 1 
                    large_list[int(float(CHG_CTOB.split("\t")[3])/int(window))][11] += 1
                elif CHG_CTOB.split("\t")[1] == "-" :
                    large_list[int(float(CHG_CTOB.split("\t")[3])/int(window))][11] += 1
    del CHG_CTOB

    ### CHH
    ### forward strand
    print "read CHH OT","chr"+str(chr),datetime.datetime.now()
    for CHH_OT in open("CHH_OT_"+str(name)+".txt", 'r') :
        if CHH_OT[0:7] != "Bismark" :
            if CHH_OT.split("\t")[2] == "chr"+str(chr) :
                if CHH_OT.split("\t")[1] == "+" :
                    large_list[int(float(CHH_OT.split("\t")[3])/int(window))][4] += 1 
                    large_list[int(float(CHH_OT.split("\t")[3])/int(window))][5] += 1
                elif CHH_OT.split("\t")[1] == "-" :
                    large_list[int(float(CHH_OT.split("\t")[3])/int(window))][5] += 1
    del CHH_OT
    print "read CHH CTOT","chr"+str(chr),datetime.datetime.now()
    for CHH_CTOT in open("CHH_CTOT_"+str(name)+".txt", 'r') :
        if CHH_CTOT[0:7] != "Bismark" :
            if CHH_CTOT.split("\t")[2] == "chr"+str(chr) :
                if CHH_CTOT.split("\t")[1] == "+" :
                    large_list[int(float(CHH_CTOT.split("\t")[3])/int(window))][4] += 1 
                    large_list[int(float(CHH_CTOT.split("\t")[3])/int(window))][5] += 1
                elif CHH_CTOT.split("\t")[1] == "-" :
                    large_list[int(float(CHH_CTOT.split("\t")[3])/int(window))][5] += 1
    del CHH_CTOT     
    ### reverse strand 
    print "read CHH OB","chr"+str(chr),datetime.datetime.now()
    for CHH_OB in open("CHH_OB_"+str(name)+".txt", 'r') :
        if CHH_OB[0:7] != "Bismark" :
            if CHH_OB.split("\t")[2] == "chr"+str(chr) :
                if CHH_OB.split("\t")[1] == "+" :
                    large_list[int(float(CHH_OB.split("\t")[3])/int(window))][12] += 1 
                    large_list[int(float(CHH_OB.split("\t")[3])/int(window))][13] += 1
                elif CHH_OB.split("\t")[1] == "-" :
                    large_list[int(float(CHH_OB.split("\t")[3])/int(window))][13] += 1
    del CHH_OB
    print "read CHH CTOB","chr"+str(chr),datetime.datetime.now()
    for CHH_CTOB in open("CHH_CTOB_"+str(name)+".txt", 'r') :
        if CHH_CTOB[0:7] != "Bismark" :
            if CHH_CTOB.split("\t")[2] == "chr"+str(chr) :
                if CHH_CTOB.split("\t")[1] == "+" :
                    large_list[int(float(CHH_CTOB.split("\t")[3])/int(window))][12] += 1 
                    large_list[int(float(CHH_CTOB.split("\t")[3])/int(window))][13] += 1
                elif CHH_CTOB.split("\t")[1] == "-" :
                    large_list[int(float(CHH_CTOB.split("\t")[3])/int(window))][13] += 1
    del CHH_CTOB

    print "all metC, metC", "chr"+str(chr),datetime.datetime.now()
    for list in large_list :
        list[6] += list[0]+list[2]+list[4]
        list[7] += list[1]+list[3]+list[5]
        list[14] += list[8]+list[10]+list[12]
        list[15] += list[9]+list[11]+list[13]
    print "write to file", "chr"+str(chr),datetime.datetime.now()
    outfile = open(str(name)+"_tile_chr"+str(chr)+"_"+str(window)+"bp.txt", 'w')
    outfile.write("Start\tEnd\tChrom\tF CpG metC\tF CpG C\tF CHG metC\tF CHG C\tF CHH metC\tF CHH C\tF all metC\tF all C\tR CpG metC\tR CpG C\tR CHG metC\tR CHG C\tR CHH metC\tR CHH C\tR all metC\tR all C\n")
    x = 0
    for list in large_list :
        if list[15] != 0 or list[7] != 0 :
            outfile.write(str(x*int(window))+"\t")
            outfile.write(str((x+1)*int(window)-1)+"\t")
            outfile.write("chr"+str(chr)+"\t")
            list = map(str, list)
            outfile.write("\t".join(list))
            outfile.write("\n")
        x += 1
    outfile.close()    

def split_by_chr(chr, type, name) : ### merged strands
    print "split chromosomes", chr, type, name, datetime.datetime.now()
    outfile = open(str(type)+"_context_"+str(name)+"_chr"+str(chr)+"_forward.out", 'w')
    for a in open(str(type)+"_OT_"+str(name)+".txt", 'r') :
        if len(a.split("\t")) > 2 and str(a.split("\t")[2]) == "chr"+str(chr) :  
            outfile.write(a)
    for a in open(str(type)+"_CTOT_"+str(name)+".txt", 'r') :
        if len(a.split("\t")) > 2 and str(a.split("\t")[2]) == "chr"+str(chr) :  
            outfile.write(a)
    outfile.close()
    outfile = open(str(type)+"_context_"+str(name)+"_chr"+str(chr)+"_reverse.out", 'w')
    for a in open(str(type)+"_OB_"+str(name)+".txt", 'r') :
        if len(a.split("\t")) > 2 and str(a.split("\t")[2]) == "chr"+str(chr) :  
            outfile.write(a)
    for a in open(str(type)+"_CTOB_"+str(name)+".txt", 'r') :
        if len(a.split("\t")) > 2 and str(a.split("\t")[2]) == "chr"+str(chr) :  
            outfile.write(a)
    outfile.close()
    print "finished split", datetime.datetime.now()
    
def pre_methylKit(chr, type, name) :
    print "methylKit prep", chr, type, name
    outfile = open(str(type)+"_context_"+str(name)+"_chr"+str(chr)+".methylKit" ,'w')
    ## save strand status
    for x in range(0, 4) : ## max set to 4 0 to 399Mbp
        large_list = []
        for y in range(0,100000000): ### 56 percent of lonestar node, 2min to make
            large_list.append([0]*4)
        print x*100000000, "to", (x+1)*100000000, datetime.datetime.now(), ": writing array"
        for a_f in open(str(type)+"_OT_"+str(name)+".fq_bismark_bt2_pe.txt", 'r') : # example: CpG_OT_11510-pool_GGTAGC_L004_R1_001_val_1.fq_bismark_bt2_pe.txt
            if len(a_f.split("\t")) > 2 and a_f.split("\t")[2] == "chr"+str(chr) and len(a_f.split("\t")) != 0 and int(a_f.split("\t")[3]) >= x*100000000 and int(a_f.split("\t")[3]) < (x+1)*100000000 :
                if a_f.split("\t")[1] == "+" :
                    large_list[int(float(a_f.split("\t")[3]))-(x*100000000)][0] += 1
                elif a_f.split("\t")[1] == "-" :
                    large_list[int(float(a_f.split("\t")[3]))-(x*100000000)][1] += 1
        for a_f in open(str(type)+"_CTOT_"+str(name)+".fq_bismark_bt2_pe.txt", 'r') : # example: 
            if len(a_f.split("\t")) > 2 and a_f.split("\t")[2] == "chr"+str(chr) and len(a_f.split("\t")) != 0 and int(a_f.split("\t")[3]) >= x*100000000 and int(a_f.split("\t")[3]) < (x+1)*100000000 :
                if a_f.split("\t")[1] == "+" :
                    large_list[int(float(a_f.split("\t")[3]))-(x*100000000)][0] += 1
                elif a_f.split("\t")[1] == "-" :
                    large_list[int(float(a_f.split("\t")[3]))-(x*100000000)][1] += 1
        
        for a_r in open(str(type)+"_OB_"+str(name)+".fq_bismark_bt2_pe.txt", 'r') : 
            if len(a_r.split("\t")) > 2 and a_r.split("\t")[2] == "chr"+str(chr) and len(a_r.split("\t")) != 0 and int(a_r.split("\t")[3]) >= x*100000000 and int(a_r.split("\t")[3]) < (x+1)*100000000 :
                if a_r.split("\t")[1] == "+" :
                    large_list[int(float(a_r.split("\t")[3]))-(x*100000000)][2] += 1
                elif a_r.split("\t")[1] == "-" :
                    large_list[int(float(a_r.split("\t")[3]))-(x*100000000)][3] += 1
        for a_r in open(str(type)+"_CTOB_"+str(name)+".fq_bismark_bt2_pe.txt", 'r') : 
            if len(a_r.split("\t")) > 2 and a_r.split("\t")[2] == "chr"+str(chr) and len(a_r.split("\t")) != 0 and int(a_r.split("\t")[3]) >= x*100000000 and int(a_r.split("\t")[3]) < (x+1)*100000000 :
                if a_r.split("\t")[1] == "+" :
                    large_list[int(float(a_r.split("\t")[3]))-(x*100000000)][2] += 1
                elif a_r.split("\t")[1] == "-" :
                    large_list[int(float(a_r.split("\t")[3]))-(x*100000000)][3] += 1
        z = 0
        for list in large_list :
            if int(list[0])+int(list[1]) != 0 :
                outfile.write("chr"+str(chr)+"."+str(int(z)+(int(x)*100000000))+"\t")
                outfile.write("chr"+str(chr)+"\t")
                outfile.write(str(int(z)+(int(x)*100000000))+"\t")
                outfile.write("F\t") ## change when later Forward and reverse strands are read in separately
                outfile.write(str(int(list[0])+int(list[1]))+"\t")
                outfile.write(str(float(list[0])/float(int(list[0])+int(list[1]))*100)+"\t")
                outfile.write(str(float(list[1])/float(int(list[0])+int(list[1]))*100)+"\n")
            elif int(list[2])+int(list[3]) != 0 :
                outfile.write("chr"+str(chr)+"."+str(int(z)+(int(x)*100000000))+"\t")
                outfile.write("chr"+str(chr)+"\t")
                outfile.write(str(int(z)+(int(x)*100000000))+"\t")
                outfile.write("R\t") ## change when later Forward and reverse strands are read in separately
                outfile.write(str(int(list[2])+int(list[3]))+"\t")
                outfile.write(str(float(list[2])/float(int(list[2])+int(list[3]))*100)+"\t")
                outfile.write(str(float(list[3])/float(int(list[2])+int(list[3]))*100)+"\n")
            z += 1
        del large_list
    outfile.close()

def prep_methylKit2(chr, type, name) :
    print "methylKit prep", chr, type, name
    outfile = open(str(type)+"_context_"+str(name)+"_chr"+str(chr)+".methylKit2" ,'w')
    ## save strand status
    for x in range(0, 4) : ## max set to 4 0 to 399Mbp
        large_list = []
        for y in range(0,100000000): ### 56 percent of lonestar node, 2min to make
            large_list.append([0]*4)
        print x*100000000, "to", (x+1)*100000000, datetime.datetime.now(), ": writing array"
        for a_f in open(str(type)+"_OT_"+str(name)+"_bt202.txt", 'r') : # example: CpG_OT_11510-pool_GGTAGC_L004_R1_001_val_1.fq_bismark_bt2_pe.txt
            if len(a_f.split("\t")) > 2 and a_f.split("\t")[2] == "chr"+str(chr) and len(a_f.split("\t")) != 0 and int(a_f.split("\t")[3]) >= x*100000000 and int(a_f.split("\t")[3]) < (x+1)*100000000 :
                if a_f.split("\t")[1] == "+" :
                    large_list[int(float(a_f.split("\t")[3]))-(x*100000000)][0] += 1
                elif a_f.split("\t")[1] == "-" :
                    large_list[int(float(a_f.split("\t")[3]))-(x*100000000)][1] += 1
        for a_f in open(str(type)+"_CTOT_"+str(name)+"_bt202.txt", 'r') : # example: 
            if len(a_f.split("\t")) > 2 and a_f.split("\t")[2] == "chr"+str(chr) and len(a_f.split("\t")) != 0 and int(a_f.split("\t")[3]) >= x*100000000 and int(a_f.split("\t")[3]) < (x+1)*100000000 :
                if a_f.split("\t")[1] == "+" :
                    large_list[int(float(a_f.split("\t")[3]))-(x*100000000)][0] += 1
                elif a_f.split("\t")[1] == "-" :
                    large_list[int(float(a_f.split("\t")[3]))-(x*100000000)][1] += 1
        
        for a_r in open(str(type)+"_OB_"+str(name)+"_bt202.txt", 'r') : 
            if len(a_r.split("\t")) > 2 and a_r.split("\t")[2] == "chr"+str(chr) and len(a_r.split("\t")) != 0 and int(a_r.split("\t")[3]) >= x*100000000 and int(a_r.split("\t")[3]) < (x+1)*100000000 :
                if a_r.split("\t")[1] == "+" :
                    large_list[int(float(a_r.split("\t")[3]))-(x*100000000)][2] += 1
                elif a_r.split("\t")[1] == "-" :
                    large_list[int(float(a_r.split("\t")[3]))-(x*100000000)][3] += 1
        for a_r in open(str(type)+"_CTOB_"+str(name)+"_bt202.txt", 'r') : 
            if len(a_r.split("\t")) > 2 and a_r.split("\t")[2] == "chr"+str(chr) and len(a_r.split("\t")) != 0 and int(a_r.split("\t")[3]) >= x*100000000 and int(a_r.split("\t")[3]) < (x+1)*100000000 :
                if a_r.split("\t")[1] == "+" :
                    large_list[int(float(a_r.split("\t")[3]))-(x*100000000)][2] += 1
                elif a_r.split("\t")[1] == "-" :
                    large_list[int(float(a_r.split("\t")[3]))-(x*100000000)][3] += 1
        z = 0
        for list in large_list :
            if int(list[0])+int(list[1]) != 0 :
                outfile.write(str(int(z)+(int(x)*100000000))+"\t")
                outfile.write(str(chr)+"\t")
                outfile.write(str(int(z)+(int(x)*100000000))+"\t")
                outfile.write("1\t") ## change when later Forward and reverse strands are read in separately
                outfile.write(str(int(list[0])+int(list[1]))+"\t")
                outfile.write(str(float(list[0])/float(int(list[0])+int(list[1]))*100)+"\t")
                outfile.write(str(float(list[1])/float(int(list[0])+int(list[1]))*100)+"\n")
            elif int(list[2])+int(list[3]) != 0 :
                outfile.write(str(int(z)+(int(x)*100000000))+"\t")
                outfile.write(str(chr)+"\t")
                outfile.write(str(int(z)+(int(x)*100000000))+"\t")
                outfile.write("0\t") ## change when later Forward and reverse strands are read in separately
                outfile.write(str(int(list[2])+int(list[3]))+"\t")
                outfile.write(str(float(list[2])/float(int(list[2])+int(list[3]))*100)+"\t")
                outfile.write(str(float(list[3])/float(int(list[2])+int(list[3]))*100)+"\n")
            z += 1
        del large_list
    outfile.close()

def methylkit(name1, name2, chrom, type) :
    os.system("Rscript /work/02114/wonaya/scripts/r_methylkit.R "+str(name1)+"/"+str(type)+"_context_"+str(name1)+"_bt202_chr"+str(chr)+".methylKit "+str(name2)+"/"+str(type)+"_context_"+str(name2)+"_bt202_chr"+str(chr)+".methylKit "+str(name1.split("_all3/")[0])+" "+str(name2.split("_all3/")[0])+" "+str(type)+" "+str(name1.split("_all3/")[0])+"_"+str(name2.split("_all3/")[0])+"_"+str(type)+"_dmr_chr"+str(chr)+".txt")

def merge_chrom(name1, name2, type) :
    #print "merge chromsomes", type, name1, name2
    ## check if there's all 10 files
    file_list = []
    for files in os.listdir("."): 
        if "_".join(files.split("_")[0:2]) == name1 :
            if "_".join(files.split("_")[2:4]) == name2 : 
                if files.split("_")[4] == type :
                    for x in range(1,11) :
                        if files.split("_")[6] == "chr"+str(x)+".txt" :
                            file_list.append(files)
    print type, name1, name2, len(file_list)
    if len(file_list) == 10 :
        print len(file_list), file_list
        outfile = open(str(name1.split("_all3/")[0])+"_"+str(name2.split("_all3/")[0])+"_"+str(type)+"_dmr_merged.txt", 'w')
        outfile.write("chrom\tstart\tend\twidth\tmean.meth.diff\tnum."+str(type)+"\tnum.DMCs\tDMR.pvalue\tDMR.qvalue\n")
        for file in file_list :
            file_line = open(file, 'r')
            lines = file_line.readlines()
            for line in lines[1:] :
                outfile.write(line.split("\t")[1].strip('"'))
                outfile.write("\t")
                outfile.write(line.split("\t")[2])
                outfile.write("\t")
                outfile.write(line.split("\t")[3])
                outfile.write("\t")
                outfile.write(line.split("\t")[4])
                outfile.write("\t")
                outfile.write(line.split("\t")[6])
                outfile.write("\t")
                outfile.write(line.split("\t")[7])
                outfile.write("\t")
                outfile.write(line.split("\t")[8])
                outfile.write("\t")
                outfile.write(line.split("\t")[9])
                outfile.write("\t")
                outfile.write(line.split("\t")[10])
        outfile.close()
        os.system("rm -Rf "+str(name1.split("_all3/")[0])+"_"+str(name2.split("_all3/")[0])+"_"+str(type)+"_dmr_chr*.txt")
    
if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("--run", default=None, help="tile merge bedgraph circos")
    parser.add_option("--name", default=None, help="select from: B73_all3, Mo17_all3, Oh43_all3, CML322_all3, Tx303_all3 (basically name of bam used for the BME)")
    parser.add_option("--name1", default=None, help="select from: B73_all3, Mo17_all3, Oh43_all3, CML322_all3, Tx303_all3 (basically name of bam used for the BME)")
    parser.add_option("--name2", default=None, help="select from: B73_all3, Mo17_all3, Oh43_all3, CML322_all3, Tx303_all3 (basically name of bam used for the BME)")
    parser.add_option("--window", default=100, help="specify window size 100, 50, 10000")
    parser.add_option("--specie", default="maize", help="select from: arabidopsis, maize")
    parser.add_option("--type", default=None, help="use with bedGraph CpG, CHH, CHG or all")
    parser.add_option("--direction", default=None, help="use with bedGraph reverse, forward or both")
    parser.add_option("--chrom", default=None, help="chromosome 1 to 10")
    options, args = parser.parse_args()
    
    if options.run == "tile" :
        jobs = []
        for chr in range(1,11):
            o = multiprocessing.Process(target=tile, args=(options.name, chr, options.window,))
            jobs.append(o)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
        
    elif options.run == "tile_fix" :
        outfile = open(str(options.name)+"_tile_merged_"+str(options.window)+"bp_fix.txt", 'w')
        for a in open(str(options.name)+"_tile_merged_"+str(options.window)+"bp.txt", 'r') :
            if a.split("\t")[0] == "Start" :
                outfile.write(a)
            else :
                outfile.write(a.split("\t")[0])
                outfile.write("\t")
                outfile.write(str(int(a.split("\t")[1])-1))
                outfile.write("\t")
                outfile.write("\t".join(map(str, a.split("\t")[2:])))
        outfile.close()
    
    elif options.run == "merge" :
        outfile = open(str(options.name)+"_tile_merged_"+str(options.window)+"bp.txt", 'w')
        outfile.write("Start\tEnd\tChrom\tF CpG metC\tF CpG C\tF CHG metC\tF CHG C\tF CHH metC\tF CHH C\tF all metC\tF all C\tR CpG metC\tR CpG C\tR CHG metC\tR CHG C\tR CHH metC\tR CHH C\tR all metC\tR all C\n")
        for chr in range(1,11) :
            print "writing chr"+str(chr)
            for a in open(str(options.name)+"_tile_chr"+str(chr)+"_"+str(options.window)+"bp.txt", 'r') :
                if a.split("\t")[0] != "Start" :
                    outfile.write(a)
        outfile.close()
        os.system("rm -Rf "+str(options.name)+"_tile_chr*")
    
    elif options.run == "bedgraph" :
        outfile = open(str(options.name)+"_tile_merged_"+str(options.type)+"_"+str(options.direction)+"_"+str(options.window)+"bp.bedGraph", 'w')
        for a in open(str(options.name)+"_tile_merged_"+str(options.window)+"bp.txt", 'r') :
            if a.split("\t")[0] != "Start" :
                if float(a.split("\t")[10]) != 0 :
                    if options.type == "CpG" :
                        if options.direction == "forward" :
                            if float(a.split("\t")[4]) != 0 :
                                outfile.write(str(a.split("\t")[2])+"\t")
                                outfile.write(str(a.split("\t")[0])+"\t")
                                outfile.write(str(a.split("\t")[1])+"\t")
                                outfile.write(str(float(float(a.split("\t")[3])/float(a.split("\t")[4]))*100))
                                outfile.write("\t+\n")
                        elif options.direction == "reverse" :
                            if float(a.split("\t")[12]) != 0 :
                                outfile.write(str(a.split("\t")[2])+"\t")
                                outfile.write(str(a.split("\t")[0])+"\t")
                                outfile.write(str(a.split("\t")[1])+"\t")
                                outfile.write(str(float(float(a.split("\t")[11])/float(a.split("\t")[12]))*100))
                                outfile.write("\t-\n")
                        elif options.direction == "both" :
                            if float(int(a.split("\t")[4])+int(a.split("\t")[12])) != 0 :
                                outfile.write(str(a.split("\t")[2])+"\t")
                                outfile.write(str(a.split("\t")[0])+"\t")
                                outfile.write(str(a.split("\t")[1])+"\t")
                                outfile.write(str(float(float(int(a.split("\t")[3])+int(a.split("\t")[11]))/float(int(a.split("\t")[4])+int(a.split("\t")[12])))*100))
                                outfile.write("\n")
                        else : print "check direction option" ; sys.exit()
                    elif options.type == "CHG" :
                        if options.direction == "forward" :
                            if float(a.split("\t")[6]) != 0 :
                                outfile.write(str(a.split("\t")[2])+"\t")
                                outfile.write(str(a.split("\t")[0])+"\t")
                                outfile.write(str(a.split("\t")[1])+"\t")
                                outfile.write(str(float(float(a.split("\t")[5])/float(a.split("\t")[6]))*100))
                                outfile.write("\t+\n")
                        elif options.direction == "reverse" :
                            if float(a.split("\t")[14]) != 0 :
                                outfile.write(str(a.split("\t")[2])+"\t")
                                outfile.write(str(a.split("\t")[0])+"\t")
                                outfile.write(str(a.split("\t")[1])+"\t")
                                outfile.write(str(float(float(a.split("\t")[13])/float(a.split("\t")[14]))*100))
                                outfile.write("\t-\n")
                        elif options.direction == "both" :
                            if float(int(a.split("\t")[6])+int(a.split("\t")[14])) != 0 :
                                outfile.write(str(a.split("\t")[2])+"\t")
                                outfile.write(str(a.split("\t")[0])+"\t")
                                outfile.write(str(a.split("\t")[1])+"\t")
                                outfile.write(str(float(float(int(a.split("\t")[5])+int(a.split("\t")[13]))/float(int(a.split("\t")[6])+int(a.split("\t")[14])))*100))
                                outfile.write("\n")
                        else : print "check direction option" ; sys.exit()
                    elif options.type == "CHH" :
                        if options.direction == "forward" :
                            if float(a.split("\t")[8]) != 0 :
                                outfile.write(str(a.split("\t")[2])+"\t")
                                outfile.write(str(a.split("\t")[0])+"\t")
                                outfile.write(str(a.split("\t")[1])+"\t")
                                outfile.write(str(float(float(a.split("\t")[7])/float(a.split("\t")[8]))*100))
                                outfile.write("\t+\n")
                        elif options.direction == "reverse" :
                            if float(a.split("\t")[16]) != 0 :
                                outfile.write(str(a.split("\t")[2])+"\t")
                                outfile.write(str(a.split("\t")[0])+"\t")
                                outfile.write(str(a.split("\t")[1])+"\t")
                                outfile.write(str(float(float(a.split("\t")[15])/float(a.split("\t")[16]))*100))
                                outfile.write("\t-\n")
                        elif options.direction == "both" :
                            if float(int(a.split("\t")[8])+int(a.split("\t")[16])) != 0 :
                                outfile.write(str(a.split("\t")[2])+"\t")
                                outfile.write(str(a.split("\t")[0])+"\t")
                                outfile.write(str(a.split("\t")[1])+"\t")
                                outfile.write(str(float(float(int(a.split("\t")[7])+int(a.split("\t")[15]))/float(int(a.split("\t")[8])+int(a.split("\t")[16])))*100))
                                outfile.write("\n")
                        else : print "check direction option" ; sys.exit()
                    elif options.type == "all" :
                        if options.direction == "forward" :
                            if float(a.split("\t")[10]) != 0 :
                                outfile.write(str(a.split("\t")[2])+"\t")
                                outfile.write(str(a.split("\t")[0])+"\t")
                                outfile.write(str(a.split("\t")[1])+"\t")
                                outfile.write(str(float(float(a.split("\t")[9])/float(a.split("\t")[10]))*100))
                                outfile.write("\t+\n")
                        elif options.direction == "reverse" :
                            if float(a.split("\t")[18]) != 0 :
                                outfile.write(str(a.split("\t")[2])+"\t")
                                outfile.write(str(a.split("\t")[0])+"\t")
                                outfile.write(str(a.split("\t")[1])+"\t")
                                outfile.write(str(float(float(a.split("\t")[17])/float(a.split("\t")[18]))*100))
                                outfile.write("\t-\n")
                        elif options.direction == "both" :
                            if float(int(a.split("\t")[10])+int(a.split("\t")[18])) != 0 :
                                outfile.write(str(a.split("\t")[2])+"\t")
                                outfile.write(str(a.split("\t")[0])+"\t")
                                outfile.write(str(a.split("\t")[1])+"\t")
                                outfile.write(str(float(float(int(a.split("\t")[9])+int(a.split("\t")[17]))/float(int(a.split("\t")[10])+int(a.split("\t")[18])))*100))
                                outfile.write("\n")
                        else : print "check direction option" ; sys.exit()
                    else :
                        print "check your options for bedgraph" ; sys.exit()
        outfile.close()
        
    elif options.run == "circos" :
        outfile = open(str(options.name)+"_tile_merged_"+str(options.type)+"_"+str(options.direction)+"_"+str(options.window)+"bp.circos", 'w')
        #for a in open(str(options.name)+"_tile_merged.txt", 'r') :
        for a in open(str(options.name)+"_tile_merged_"+str(options.type)+"_"+str(options.direction)+"_"+str(options.window)+"bp.bedGraph", 'r') :
            outfile.write("zm"+str(a.split("\t")[0][3:])+" ")
            outfile.write(str(a.split("\t")[1])+" "+str(a.split("\t")[2])+" "+str(a.split("\t")[3].strip("\n")))
            outfile.write("\n")
        outfile.close()
    
    elif options.run == "igvtools" :
        os.system("/opt/apps/igvtools/2.1.7/igvtools toTDF "+str(options.name)+"_tile_merged_CpG_"+str(options.direction)+"_"+str(options.window)+"bp.bedGraph "+str(options.name)+"_tile_merged_CpG_"+str(options.direction)+"_"+str(options.window)+"bp.tdf B73") 
        os.system("/opt/apps/igvtools/2.1.7/igvtools toTDF "+str(options.name)+"_tile_merged_CHH_"+str(options.direction)+"_"+str(options.window)+"bp.bedGraph "+str(options.name)+"_tile_merged_CHH_"+str(options.direction)+"_"+str(options.window)+"bp.tdf B73") 
        os.system("/opt/apps/igvtools/2.1.7/igvtools toTDF "+str(options.name)+"_tile_merged_CHG_"+str(options.direction)+"_"+str(options.window)+"bp.bedGraph "+str(options.name)+"_tile_merged_CHG_"+str(options.direction)+"_"+str(options.window)+"bp.tdf B73") 
        os.system("/opt/apps/igvtools/2.1.7/igvtools toTDF "+str(options.name)+"_tile_merged_all_"+str(options.direction)+"_"+str(options.window)+"bp.bedGraph "+str(options.name)+"_tile_merged_all_"+str(options.direction)+"_"+str(options.window)+"bp.tdf B73") 
    
    elif options.run == "split_chr" :
        split_by_chr(options.chrom, options.type, options.name)
        #os.system("/work/02114/wonaya/software/bismark_v0.10.0/bismark2bedGraph --counts --CX_context -o "+str(options.type)+"_context_"+str(options.name)+"_chr"+str(options.chrom)+".bedGraph "+str(options.type)+"_context_"+str(options.name)+"_chr"+str(options.chrom)+".out")
    
    elif options.run == "test" :
        #names = ['11147-pool_ATGAGC_L003', '11456-pool_CAAAAG_L003', '11457-pool_CAACTA_L003','11471-pool_CACCGG_L003','11472-pool_CACTCA_L003','11477-pool_CAGGCG_L003','11481-pool_CATGGC_L003','11489-pool_GCGCTA_L003','11493-pool_GTGCTA_L004','11510-pool_GGTAGC_L004','11512-pool_CAACTA_L004','KM-mop1_CAAAAG_L004','KM-tgr1_ATGAGC_L004','KM-tgr9_TATAAT_L004','UM87-pool_CACGAT_L004']
        #names = ['Mo17_all3','Oh43_all3','Tx303_all3','CML322_all3','B73_all3']
        #names = ['B73_all3']
        #for name in names :
        #    conversion_rate(name)
        prep_methylKit2(1m, "CHH", "Tx303_all3")
        
    elif options.run == "prep" :
        jobs = []
        if options.specie == "maize" :
            chrmax = 11
        for name in options.name.split(",") :
            for met in options.type.split(",") :
                for chr in range(10,11) :
                    s = multiprocessing.Process(target=pre_methylKit, args=(chr, met, name, ))
                    jobs.append(s)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
    
    elif options.run == "prep_seq" :
        for chr in range(1,11) :
            print chr
            pre_methylKit(chr, options.type, options.name)
    elif options.run == "methylkit" :
        for chr in range(6,11) :
            methylkit(options.name1, options.name2, chr, options.type)
    elif options.run == "merge_methylkit" :
        names = ['Mo17_all3','Oh43_all3','Tx303_all3','CML322_all3','B73_all3']
        types = ['CHG','CHH','CpG']
        for name1 in names :
            for name2 in names :
                if name1 != name2 :
                    for type in types :
                        merge_chrom(name1, name2, type)