### use largemem

import os,sys
import numpy
import multiprocessing

def merge(chr) :
    print "###generating temp DMR list for chr", chr
    """
    outfile = open("temp_alldmr_"+str(chr)+".txt" , 'w')
    dmr_loc1 = []
    if os.path.isfile("B73_Mo17_CpG_dmr_chr"+str(chr)+".txt") :
        for line1 in open("B73_Mo17_CpG_dmr_chr"+str(chr)+".txt", 'r') :
            if line1.split("\t")[0] != '"seqnames"' and line1.split("\t")[1] == '"'+str(chr)+'"' :
                dmr_loc1.append(str(line1.split("\t")[2])+"-"+str(int(line1.split("\t")[3])+1))
                outfile.write("B73\tMo17\t"+str(line1.split("\t")[2])+"\t"+str(int(line1.split("\t")[3]))+"\n")
    dmr_loc2 = []                    
    if os.path.isfile("B73_Oh43_CpG_dmr_chr"+str(chr)+".txt") :
        for line2 in open("B73_Oh43_CpG_dmr_chr"+str(chr)+".txt", 'r') :
            if line2.split("\t")[0] != '"seqnames"' and line2.split("\t")[1] == '"'+str(chr)+'"' :
                dmr_loc2.append(str(line2.split("\t")[2])+"-"+str(int(line2.split("\t")[3])+1))
                outfile.write("B73\tOh43\t"+str(line2.split("\t")[2])+"\t"+str(int(line2.split("\t")[3]))+"\n")
    new_dmr = []
    no_overlap1 = []    
    no_overlap2 = []  
    with_overlap1 = []
    with_overlap2 = []  
    print "B73_Mo17,Oh43_CpG", len(dmr_loc1), len(dmr_loc2)   
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
    if os.path.isfile("B73_CML322_CpG_dmr_chr"+str(chr)+".txt") :
        for line3 in open("B73_CML322_CpG_dmr_chr"+str(chr)+".txt", 'r') :
            if line3.split("\t")[0] != '"seqnames"' and line3.split("\t")[1] == '"'+str(chr)+'"' :
                dmr_loc2.append(str(line3.split("\t")[2])+"-"+str(int(line3.split("\t")[3])+1))
                outfile.write("B73\tCML322\t"+str(line3.split("\t")[2])+"\t"+str(int(line3.split("\t")[3]))+"\n")
    print "B73_Mo17,Oh43,CML322_CpG", len(dmr_loc1), len(dmr_loc2)
    
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
    if os.path.isfile("B73_Tx303_CpG_dmr_chr"+str(chr)+".txt") :
        for line1 in open("B73_Tx303_CpG_dmr_chr"+str(chr)+".txt", 'r') :
            if line1.split("\t")[0] != '"seqnames"' and line1.split("\t")[1] == '"'+str(chr)+'"' :
                dmr_loc2.append(str(line1.split("\t")[2])+"-"+str(int(line1.split("\t")[3])+1))
                outfile.write("B73\tTx303\t"+str(line1.split("\t")[2])+"\t"+str(int(line1.split("\t")[3]))+"\n")
    
    print "B73_Mo17,Oh43,CML322,Tx303_CpG", len(dmr_loc1), len(dmr_loc2)
    
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
    #dmr_loc1 = new_dmr
    
    dmr_loc2 = []
    if os.path.isfile("B73_CML322_CHG_dmr_chr"+str(chr)+".txt") :
        for line2 in open("B73_CML322_CHG_dmr_chr"+str(chr)+".txt", 'r') :
            if line2.split("\t")[0] != '"seqnames"' and line2.split("\t")[1] == '"'+str(chr)+'"' :
                dmr_loc2.append(str(line2.split("\t")[2])+"-"+str(int(line2.split("\t")[3])+1))
                
    print "B73_Mo17,CML322_CpG_CHG", len(dmr_loc1), len(dmr_loc2)
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
    if os.path.isfile("B73_CML322_CHH_dmr_chr"+str(chr)+".txt") :
        for line3 in open("B73_CML322_CHH_dmr_chr"+str(chr)+".txt", 'r') :
            if line3.split("\t")[0] != '"seqnames"' and line3.split("\t")[1] == '"'+str(chr)+'"' :
                dmr_loc2.append(str(line3.split("\t")[2])+"-"+str(int(line3.split("\t")[3])+1))
    
    print "B73_Mo17,CML322_CpG_CHG_CHH", len(dmr_loc1), len(dmr_loc2)
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
    if os.path.isfile("B73_Oh43_CpG_dmr_chr"+str(chr)+".txt") :
        for line3 in open("B73_Oh43_CpG_dmr_chr"+str(chr)+".txt", 'r') :
            if line3.split("\t")[0] != '"seqnames"' and line3.split("\t")[1] == '"'+str(chr)+'"' :
                dmr_loc2.append(str(line3.split("\t")[2])+"-"+str(int(line3.split("\t")[3])+1))
    
    print "B73_Mo17_CML322, Oh43_CpG", len(dmr_loc1), len(dmr_loc2)
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
    if os.path.isfile("B73_Oh43_CHG_dmr_chr"+str(chr)+".txt") :
        for line3 in open("B73_Oh43_CHG_dmr_chr"+str(chr)+".txt", 'r') :
            if line3.split("\t")[0] != '"seqnames"' and line3.split("\t")[1] == '"'+str(chr)+'"' :
                dmr_loc2.append(str(line3.split("\t")[2])+"-"+str(int(line3.split("\t")[3])+1))
    
    print "B73_Mo17_CML322, Oh43_CpG_CHG", len(dmr_loc1), len(dmr_loc2)
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
    if os.path.isfile("B73_Oh43_CHH_dmr_chr"+str(chr)+".txt") :
        for line3 in open("B73_Oh43_CHH_dmr_chr"+str(chr)+".txt", 'r') :
            if line3.split("\t")[0] != '"seqnames"' and line3.split("\t")[1] == '"'+str(chr)+'"' :
                dmr_loc2.append(str(line3.split("\t")[2])+"-"+str(int(line3.split("\t")[3])+1))
    
    print "B73_Mo17_CML322, Oh43_CpG_CHG_CHH", len(dmr_loc1), len(dmr_loc2)
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
    if os.path.isfile("B73_Tx303_CpG_dmr_chr"+str(chr)+".txt") :
        for line3 in open("B73_Tx303_CpG_dmr_chr"+str(chr)+".txt", 'r') :
            if line3.split("\t")[0] != '"seqnames"' and line3.split("\t")[1] == '"'+str(chr)+'"' :
                dmr_loc2.append(str(line3.split("\t")[2])+"-"+str(int(line3.split("\t")[3])+1))
    
    print "B73_Mo17_CML322_Oh43, Tx303_CpG", len(dmr_loc1), len(dmr_loc2)
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
    if os.path.isfile("B73_Tx303_CHG_dmr_chr"+str(chr)+".txt") :
        for line3 in open("B73_Tx303_CHG_dmr_chr"+str(chr)+".txt", 'r') :
            if line3.split("\t")[0] != '"seqnames"' and line3.split("\t")[1] == '"'+str(chr)+'"' :
                dmr_loc2.append(str(line3.split("\t")[2])+"-"+str(int(line3.split("\t")[3])+1))
    
    print "B73_Mo17_CML322_Oh43, Tx303_CpG_CHG", len(dmr_loc1), len(dmr_loc2)
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
    if os.path.isfile("B73_Tx303_CHH_dmr_chr"+str(chr)+".txt") :
        for line3 in open("B73_Tx303_CHH_dmr_chr"+str(chr)+".txt", 'r') :
            if line3.split("\t")[0] != '"seqnames"' and line3.split("\t")[1] == '"'+str(chr)+'"' :
                dmr_loc2.append(str(line3.split("\t")[2])+"-"+str(int(line3.split("\t")[3])+1))
    
    print "B73_Mo17_CML322_Oh43, Tx303_CpG_CHG_CHH", len(dmr_loc1), len(dmr_loc2)
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
    
    print "B73_Mo17_CML322_Oh43_Tx303, CpG", chr, len(sorted_dmr), len(dmr_dict)
    print "###calculating reads mapped for DMR regions"
    ### get coverage
    outfile = open("temp_"+str(chr)+".bed", 'w')
    for dmr in sorted_dmr :
        outfile.write("chr"+str(chr)+"\t"+str(dmr.split("-")[0])+"\t"+str(dmr.split("-")[1])+"\t0\n")
    outfile.close()
    
    os.system("/opt/apps/bedtools/2.17.0/bin/multiBamCov -bams B73_all3/B73_all3_R1_bt202_sorted.bam Mo17_all3/Mo17_all3_bt202_sorted.bam Oh43_all3/Oh43_all3_bt202_sorted.bam CML322_all3/CML322_all3_bt202_sorted.bam Tx303_all3/Tx303_all3_bt202_sorted.bam -bed temp_"+str(chr)+".bed > temp_bam_count_"+str(chr)+".txt")
    
    print "###calculating length of DMR"
    ### calculate length
    outfile = open("temp_"+str(chr)+".length", 'w')
    for a in open("temp_"+str(chr)+".bed", 'r') :
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(int(a.split("\t")[2])-int(a.split("\t")[1]))+"\n")
    outfile.close()
    
    ### calculate non-covered C
    print "###calculating total number of possible methylation sites"
    outfile = open("temp_"+str(chr)+".allC", 'w')
    genome_list = []
    whole_list = []
    chr_list = []
    print "reading in genome"
    for genome in open("/work/02114/wonaya/genome/Zea_mays.AGPv2.14/Zea_mays.AGPv2.14.fa", 'r') :
        if genome[0] == ">" :
            if len(chr_list) == 0 :
                chr_list.append(genome.strip("\n"))
            else :
                whole_list.append("".join(genome_list))
                chr_list.append(genome.strip("\n"))
                genome_list = []
        else :
            genome_list.append(genome.strip("\n"))
    whole_list.append("".join(genome_list))
    
    print chr_list[int(chr)-1]
    
    chh_list = ["CAA", "CAC", "CAT", "CCA", "CCC", "CCT", "CTA", "CTC", "CTT"]
    chh_rev_list = ["TTG", "GTG", "ATG", "TGG", "GGG", "AGG", "TAG", "GAG", "AAG"]
    chg_list = ["CAG", "CTG", "CCG"]
    chg_rev_list = ["CTG", "CAG", "CGG"]
    for a in open("temp_"+str(chr)+".bed", 'r') :
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
        count_cpg += whole_list[int(chr)-1][int(a.split("\t")[1])-1:int(a.split("\t")[2])].count("GC")
        outfile.write("\t".join(a.split("\t")[:3])+"\t"+str(count_cpg)+"\t"+str(count_chg)+"\t"+str(count_chh)+"\n")
    outfile.close()
    """
    ### B73 % methylation
    print "###getting methylation rates for individual DMRs"
    print "B73"
    large_list = []
    #for x in range(0,310000000):
    #   large_list.append([0]*3)
    cpg_list = []
    for a in open("Oh43_all3/CpG_context_Oh43_all3_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        if int(a.split("\t")[1]) >= 86795 and int(a.split("\t")[1]) <= 87511 :
            if chr == 1 :
                cpg_list.append(int(a.split("\t")[1]))
        if int(a.split("\t")[1]) >= 90000 :
            break
    print "Finished"
    chg_list = []
    for a in open("Oh43_all3/CHG_context_Oh43_all3_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        if int(a.split("\t")[1]) >= 86795 and int(a.split("\t")[1]) <= 87511 :
            if chr == 1 :
                chg_list.append(int(a.split("\t")[1]))
        if int(a.split("\t")[1]) >= 90000 :
            break
    print len(cpg_list), len(chg_list)
    print len(set(cpg_list)&set(chg_list))            
    sys.exit()
    
    for a in open("B73_all3/CpG_context_B73_all3_R1_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
       
    outfile = open("temp_b73_"+str(chr)+".cov", 'w')
    for a in open("temp_"+str(chr)+".bed", 'r') :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for region in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
            countmet += region[0]
            countunmet += region[1]
            uniquec += region[2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "Mo17"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    
    for a in open("Mo17_all3/CpG_context_Mo17_all3_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
       
    outfile = open("temp_mo17_"+str(chr)+".cov", 'w')
    for a in open("temp_"+str(chr)+".bed", 'r') :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for region in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
            countmet += region[0]
            countunmet += region[1]
            uniquec += region[2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "CML322"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    
    for a in open("CML322_all3/CpG_context_CML322_all3_bt202_sortedn_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
       
    outfile = open("temp_cml322_"+str(chr)+".cov", 'w')
    for a in open("temp_"+str(chr)+".bed", 'r') :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for region in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
            countmet += region[0]
            countunmet += region[1]
            uniquec += region[2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "Oh43"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    
    for a in open("Oh43_all3/CpG_context_Oh43_all3_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
       
    outfile = open("temp_oh43_"+str(chr)+".cov", 'w')
    for a in open("temp_"+str(chr)+".bed", 'r') :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for region in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
            countmet += region[0]
            countunmet += region[1]
            uniquec += region[2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "Tx303"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    
    for a in open("Tx303_all3/CpG_context_Tx303_all3_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
       
    outfile = open("temp_tx303_"+str(chr)+".cov", 'w')
    for a in open("temp_"+str(chr)+".bed", 'r') :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for region in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
            countmet += region[0]
            countunmet += region[1]
            uniquec += region[2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    
    print "###start meDIP analysis"
    ### meDIP
    outfile = open("temp_meDIP_"+str(chr)+".txt", 'w')
    for a in open("temp_"+str(chr)+".bed", 'r') :
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
    
    
    ### B73 CHG % methylation
    print "### calculate CHG methylation rate"
    print "B73"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    
    for a in open("B73_all3/CHG_context_B73_all3_R1_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
       
    outfile = open("temp_chg_b73_"+str(chr)+".cov", 'w')
    for a in open("temp_"+str(chr)+".bed", 'r') :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for region in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
            countmet += region[0]
            countunmet += region[1]
            uniquec += region[2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "Mo17"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    
    for a in open("Mo17_all3/CHG_context_Mo17_all3_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
       
    outfile = open("temp_chg_mo17_"+str(chr)+".cov", 'w')
    for a in open("temp_"+str(chr)+".bed", 'r') :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for region in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
            countmet += region[0]
            countunmet += region[1]
            uniquec += region[2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "CML322"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    
    for a in open("CML322_all3/CHG_context_CML322_all3_bt202_sortedn_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
       
    outfile = open("temp_chg_cml322_"+str(chr)+".cov", 'w')
    for a in open("temp_"+str(chr)+".bed", 'r') :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for region in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
            countmet += region[0]
            countunmet += region[1]
            uniquec += region[2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "Oh43"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    
    for a in open("Oh43_all3/CHG_context_Oh43_all3_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
       
    outfile = open("temp_chg_oh43_"+str(chr)+".cov", 'w')
    for a in open("temp_"+str(chr)+".bed", 'r') :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for region in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
            countmet += region[0]
            countunmet += region[1]
            uniquec += region[2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "Tx303"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    
    for a in open("Tx303_all3/CHG_context_Tx303_all3_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
       
    outfile = open("temp_chg_tx303_"+str(chr)+".cov", 'w')
    for a in open("temp_"+str(chr)+".bed", 'r') :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for region in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
            countmet += region[0]
            countunmet += region[1]
            uniquec += region[2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    ### B73 CHH % methylation
    print "### calculate CHH methylation rate"
    print "B73"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    
    for a in open("B73_all3/CHH_context_B73_all3_R1_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
       
    outfile = open("temp_chh_b73_"+str(chr)+".cov", 'w')
    for a in open("temp_"+str(chr)+".bed", 'r') :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for region in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
            countmet += region[0]
            countunmet += region[1]
            uniquec += region[2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "Mo17"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    
    for a in open("Mo17_all3/CHH_context_Mo17_all3_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
       
    outfile = open("temp_chh_mo17_"+str(chr)+".cov", 'w')
    for a in open("temp_"+str(chr)+".bed", 'r') :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for region in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
            countmet += region[0]
            countunmet += region[1]
            uniquec += region[2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "CML322"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    
    for a in open("CML322_all3/CHH_context_CML322_all3_bt202_sortedn_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
       
    outfile = open("temp_chh_cml322_"+str(chr)+".cov", 'w')
    for a in open("temp_"+str(chr)+".bed", 'r') :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for region in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
            countmet += region[0]
            countunmet += region[1]
            uniquec += region[2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "Oh43"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    
    for a in open("Oh43_all3/CHH_context_Oh43_all3_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
       
    outfile = open("temp_chh_oh43_"+str(chr)+".cov", 'w')
    for a in open("temp_"+str(chr)+".bed", 'r') :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for region in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
            countmet += region[0]
            countunmet += region[1]
            uniquec += region[2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    
    print "Tx303"
    large_list = []
    for x in range(0,310000000):
        large_list.append([0]*3)
    
    for a in open("Tx303_all3/CHH_context_Tx303_all3_bt202_chr"+str(chr)+".bismark.cov", 'r') :
        large_list[int(a.split("\t")[1])][0] += int(a.split("\t")[4])
        large_list[int(a.split("\t")[1])][1] += int(a.split("\t")[5].strip("\n"))
        large_list[int(a.split("\t")[1])][2] += 1
       
    outfile = open("temp_chh_tx303_"+str(chr)+".cov", 'w')
    for a in open("temp_"+str(chr)+".bed", 'r') :
        countmet = 0
        countunmet = 0
        uniquec = 0
        for region in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
            countmet += region[0]
            countunmet += region[1]
            uniquec += region[2]
        outfile.write("\t".join(a.split("\t")[:2])+"\t"+str(countmet)+"\t"+str(countunmet)+"\t"+str(uniquec)+"\n")
    outfile.close()
    del large_list
    """
    
    outfile = open("CpG_DMR_5genos_chr"+str(chr)+".txt", 'w')
    a = open("temp_"+str(chr)+".bed", 'r')
    b1 = open("temp_b73_"+str(chr)+".cov", 'r')
    b2 = open("temp_mo17_"+str(chr)+".cov", 'r')
    b3 = open("temp_oh43_"+str(chr)+".cov", 'r')
    b4 = open("temp_cml322_"+str(chr)+".cov", 'r')
    b5 = open("temp_tx303_"+str(chr)+".cov", 'r')
    c = open("temp_"+str(chr)+".allC", 'r')
    d = open("temp_bam_count_"+str(chr)+".txt", 'r')
    e = open("temp_meDIP_"+str(chr)+".txt", 'r')
    f1 = open("temp_chg_b73_"+str(chr)+".cov", 'r')
    f2 = open("temp_chg_mo17_"+str(chr)+".cov", 'r')
    f3 = open("temp_chg_oh43_"+str(chr)+".cov", 'r')
    f4 = open("temp_chg_cml322_"+str(chr)+".cov", 'r')
    f5 = open("temp_chg_tx303_"+str(chr)+".cov", 'r')
    g1 = open("temp_chh_b73_"+str(chr)+".cov", 'r')
    g2 = open("temp_chh_mo17_"+str(chr)+".cov", 'r')
    g3 = open("temp_chh_oh43_"+str(chr)+".cov", 'r')
    g4 = open("temp_chh_cml322_"+str(chr)+".cov", 'r')
    g5 = open("temp_chh_tx303_"+str(chr)+".cov", 'r')
    h = open("temp_alldmr_"+str(chr)+".txt", 'r')
    alines = a.readlines()
    b1lines = b1.readlines()
    b2lines = b2.readlines()
    b3lines = b3.readlines()
    b4lines = b4.readlines()
    b5lines = b5.readlines()
    clines = c.readlines()
    dlines = d.readlines()
    elines = e.readlines()
    f1lines = f1.readlines()
    f2lines = f2.readlines()
    f3lines = f3.readlines()
    f4lines = f4.readlines()
    f5lines = f5.readlines()
    g1lines = g1.readlines()
    g2lines = g2.readlines()
    g3lines = g3.readlines()
    g4lines = g4.readlines()
    g5lines = g5.readlines()
    hlines = h.readlines()
    
    print chr, len(alines), len(b1lines), len(b2lines), len(b3lines), len(b4lines), len(b5lines), len(clines), len(dlines), len(elines)
    outfile.write("chr\tstart\tend\tlength\tContributing genotypes\tB73 cov\tMo17 cov\tOh43 cov\tCML322 cov\tTx303 cov\tB73 meth%\tMo17 meth%\tOh43 meth%\tCML322 meth%\tTx303 meth%\tB73 CHG meth%\tMo17 CHG meth%\tOh43 CHG meth%\tCML322 CHG meth%\tTx303 CHG meth%\tB73 CHH meth%\tMo17 CHH meth%\tOh43 CHH meth%\tCML322 CHH meth%\tTx303 CHH meth%\tB73 C mapped\tMo17 C mapped\tOh43 C mapped\tCML322 C mapped\tTx303 C mapped\tNo. of probe +/-300bp\tB73 mean probe val\tMo17 mean probe val\tOh43 mean probe val\tCML322 mean probe val\tTx303 mean probe val\tTotal CpG\tTotal CHG\tTotal CHH\n")
    for aline in alines :
        coordinates = aline.split("\t")[:3]
        length = str(int(aline.split("\t")[2])-int(aline.split("\t")[1]))
        ## get dmr genos
        contrib_geno_list = []
        for hline in hlines :
            list1 = range(int(hline.split("\t")[2]),int(hline.split("\t")[3].strip("\n"))+1)
            list2 = range(int(aline.split("\t")[1]),int(aline.split("\t")[2])+1)
            if len(set(list1)&set(list2)) > 0 : 
                if hline.split("\t")[1] not in contrib_geno_list :
                    contrib_geno_list.append(hline.split("\t")[1])
        if len(contrib_geno_list) == 0 :
            contrib_geno = "None"
        else :
            contrib_geno = str(",".join(contrib_geno_list))
        b73_uniquenoofC = str(b1lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(b1lines[alines.index(aline)].split("\t")[2])+float(b1lines[alines.index(aline)].split("\t")[3])) == 0 :
            b73_metrate = "NA"
        else : 
            b73_metrate  = str((float(b1lines[alines.index(aline)].split("\t")[2]))/(float(b1lines[alines.index(aline)].split("\t")[2])+float(b1lines[alines.index(aline)].split("\t")[3])))
        mo17_uniquenoofC = str(b2lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(b2lines[alines.index(aline)].split("\t")[2])+float(b2lines[alines.index(aline)].split("\t")[3])) == 0 :
            mo17_metrate = "NA"
        else : 
            mo17_metrate  = str((float(b2lines[alines.index(aline)].split("\t")[2]))/(float(b2lines[alines.index(aline)].split("\t")[2])+float(b2lines[alines.index(aline)].split("\t")[3])))
        oh43_uniquenoofC = str(b3lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(b3lines[alines.index(aline)].split("\t")[2])+float(b3lines[alines.index(aline)].split("\t")[3])) == 0 :
            oh43_metrate = "NA"
        else : 
            oh43_metrate  = str((float(b3lines[alines.index(aline)].split("\t")[2]))/(float(b3lines[alines.index(aline)].split("\t")[2])+float(b3lines[alines.index(aline)].split("\t")[3])))
        cml322_uniquenoofC = str(b4lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(b4lines[alines.index(aline)].split("\t")[2])+float(b4lines[alines.index(aline)].split("\t")[3])) == 0 :
            cml322_metrate = "NA"
        else : 
            cml322_metrate  = str((float(b4lines[alines.index(aline)].split("\t")[2]))/(float(b4lines[alines.index(aline)].split("\t")[2])+float(b4lines[alines.index(aline)].split("\t")[3])))
        tx303_uniquenoofC = str(b5lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(b5lines[alines.index(aline)].split("\t")[2])+float(b5lines[alines.index(aline)].split("\t")[3])) == 0 :
            tx303_metrate = "NA"
        else : 
            tx303_metrate  = str((float(b5lines[alines.index(aline)].split("\t")[2]))/(float(b5lines[alines.index(aline)].split("\t")[2])+float(b5lines[alines.index(aline)].split("\t")[3])))
        
        ## CHG
        b73_chg_uniquenoofC = str(f1lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(f1lines[alines.index(aline)].split("\t")[2])+float(f1lines[alines.index(aline)].split("\t")[3])) == 0 :
            b73_chg_metrate = "NA"
        else : 
            b73_chg_metrate  = str((float(f1lines[alines.index(aline)].split("\t")[2]))/(float(f1lines[alines.index(aline)].split("\t")[2])+float(f1lines[alines.index(aline)].split("\t")[3])))
        mo17_chg_uniquenoofC = str(f2lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(f2lines[alines.index(aline)].split("\t")[2])+float(f2lines[alines.index(aline)].split("\t")[3])) == 0 :
            mo17_chg_metrate = "NA"
        else : 
            mo17_chg_metrate  = str((float(f2lines[alines.index(aline)].split("\t")[2]))/(float(f2lines[alines.index(aline)].split("\t")[2])+float(f2lines[alines.index(aline)].split("\t")[3])))
        oh43_chg_uniquenoofC = str(f3lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(f3lines[alines.index(aline)].split("\t")[2])+float(f3lines[alines.index(aline)].split("\t")[3])) == 0 :
            oh43_chg_metrate = "NA"
        else : 
            oh43_chg_metrate  = str((float(f3lines[alines.index(aline)].split("\t")[2]))/(float(f3lines[alines.index(aline)].split("\t")[2])+float(f3lines[alines.index(aline)].split("\t")[3])))
        cml322_chg_uniquenoofC = str(f4lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(f4lines[alines.index(aline)].split("\t")[2])+float(f4lines[alines.index(aline)].split("\t")[3])) == 0 :
            cml322_chg_metrate = "NA"
        else : 
            cml322_chg_metrate  = str((float(f4lines[alines.index(aline)].split("\t")[2]))/(float(f4lines[alines.index(aline)].split("\t")[2])+float(f4lines[alines.index(aline)].split("\t")[3])))
        tx303_chg_uniquenoofC = str(f5lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(f5lines[alines.index(aline)].split("\t")[2])+float(f5lines[alines.index(aline)].split("\t")[3])) == 0 :
            tx303_chg_metrate = "NA"
        else : 
            tx303_chg_metrate  = str((float(f5lines[alines.index(aline)].split("\t")[2]))/(float(f5lines[alines.index(aline)].split("\t")[2])+float(f5lines[alines.index(aline)].split("\t")[3])))
        
        ## CHH
        b73_chh_uniquenoofC = str(g1lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(g1lines[alines.index(aline)].split("\t")[2])+float(g1lines[alines.index(aline)].split("\t")[3])) == 0 :
            b73_chh_metrate = "NA"
        else : 
            b73_chh_metrate  = str((float(g1lines[alines.index(aline)].split("\t")[2]))/(float(g1lines[alines.index(aline)].split("\t")[2])+float(g1lines[alines.index(aline)].split("\t")[3])))
        mo17_chh_uniquenoofC = str(g2lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(g2lines[alines.index(aline)].split("\t")[2])+float(g2lines[alines.index(aline)].split("\t")[3])) == 0 :
            mo17_chh_metrate = "NA"
        else : 
            mo17_chh_metrate  = str((float(g2lines[alines.index(aline)].split("\t")[2]))/(float(g2lines[alines.index(aline)].split("\t")[2])+float(g2lines[alines.index(aline)].split("\t")[3])))
        oh43_chh_uniquenoofC = str(g3lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(g3lines[alines.index(aline)].split("\t")[2])+float(g3lines[alines.index(aline)].split("\t")[3])) == 0 :
            oh43_chh_metrate = "NA"
        else : 
            oh43_chh_metrate  = str((float(g3lines[alines.index(aline)].split("\t")[2]))/(float(g3lines[alines.index(aline)].split("\t")[2])+float(g3lines[alines.index(aline)].split("\t")[3])))
        cml322_chh_uniquenoofC = str(g4lines[alines.index(aline)].split("\t")[4].strip("\n"))
        if (float(g4lines[alines.index(aline)].split("\t")[2])+float(g4lines[alines.index(aline)].split("\t")[3])) == 0 :
            cml322_chh_metrate = "NA"
        else : 
            cml322_chh_metrate  = str((float(g4lines[alines.index(aline)].split("\t")[2]))/(float(g4lines[alines.index(aline)].split("\t")[2])+float(g4lines[alines.index(aline)].split("\t")[3])))
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
        oh43_cov = str(dlines[alines.index(aline)].split("\t")[7])
        cml322_cov = str(dlines[alines.index(aline)].split("\t")[6])
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
        outfile.write("\t".join(coordinates)+"\t"+length+"\t"+str(contrib_geno)+"\t"+b73_cov+"\t"+mo17_cov+"\t"+oh43_cov+"\t"+cml322_cov+"\t"+tx303_cov+"\t"+b73_metrate+"\t"+mo17_metrate+"\t"+oh43_metrate+"\t"+cml322_metrate+"\t"+tx303_metrate+"\t"+b73_chg_metrate+"\t"+mo17_chg_metrate+"\t"+oh43_chg_metrate+"\t"+cml322_chg_metrate+"\t"+tx303_chg_metrate+"\t"+b73_chh_metrate+"\t"+mo17_chh_metrate+"\t"+oh43_chh_metrate+"\t"+cml322_chh_metrate+"\t"+tx303_chh_metrate+"\t"+b73_uniquenoofC+"\t"+mo17_uniquenoofC+"\t"+oh43_uniquenoofC+"\t"+cml322_uniquenoofC+"\t"+tx303_uniquenoofC+"\t"+probe_count+"\t"+b73_probe_val+"\t"+mo17_probe_val+"\t"+oh43_probe_val+"\t"+cml322_probe_val+"\t"+tx303_probe_val+"\t"+totalCpG+"\t"+totalCHG+"\t"+totalCHH+"\n")
        
    outfile.close()
    #os.system("rm -Rf temp*"+str(chr)+".*")
"""
jobs = []
for chr in range(1,2):
    s1 = multiprocessing.Process(target=merge, args=(chr, ))
    jobs.append(s1)
    s1.start()
[x.join() for x in jobs]
sys.exit()
outfile = open("CpG_DMR_5genos_merged.txt", 'w')
for chr in range(1,11) :
    print chr
    if chr == 1 :
        for a in open("CpG_DMR_5genos_chr"+str(chr)+".txt", 'r') :
            outfile.write(a)
    else :
        a = open("CpG_DMR_5genos_chr"+str(chr)+".txt", 'r')
        alines = a.readlines()
        for aline in alines[1:] :
            outfile.write(aline)
outfile.close()  
sys.exit()     