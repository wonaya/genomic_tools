import os,sys
import datetime
import random
import multiprocessing
import subprocess
import optparse
import numpy as np

def tile(name, chr, window) :
    ## define genotype and chromosome
     
    ## make array of 8 columns and put 0
    ## CpG/TotC/CHG/TotC/CHH/TotC/all/TotC forward (first 8 columns) and reverse
    
    print "make array", "chr"+str(chr), datetime.datetime.now()
    large_list = np.zeros((5000000, 16))
    
    ### CpG
    ### forward strand
    print "read CpG OT", "chr"+str(chr), datetime.datetime.now()
    for CpG_OT in open("CpG_OT_"+str(name)+".txt", 'r') :
        if CpG_OT[0:7] != "Bismark" :
            if CpG_OT.split("\t")[2] == str(chr) :
                if CpG_OT.split("\t")[1] == "+" :
                    large_list[int(float(CpG_OT.split("\t")[3])/int(window))][0] += 1 
                    large_list[int(float(CpG_OT.split("\t")[3])/int(window))][1] += 1
                elif CpG_OT.split("\t")[1] == "-" :
                    large_list[int(float(CpG_OT.split("\t")[3])/int(window))][1] += 1
    del CpG_OT
    """print "read CpG CTOT", "chr"+str(chr), datetime.datetime.now()
    for CpG_CTOT in open("CpG_CTOT_"+str(name)+".txt", 'r') :
        if CpG_CTOT[0:7] != "Bismark" :
            if CpG_CTOT.split("\t")[2] == "chr"+str(chr) :
                if CpG_CTOT.split("\t")[1] == "+" :
                    large_list[int(float(CpG_CTOT.split("\t")[3])/int(window))][0] += 1 
                    large_list[int(float(CpG_CTOT.split("\t")[3])/int(window))][1] += 1
                elif CpG_CTOT.split("\t")[1] == "-" :
                    large_list[int(float(CpG_CTOT.split("\t")[3])/int(window))][1] += 1
    del CpG_CTOT     
    """
    ### reverse strand 
    print "read CpG OB","chr"+str(chr),datetime.datetime.now()
    for CpG_OB in open("CpG_OB_"+str(name)+".txt", 'r') :
        if CpG_OB[0:7] != "Bismark" :
            if CpG_OB.split("\t")[2] == str(chr) :
                if CpG_OB.split("\t")[1] == "+" :
                    large_list[int(float(CpG_OB.split("\t")[3])/int(window))][8] += 1 
                    large_list[int(float(CpG_OB.split("\t")[3])/int(window))][9] += 1
                elif CpG_OB.split("\t")[1] == "-" :
                    large_list[int(float(CpG_OB.split("\t")[3])/int(window))][9] += 1
    del CpG_OB
    """
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
    """
    ### CHG
    ### forward strand
    print "read CHG OT","chr"+str(chr),datetime.datetime.now()
    for CHG_OT in open("CHG_OT_"+str(name)+".txt", 'r') :
        if CHG_OT[0:7] != "Bismark" :
            if CHG_OT.split("\t")[2] == str(chr) :
                print CHG_OT
                if CHG_OT.split("\t")[1] == "+" :
                    large_list[int(float(CHG_OT.split("\t")[3])/int(window))][2] += 1 
                    large_list[int(float(CHG_OT.split("\t")[3])/int(window))][3] += 1
                elif CHG_OT.split("\t")[1] == "-" :
                    large_list[int(float(CHG_OT.split("\t")[3])/int(window))][3] += 1
    del CHG_OT
    """print "read CHG CTOT","chr"+str(chr),datetime.datetime.now()
    for CHG_CTOT in open("CHG_CTOT_"+str(name)+".txt", 'r') :
        if CHG_CTOT[0:7] != "Bismark" :
            if CHG_CTOT.split("\t")[2] == "chr"+str(chr) :
                if CHG_CTOT.split("\t")[1] == "+" :
                    large_list[int(float(CHG_CTOT.split("\t")[3])/int(window))][2] += 1 
                    large_list[int(float(CHG_CTOT.split("\t")[3])/int(window))][3] += 1
                elif CHG_CTOT.split("\t")[1] == "-" :
                    large_list[int(float(CHG_CTOT.split("\t")[3])/int(window))][3] += 1
    del CHG_CTOT    
    """### reverse strand 
    print "read CHG OB","chr"+str(chr),datetime.datetime.now()
    for CHG_OB in open("CHG_OB_"+str(name)+".txt", 'r') :
        if CHG_OB[0:7] != "Bismark" :
            if CHG_OB.split("\t")[2] == str(chr) :
                if CHG_OB.split("\t")[1] == "+" :
                    large_list[int(float(CHG_OB.split("\t")[3])/int(window))][10] += 1 
                    large_list[int(float(CHG_OB.split("\t")[3])/int(window))][11] += 1
                elif CHG_OB.split("\t")[1] == "-" :
                    large_list[int(float(CHG_OB.split("\t")[3])/int(window))][11] += 1
    del CHG_OB
    """print "read CHG CTOB","chr"+str(chr),datetime.datetime.now()
    for CHG_CTOB in open("CHG_CTOB_"+str(name)+".txt", 'r') :
        if CHG_CTOB[0:7] != "Bismark" :
            if CHG_CTOB.split("\t")[2] == "chr"+str(chr) :
                if CHG_CTOB.split("\t")[1] == "+" :
                    large_list[int(float(CHG_CTOB.split("\t")[3])/int(window))][10] += 1 
                    large_list[int(float(CHG_CTOB.split("\t")[3])/int(window))][11] += 1
                elif CHG_CTOB.split("\t")[1] == "-" :
                    large_list[int(float(CHG_CTOB.split("\t")[3])/int(window))][11] += 1
    del CHG_CTOB
    """
    ### CHH
    ### forward strand
    print "read CHH OT","chr"+str(chr),datetime.datetime.now()
    for CHH_OT in open("CHH_OT_"+str(name)+".txt", 'r') :
        if CHH_OT[0:7] != "Bismark" :
            if CHH_OT.split("\t")[2] == str(chr) :
                if CHH_OT.split("\t")[1] == "+" :
                    large_list[int(float(CHH_OT.split("\t")[3])/int(window))][4] += 1 
                    large_list[int(float(CHH_OT.split("\t")[3])/int(window))][5] += 1
                elif CHH_OT.split("\t")[1] == "-" :
                    large_list[int(float(CHH_OT.split("\t")[3])/int(window))][5] += 1
    del CHH_OT
    """print "read CHH CTOT","chr"+str(chr),datetime.datetime.now()
    for CHH_CTOT in open("CHH_CTOT_"+str(name)+".txt", 'r') :
        if CHH_CTOT[0:7] != "Bismark" :
            if CHH_CTOT.split("\t")[2] == "chr"+str(chr) :
                if CHH_CTOT.split("\t")[1] == "+" :
                    large_list[int(float(CHH_CTOT.split("\t")[3])/int(window))][4] += 1 
                    large_list[int(float(CHH_CTOT.split("\t")[3])/int(window))][5] += 1
                elif CHH_CTOT.split("\t")[1] == "-" :
                    large_list[int(float(CHH_CTOT.split("\t")[3])/int(window))][5] += 1
    del CHH_CTOT     
    """### reverse strand 
    
    print "read CHH OB","chr"+str(chr),datetime.datetime.now()
    for CHH_OB in open("CHH_OB_"+str(name)+".txt", 'r') :
        if CHH_OB[0:7] != "Bismark" :
            if CHH_OB.split("\t")[2] == str(chr) :
                if CHH_OB.split("\t")[1] == "+" :
                    large_list[int(float(CHH_OB.split("\t")[3])/int(window))][12] += 1 
                    large_list[int(float(CHH_OB.split("\t")[3])/int(window))][13] += 1
                elif CHH_OB.split("\t")[1] == "-" :
                    large_list[int(float(CHH_OB.split("\t")[3])/int(window))][13] += 1
    del CHH_OB
    
    """print "read CHH CTOB","chr"+str(chr),datetime.datetime.now()
    for CHH_CTOB in open("CHH_CTOB_"+str(name)+".txt", 'r') :
        if CHH_CTOB[0:7] != "Bismark" :
            if CHH_CTOB.split("\t")[2] == "chr"+str(chr) :
                if CHH_CTOB.split("\t")[1] == "+" :
                    large_list[int(float(CHH_CTOB.split("\t")[3])/int(window))][12] += 1 
                    large_list[int(float(CHH_CTOB.split("\t")[3])/int(window))][13] += 1
                elif CHH_CTOB.split("\t")[1] == "-" :
                    large_list[int(float(CHH_CTOB.split("\t")[3])/int(window))][13] += 1
    del CHH_CTOB
    """
    ## 0: CpG met top 1: CpG tot top 2: CHG met top 3: CHG tot top 8: CpG met bot 9: CpG tot bot 10: CHG met bot 11: CHG tot bot
    print "all metC, metC", "chr"+str(chr),datetime.datetime.now()
    large_list_2 = large_list.tolist()
    del large_list
    for list in large_list_2 :
        list[6] += list[0]+list[2]+list[4] ## top tot met
        list[7] += list[1]+list[3]+list[5] ## top total C
        list[14] += list[8]+list[10]+list[12] ## bot tot met
        list[15] += list[9]+list[11]+list[13] ## bot total C
    print "write to file", "chr"+str(chr),datetime.datetime.now()
    outfile = open(str(name)+"_tile_chr"+str(chr)+"_"+str(window)+"bp.txt", 'w')
    outfile.write("Start\tEnd\tChrom\tF CpG metC\tF CpG C\tF CHG metC\tF CHG C\tF CHH metC\tF CHH C\tF all metC\tF all C\tR CpG metC\tR CpG C\tR CHG metC\tR CHG C\tR CHH metC\tR CHH C\tR all metC\tR all C\n")
    x = 0
    for list in large_list_2 :
        if list[15] != 0 or list[7] != 0 :
            outfile.write(str(x*int(window))+"\t")
            outfile.write(str((x+1)*int(window)-1)+"\t")
            outfile.write("chr"+str(chr)+"\t")
            list = map(int, list)
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

if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("--run", default=None, help="tile merge bedgraph circos")
    parser.add_option("--name", default=None, help="select from: B73_all3, Mo17_all3, Oh43_all3, CML322_all3, Tx303_all3 (basically name of bam used for the BME)")
    parser.add_option("--window", default=100, help="specify window size 100, 50, 10000")
    #parser.add_option("--type", default=None, help="use with bedGraph CpG, CHH, CHG or all")
    #parser.add_option("--direction", default=None, help="use with bedGraph reverse, forward or both")
    #parser.add_option("--chrom", default=None, help="chromosome 1 to 10")
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
    elif options.run == "split_chr" :
        split_by_chr(options.chrom, options.type, options.name)
        #os.system("/work/02114/wonaya/software/bismark_v0.10.0/bismark2bedGraph --counts --CX_context -o "+str(options.type)+"_context_"+str(options.name)+"_chr"+str(options.chrom)+".bedGraph "+str(options.type)+"_context_"+str(options.name)+"_chr"+str(options.chrom)+".out")
    
    elif options.run == "test" :
        pre_methylKit(9, "CpG", options.name)
        
    
    elif options.run == "prep" :
        jobs = []
        if options.specie == "maize" :
            chrmax = 11
        for name in options.name.split(",") :
            for met in options.type.split(",") :
                for chr in range(1,chrmax) :
                    s = multiprocessing.Process(target=pre_methylKit, args=(chr, met, name, ))
                    jobs.append(s)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
