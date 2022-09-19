import multiprocessing
import os,sys
import datetime
import subprocess
import optparse
def bismark(name) :
    print "running bismark/methylation extractor"
    print subprocess.Popen("module load bowtie/2.1.0", stdout=subprocess.PIPE, shell=True, executable="/bin/bash").stdout.read()
    print subprocess.Popen("module load bismark/0.7.12", stdout=subprocess.PIPE, shell=True, executable="/bin/bash").stdout.read()
    os.system("/opt/apps/bismark/0.7.12/methylation_extractor  -p --no_overlap --comprehensive --counts --report "+str(name)+".sam")
    
def split(chr, type, name):
    print "split chromosomes", chr, type, name
    outfile = open(str(type)+"_context_"+str(name)+"_chr"+str(chr)+".out", 'w')
    infile = open(str(type)+"_context_"+str(name)+".txt", 'r')
    lines = infile.readlines()
    for a in lines :
        if len(a.split("\t")) > 2 and str(a.split("\t")[2]) == "chr"+str(chr) :  
            outfile.write(a) 
    outfile.close()
    return

def split_by_chr(chr, type, name) :
    print datetime.datetime.now()
    print "split chromosomes", chr, type, name
    list = []
    #infile = open(str(type)+"_context_"+str(name)+".txt", 'r')
    #lines = infile.readlines()
    #for a in lines :
    for a in open(str(type)+"_context_"+str(name)+".txt", 'r') :
        if len(a.split("\t")) > 2 and str(a.split("\t")[2]) == "chr"+str(chr) :  
            list.append(a)
    print datetime.datetime.now()
    outfile = open(str(type)+"_context_"+str(name)+"_chr"+str(chr)+".out", 'w')
    for line in list :
        outfile.write(line)
    outfile.close()
    
def combine_type(chr, name):
    print "combine contexts", chr, name
    if os.path.isfile("All_context_"+str(name)+"_chr"+str(chr)+"_multithread.out") :
        print "file exists, remove"
        os.system("rm -Rf All_context_"+str(name)+"_chr"+str(chr)+"_multithread.out")
    list = ['CpG', 'CHH', 'CHG']
    for type in list :
        os.system("cat "+str(type)+"_context_"+str(name)+"_chr"+str(chr)+"_multithread.out >> All_context_"+str(name)+"_chr"+str(chr)+"_multithread.out")
    return
    
def bis2bed(chr, type, name):
    print "bismark to bedGraph", chr, type, name
    os.system("/work/02114/wonaya/software/bismark2bedGraph --buffer_size 50% --counts --CX_context -o "+str(type)+"_context_"+str(name)+"_chr"+str(chr)+".bedGraph "+str(type)+"_context_"+str(name)+"_chr"+str(chr)+".out")
    return

def pre_methylKit(chr, type, name) :
    print "methylKit prep", chr, type, name
    pos = []
    neg = []
    dict = {}
    ## save strand status
    for a in open(str(type)+"_context_"+str(name)+"_chr"+str(chr)+".out", 'r') : # example: CpG_context_B73_ACAGTG_merged.sorted_chr10.out bsseq/methylKit/swDMR
        dict[a.split("\t")[3]] = a.split("\t")[1]
    outfile = open(str(type)+"_context_"+str(name)+"_chr"+str(chr)+".methylKit", 'w')
    for b in open(str(type)+"_context_"+str(name)+"_chr"+str(chr)+".bedGraph", 'r') :
        if  dict[str(int(b.split("\t")[1])+1)] == "+" :
            strand = "F"
        elif  dict[str(int(b.split("\t")[1])+1)] == "-" :
            strand = "R"
        total = float(int(b.split("\t")[4])+int(b.split("\t")[5].strip("\n")))
        freqc = float(int(b.split("\t")[4]))
        freqt = float(int(b.split("\t")[5].strip("\n")))
        percc = freqc/total*100
        perct = freqt/total*100
        outfile.write(str(b.split("\t")[0])+"."+str(b.split("\t")[1]))
        outfile.write("\t")
        outfile.write(str(b.split("\t")[0])+"\t"+str(b.split("\t")[1]))
        outfile.write("\t")
        outfile.write(str(strand)+"\t"+str(int(total))+"\t"+str(percc)+"\t"+str(perct)+"\n")
    outfile.close()

def merge_strands(chr, type, name) :
    list_unmerge = []
    for a in open(str(type)+"_context_"+str(name)+"_chr"+str(chr)+".bedGraph", 'r') :
        list_unmerge.append(int(a.split("\t")[1]))
    b = open(str(type)+"_context_"+str(name)+"_chr"+str(chr)+".bedGraph", 'r') 
    blines = b.readlines()
    end = int(blines[-1].split("\t")[1])
    list_merge = [0]
    for coord in list_unmerge :
        if coord != end :
            if list_unmerge[list_unmerge.index(coord)+1] == coord+1 and list_merge[-1] != coord-1 :
                list_unmerge.remove(list_unmerge[list_unmerge.index(coord)+1])
                list_merge.append(coord)
    outfile = open(str(type)+"_context_"+str(name)+"_chr"+str(chr)+"_merged.bedGraph", 'w')
    for bline in blines :
        if int(bline.split("\t")[1]) != end :
            if int(bline.split("\t")[1]) in list_merge :
                met = int(bline.split("\t")[4])+int(blines[blines.index(bline)+1].split("\t")[4]) 
                unmet = int(bline.split("\t")[5].strip("\n"))+int(blines[blines.index(bline)+1].split("\t")[5].strip("\n")) 
                outfile.write(str(bline.split("\t")[0])+"\t"+str(bline.split("\t")[1])+"\t"+str(bline.split("\t")[2])+"\t"+str(float(met)/(float(met)+float(unmet))*100)+"\t"+str(met)+"\t"+str(unmet)+"\n")
                blines.remove(blines[blines.index(bline)+1])
            else :
                outfile.write(bline)
    outfile.close()    
            
def methylkit(chr, type, name1, name2) :
    print "methylKit run", chr, type, name1, name2
    os.system("/opt/apps/R/2.15.3/bin/Rscript /work/02114/wonaya/scripts/r_test2.R "+str(type)+"_context_"+str(name1)+"_chr"+str(chr)+".methylKit "+str(type)+"_context_"+str(name2)+"_chr"+str(chr)+".methylKit "+str(type)+" "+str(type)+"_context_chr"+str(chr)+".out")

def merge_chrom(type) :
    print "merge chromsomes", type
    outfile = open(str(type)+"_context_merged_DMR.out", 'w')
    outfile.write("chrom\tstart\tend\twidth\tmean.meth.diff\tnum."+str(type)+"\tnum.DMCs\tDMR.pvalue\tDMR.qvalue\n")
    for x in range(1,11) :
        infile = open(str(type)+"_context_chr"+str(x)+".txt", 'r') 
        lines = infile.readlines()
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

def wmet(chr, type, name) :
    print "wMet calculation", name, type, chr
    dmr_loc = {}
    dmr_loc_start_val = []
    file = open(str(type)+"_context_chr"+str(chr)+".txt", 'r')
    lines = file.readlines()   
    for c in lines[1:] : 
        dmr_loc[int(c.split("\t")[2])] = int(c.split("\t")[3])
        dmr_loc_start_val.append(int(c.split("\t")[2]))     
    dmr_loc_start_val.sort()   
    
    ## calculate weighted methylation
    a = open(str(type)+"_context_"+str(name)+"_chr"+str(chr)+"_multithread.out", 'r') 
    alines = a.readlines()
    outfile = open(str(type)+"_context_"+str(name)+"_chr"+str(chr)+"_wMet.out", 'w')
    outfile.write("chrom\tstart\ttotal C\ttotal Read\ttotal valid C\ttotal valid reads\tFrag\tMean\tWeighted\n") 
    for dmr in dmr_loc_start_val :
        by_dmr = []
        pos = []
        for aline in alines :
            if int(aline.split("\t")[3]) >= int(dmr) and int(aline.split("\t")[3]) <= int(dmr_loc[dmr]) :
                by_dmr.append(aline)
                if int(aline.split("\t")[3]) not in pos :
                    pos.append(int(aline.split("\t")[3]))
        pos.sort()
        length = int(dmr_loc[dmr])-int(dmr)
        total_c = len(pos)
        total_read = len(by_dmr)
        val_met = 0
        val_met_val = 0
        weighted_top = 0
        for position in pos :
            count_pos = 0
            count_all = 0
            for indiv_read in by_dmr :
                if int(indiv_read.split("\t")[3]) == position :
                    count_all += 1
                    if indiv_read.split("\t")[1] == "+" :
                        count_pos += 1
            if float(count_pos)/float(count_all)  >= 0.126 :
                val_met += 1
                val_met_val += float(count_pos)/float(count_all)
                weighted_top += float(count_pos)
        frag = float(val_met)/float(total_c)
        mean = float(val_met_val)/float(total_c)
        weighted = float(weighted_top)/total_read
    
        ## write
        outfile.write("chr"+str(chr)+"\t")
        outfile.write(str(dmr)+"\t")
        outfile.write(str(dmr_loc[dmr])+"\t")
        outfile.write(str(total_c)+"\t")
        outfile.write(str(total_read)+"\t")
        outfile.write(str(val_met)+"\t")
        outfile.write(str(weighted_top)+"\t")
        outfile.write(str(frag)+"\t"+str(mean)+"\t"+str(weighted)+"\n")
    outfile.close()

def merge_wmet(type, name) :
    print "merge wmet", type, name
    outfile = open(str(type)+"_context_"+str(name)+"_merged_wMet.out", 'w') 
    outfile.write("chrom\tstart\ttotal C\ttotal Read\ttotal valid C\ttotal valid reads\tFrag\tMean\tWeighted\n")
    for x in range(1,11) :
        infile = open(str(type)+"_context_"+str(name)+"_chr"+str(x)+"_wMet.out", 'r')
        lines = infile.readlines()
        for line in lines[1:] :
            outfile.write(line)
    outfile.close()
            
def tile(name, chr, index) :#sort methylation extracted file by coordinates
    #read in files CpG, CHH, CHG 
    if not os.path.isfile("CpG_context_"+str(name)+"_chr"+str(chr)+"_sorted.out") :
        os.system("sort -nk4 CpG_context_"+str(name)+"_chr"+str(chr)+".out > CpG_context_"+str(name)+"_chr"+str(chr)+"_sorted.out")
    print "tile : reading CpG files"
    
    #divide region by 1M Coordinates (approx 100 parallel jobs)
    x = int(index)
    
    cpgfile = open("CpG_context_"+str(name)+"_chr"+str(chr)+"_sorted.out", 'r')
    cpglines = cpgfile.readlines()
    cpg_split_lines = []
    for cpgline in cpglines :
        if int(cpgline.split("\t")[3]) >= x*1000000 and int(cpgline.split("\t")[3]) < (x+1)*1000000 :
            cpg_split_lines.append(cpgline)
    del cpglines
    
    if not os.path.isfile("CHH_context_"+str(name)+"_chr"+str(chr)+"_sorted.out") :
        os.system("sort -nk4 CHH_context_"+str(name)+"_chr"+str(chr)+".out > CHH_context_"+str(name)+"_chr"+str(chr)+"_sorted.out")
    print "tile : reading CHH files"
    
    chhfile = open("CHH_context_"+str(name)+"_chr"+str(chr)+"_sorted.out", 'r')
    chhlines = chhfile.readlines()
    chh_split_lines = []
    for chhline in chhlines :
        if int(chhline.split("\t")[3]) >= x*1000000 and int(chhline.split("\t")[3]) < (x+1)*1000000 :
            chh_split_lines.append(chhline)
    del chhlines
    
    if not os.path.isfile("CHG_context_"+str(name)+"_chr"+str(chr)+"_sorted.out") :
        os.system("sort -nk4 CHG_context_"+str(name)+"_chr"+str(chr)+".out > CHG_context_"+str(name)+"_chr"+str(chr)+"_sorted.out")
    print "tile : reading CHG files"
    chgfile = open("CHG_context_"+str(name)+"_chr"+str(chr)+"_sorted.out", 'r')
    chglines = chgfile.readlines()
    chg_split_lines = []
    for chgline in chglines :
        if int(chgline.split("\t")[3]) >= x*1000000 and int(chgline.split("\t")[3]) < (x+1)*1000000 :
            chg_split_lines.append(chgline)
    del chglines
    print datetime.datetime.now()

    #for each region, start from 0 or n*1000,000 in increment of 100 to capture methylations
    outfile = open(str(name)+"_chr"+str(chr)+"_tile_"+str(x)+".out", 'w')
    for y in range(x*10000, (x+1)*10000) :
        metCpG = 0
        totCpG = 0
        for cpg_split_line in cpg_split_lines :
            if int(cpg_split_line.split("\t")[3]) >= y*100 and int(cpg_split_line.split("\t")[3]) < (y+1)*100 :
                if cpg_split_line.split("\t")[1] == "+" :
                    metCpG += 1
                    totCpG += 1
                elif cpg_split_line.split("\t")[1] == "-" :
                    totCpG += 1
        metCHH = 0 
        totCHH = 0
        for chh_split_line in chh_split_lines :
            if int(chh_split_line.split("\t")[3]) >= y*100 and int(chh_split_line.split("\t")[3]) < (y+1)*100 :
                if chh_split_line.split("\t")[1] == "+" :
                    metCHH += 1
                    totCHH += 1
                elif chh_split_line.split("\t")[1] == "-" :
                    totCHH += 1
    
        metCHG = 0 
        totCHG = 0
        for chg_split_line in chg_split_lines :
            if int(chg_split_line.split("\t")[3]) >= y*100 and int(chg_split_line.split("\t")[3]) < (y+1)*100 :
                if chg_split_line.split("\t")[1] == "+" :
                    metCHG += 1
                    totCHG += 1
                elif chg_split_line.split("\t")[1] == "-" :
                    totCHG += 1
    
        if totCpG != 0 and totCHH != 0 and totCHG != 0 :
            outfile.write(str(y*100)+"\t"+str((y+1)*100)+"\t"+str(metCpG)+"\t"+str(totCpG)+"\t"+str(metCHH)+"\t"+str(totCHH)+"\t"+str(metCHG)+"\t"+str(totCHG)+"\t"+str(metCpG+metCHH+metCHG)+"\t"+str(totCpG+totCHH+totCHG)+"\n")
    outfile.close()

if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("--run", default=None, help="select from: split, bismark, combine_all")
    parser.add_option("--name", default=None, help="select from: B73_all3, Mo17_all3, Oh43_all3, CML322_all3, Tx303_all3")
    parser.add_option("--specie", default="maize", help="select from: arabidopsis, maize")
    parser.add_option("--met_type", default=None, help="select from: CpG, CHH, CHG, All")
    
    options, args = parser.parse_args()
    #met_type = ['CHG', 'CpG', 'CHH', 'All']
    list_of_defs = ["bismark", "split", "split_argv", "split_chr", "combine_type"]
    if sys.argv[1] not in list_of_defs :
        print "please either select bismark, split, split_argv, split_chr, combine_type"
        sys.exit()
        
    if sys.argv[1] == "bismark" :
        jobs = []
        names = ['Mo17_all3', 'Oh43_all3']
        for name in mames :
            o = multiprocessing.Process(target=bismark, args=(name, ))
            jobs.append(o)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
                
    elif sys.argv[1] == "split" :
        jobs = []
        for name in names :
            for met in met_type :
                for chr in range(1,11): ## number of jobs
                #for chr in [4,2,7,10] :
                    p = multiprocessing.Process(target=split, args=(chr, met, name,))
                    jobs.append(p)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
    
    elif sys.argv[1] == "split_argv" :
        jobs = []
        for chr in range(1,11): ## number of jobs
            p = multiprocessing.Process(target=split_by_chr, args=(chr, sys.argv[2], sys.argv[3],))
            jobs.append(p)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
    
    elif sys.argv[1] == "split_chr" :
        split_by_chr(sys.argv[2], sys.argv[3], sys.argv[4])
        
    elif sys.argv[1] == "combinetype" :
        for chr in range(10,11):
            for name in names :
                q = multiprocessing.Process(target=combine_type, args=(chr, name, ))
                jobs.append(q)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
    
    elif sys.argv[1] == "bis2bed" :
        jobs = []
        met_type = ['All']
        for name in names :
            for met in met_type :
                #for chr in range(1,11) :
                for chr in [1,4,10] :
                    r = multiprocessing.Process(target=bis2bed, args=(chr, met, name, ))
                    jobs.append(r)
        [x.start() for x in jobs]
        [x.join()  for x in jobs]
        
    elif sys.argv[1] == "prep" :
        jobs = []
        for name in names :
            for met in met_type :
                for chr in range(1,11) :
                    s = multiprocessing.Process(target=pre_methylKit, args=(chr, met, name, ))
                    jobs.append(s)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
        
    elif sys.argv[1] == "methylKit" :
        print datetime.datetime.now()
        jobs = []
        for met in met_type :
            for chr in range(10,11) :
                t = multiprocessing.Process(target=methylkit, args=(chr, met, names[0], names[1], ))
                jobs.append(t)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
        print datetime.datetime.now()
        
    elif sys.argv[1] == "merge_strands" :
        jobs = []
        for name in names :
            u = multiprocessing.Process(target=merge_strands, args=("10", "CpG", name,))
            jobs.append(u)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
    
    elif sys.argv[1] == "merge_chrom" :
        jobs = []
        for met in met_type :
            v = multiprocessing.Process(target=merge_chrom, args=(met, ))
            jobs.append(v)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
                    
    elif sys.argv[1] == "wmet" :
        jobs = []
        for name in names :
            for met in met_type :
                for chr in range(1,11) :
                    w = multiprocessing.Process(target=wmet, args=(chr, met, name, ))
                    jobs.append(w)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
    elif sys.argv[1] == "merge_wmet" :
        jobs = []
        for name in names :
            for met in met_type :
                z = multiprocessing.Process(target=merge_wmet, args=(met,name,))
                jobs.append(z)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
    
    elif sys.argv[1] == "tile" :
        names = ['Mo17_all3']
        jobs = []
        for name in names :
            for chr in range(10,11) :
                for index in range(0,1) :
                    y = multiprocessing.Process(target=tile, args=(name, chr, index,))
                    jobs.append(y)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
    elif sys.argv[1] == "tile_argv" :
        tile(sys.argv[2], sys.argv[3], sys.argv[4])       