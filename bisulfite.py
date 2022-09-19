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

def max_coord(name, chr) :
    # determine max coordinate of CpG, CHH, CHG
    print "sorting file"
    if not os.path.isfile("CpG_context_"+str(name)+"_chr"+str(chr)+"_sorted.out") :
        os.system("sort -nk4 CpG_context_"+str(name)+"_chr"+str(chr)+".out > CpG_context_"+str(name)+"_chr"+str(chr)+"_sorted.out")
    if not os.path.isfile("CHH_context_"+str(name)+"_chr"+str(chr)+"_sorted.out") :
        os.system("sort -nk4 CHH_context_"+str(name)+"_chr"+str(chr)+".out > CHH_context_"+str(name)+"_chr"+str(chr)+"_sorted.out")
    if not os.path.isfile("CHG_context_"+str(name)+"_chr"+str(chr)+"_sorted.out") :
        os.system("sort -nk4 CHG_context_"+str(name)+"_chr"+str(chr)+".out > CHG_context_"+str(name)+"_chr"+str(chr)+"_sorted.out")
    last = []
    cpgfile = open("CpG_context_"+str(name)+"_chr"+str(chr)+"_sorted.out", 'r')
    cpglines = cpgfile.readlines()
    last.append(int(cpglines[-1].split("\t")[3]))
    del cpglines
    chgfile = open("CHG_context_"+str(name)+"_chr"+str(chr)+"_sorted.out", 'r')
    chglines = chgfile.readlines()
    last.append(int(chglines[-1].split("\t")[3]))
    del chglines
    chhfile = open("CHH_context_"+str(name)+"_chr"+str(chr)+"_sorted.out", 'r')
    chhlines = chhfile.readlines()
    last.append(int(chhlines[-1].split("\t")[3]))
    del chhlines
    print "done sorting file"
    return max(last)   
    
def tile(name, chr, index) :#sort methylation extracted file by coordinates
    #read in files CpG, CHH, CHG 
    print "tile : reading CpG files", index
    
    #divide region by 1M Coordinates (approx 100 parallel jobs)
    x = int(index)
    
    #cpgfile = open("CpG_context_"+str(name)+"_chr"+str(chr)+"_sorted.out", 'r')
    #cpglines = cpgfile.readlines()
    cpg_split_lines = []
    for cpgline in open("CpG_context_"+str(name)+"_chr"+str(chr)+"_sorted.out", 'r')
        print cpgline
        sys.exit()
        if int(cpgline.split("\t")[3]) >= x*1000000 and int(cpgline.split("\t")[3]) < (x+1)*1000000 :
            cpg_split_lines.append(cpgline)
    del cpglines
    
    print "tile : reading CHH files", index
    
    chhfile = open("CHH_context_"+str(name)+"_chr"+str(chr)+"_sorted.out", 'r')
    chhlines = chhfile.readlines()
    chh_split_lines = []
    for chhline in chhlines :
        if int(chhline.split("\t")[3]) >= x*1000000 and int(chhline.split("\t")[3]) < (x+1)*1000000 :
            chh_split_lines.append(chhline)
    del chhlines
    
    print "tile : reading CHG files", index
    chgfile = open("CHG_context_"+str(name)+"_chr"+str(chr)+"_sorted.out", 'r')
    chglines = chgfile.readlines()
    chg_split_lines = []
    for chgline in chglines :
        if int(chgline.split("\t")[3]) >= x*1000000 and int(chgline.split("\t")[3]) < (x+1)*1000000 :
            chg_split_lines.append(chgline)
    del chglines
    print "finished reading", name, index, datetime.datetime.now()

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
    print "finished writing", datetime.datetime.now()

def samtobam(file):
    os.system("/opt/apps/samtools/0.1.18/samtools view -bS "+file+".sam > "+file+".bam")
def mergebam(name):
    file_list = []
    for file in os.listdir(name) :
        if str(file).split(".")[1] == "bam" :
            file_list.append(file)
    os.chdir(name)
    file_list_scr = " ".join(file_list)
    print "/opt/apps/samtools/0.1.18/samtools merge "+name+"_bt202.bam "+str(file_list_scr)
    os.system("/opt/apps/samtools/0.1.18/samtools merge "+name+"_bt202.bam "+str(file_list_scr))
def sortbam(name):
    os.chdir(name)
    print "/opt/apps/samtools/0.1.18/samtools sort "+name+"_bt202.bam "+name+"_bt202_sorted"
    os.system("/opt/apps/samtools/0.1.18/samtools sort "+name+"_bt202.bam "+name+"_bt202_sorted")
    
if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("--run", default=None, help="select from: split, bismark, combine_all, tile, bis2bed, prep_methylKit, methylKit")
    parser.add_option("--name", default=None, help="select from: B73_all3, Mo17_all3, Oh43_all3, CML322_all3, Tx303_all3 if multiple separate by ,")
    parser.add_option("--specie", default="maize", help="select from: arabidopsis, maize")
    parser.add_option("--met_type", default=None, help="select from: CpG, CHH, CHG, All")
    parser.add_option("--chrom", default=None, help="select from: 1-5, 1-10, 10")
    parser.add_option("--numcore", default=None, help="select from: 12, 24")
    parser.add_option("--index", default=None, help="if 0, go from 0 to 9, if 1, go from 10 to 19 etc..")
    options, args = parser.parse_args()
    
    if options.run == "bismark" :
        jobs = []
        for name in options.name.split(",") :
            o = multiprocessing.Process(target=bismark, args=(name, ))
            jobs.append(o)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
                
    elif options.run == "split" :
        if options.specie == "maize" :
            chrmax = 11
        jobs = []
        for name in options.name.split(",") :
            for met in options.met_type.split(",") :
                for chr in range(1,chrmax): 
                    p = multiprocessing.Process(target=split_by_chr, args=(chr, met, name,))
                    jobs.append(p)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
        
    elif options.run == "combine_all" :
        jobs = []
        if options.specie == "maize" :
            chrmax = 11
        for chr in range(1,chrmax):
            for name in options.name.split(",") :
                q = multiprocessing.Process(target=combine_type, args=(chr, name, ))
                jobs.append(q)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
    
    elif options.run == "bis2bed" :
        jobs = []
        if options.specie == "maize" :
            chrmax = 11
        for name in options.name.split(",") :
            for met in options.met_type.split(",") :
                for chr in range(1,chrmax) :
                    r = multiprocessing.Process(target=bis2bed, args=(chr, met, name, ))
                    jobs.append(r)
        [x.start() for x in jobs]
        [x.join()  for x in jobs]
        
    elif options.run == "prep_methylKit" :
        jobs = []
        if options.specie == "maize" :
            chrmax = 11
        for name in options.name.split(",") :
            for met in options.met_type.split(",") :
                for chr in range(1,chrmax) :
                    s = multiprocessing.Process(target=pre_methylKit, args=(chr, met, name, ))
                    jobs.append(s)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
        
    elif options.run == "methylKit" :
        jobs = []
        if options.specie == "maize" :
            chrmax = 11
        if len(options.name.split(",")) != 2 :
            print "enter two genotypes separated by comma on --name option"
            sys.exit()
        for met in options.met_type.split(",") :
            for chr in range(1,chrmax) :
                t = multiprocessing.Process(target=methylkit, args=(chr, met, options.name.split(",")[0],  options.name.split(",")[1], ))
                jobs.append(t)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
        
    elif options.run == "merge_strands" :
        jobs = []
        for name in names :
            u = multiprocessing.Process(target=merge_strands, args=("10", "CpG", name,))
            jobs.append(u)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
    
    elif options.run == "merge_chrom" :
        jobs = []
        for met in met_type :
            v = multiprocessing.Process(target=merge_chrom, args=(met, ))
            jobs.append(v)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
                    
    elif options.run == "wmet" :
        jobs = []
        for name in names :
            for met in met_type :
                for chr in range(1,11) :
                    w = multiprocessing.Process(target=wmet, args=(chr, met, name, ))
                    jobs.append(w)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
    elif options.run == "merge_wmet" :
        jobs = []
        for name in names :
            for met in met_type :
                z = multiprocessing.Process(target=merge_wmet, args=(met,name,))
                jobs.append(z)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
    elif options.run == "samtobam" :
        jobs = []
        for file in os.listdir(".") :
            if str(file).split('.')[1] == "sam" :
                filename = str(file).split('.')[0]
                w = multiprocessing.Process(target=samtobam, args=(filename,))
                jobs.append(w)
        [x.start() for x in jobs]
        [x.join() for x in jobs]        
    elif options.run == "mergebam" :
        jobs = []
        names = ["Mo17_all3", "Oh43_all3", "Tx303_all3", "CML322_all3"]
        for name in names :
            w = multiprocessing.Process(target=mergebam, args=(name,))
            jobs.append(w)
        [x.start() for x in jobs]
        [x.join() for x in jobs]   
    elif options.run == "sortbam" :
        jobs = []
        names = ["Mo17_all3", "Oh43_all3", "Tx303_all3", "CML322_all3"]
        for name in names :
            w = multiprocessing.Process(target=sortbam, args=(name,))
            jobs.append(w)
        [x.start() for x in jobs]
        [x.join() for x in jobs]            
    elif options.run == "tile" :
        jobs = []
        for name in options.name.split(",") : 
            for chr in options.chrom.split(",") :
                max_index = int(float(max_coord(name, chr))/1000000)+1
                large_list = []
                small_list = []
                for x in range(0, max_index+1) :
                    if len(small_list) < int(options.numcore) :
                        small_list.append(x)
                    elif len(small_list) == int(options.numcore) :
                        large_list.append(small_list)
                        small_list = []
                        small_list.append(x)
                if len(small_list) > 0 :
                    large_list.append(small_list)
                del small_list
                for small_list in large_list[0:1] :
                    jobs = []
                    print small_list, name, chr
                    sys.exit()
                    for index in small_list :
                        y = multiprocessing.Process(target=tile, args=(name, chr, index,))
                        jobs.append(y)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
    
    elif options.run == "tile_all" :
        jobs = []
        ### after methylation extraction split 
        if options.specie == "maize" :
            chrmax = 11
        jobs = []
        for name in options.name.split(",") :
            for met in options.met_type.split(",") :
                for chr in range(1,chrmax): 
                    p = multiprocessing.Process(target=split_by_chr, args=(chr, met, name,))
                    jobs.append(p)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
        ### tile
        jobs = []
        for name in options.name.split(",") :
            for chr in range(10,11) :
                for index in range(1,10) :
                    y = multiprocessing.Process(target=tile, args=(name, chr, index,))
                    jobs.append(y)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
        
        ### merge files and delete tmp files
        
    elif sys.argv[1] == "tile_argv" :
        tile(sys.argv[2], sys.argv[3], sys.argv[4])       