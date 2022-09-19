import os,sys
import threading

exitFlag = 0

def make_split_input(type, name, chr) :
    print "split chromosomes"
    chr_list = ['chr1','chr2','chr3','chr4','chr5','chr6', 'chr7', 'chr8','chr9','chr10']
    outfile = open(str(type)+"_context_"+str(name).split(".")[0]+"_chr"+str(chr)+".out", 'w')
    for a in open(str(type)+"_context_"+str(name).split(".")[0]+".sorted.txt", 'r') :
        if len(a.split("\t")) > 2 and str(a.split("\t")[2]) == "chr"+str(chr) :  
            outfile.write(a) 
    outfile.close()

def combine_context(name,chr) :
    print "combine contexts"
    if os.path.isfile("All_context_"+str(name).split(".")[0]+"_chr"+str(int(chr)+1)+".out") :
        print "file exists, exit"
        sys.exit()
    list = ['CpG', 'CHH', 'CHG']
    for type in list :
        os.system("cat "+str(type)+"_context_"+str(name).split(".")[0]+"_chr"+str(int(chr)+1)+".out >> All_context_"+str(name).split(".")[0]+"_chr"+str(int(chr)+1)+".out")

def bis2bed(type,name,chr) :
    os.system("/work/02114/wonaya/software/bismark2bedGraph --counts --CX_context -o "+str(type)+"_context_"+str(name).split(".")[0]+"_chr"+str(int(chr)+1)+".bedGraph "+str(type)+"_context_"+str(name).split(".")[0]+"_chr"+str(int(chr)+1)+".out")
    
def pre_methylKit(type, name, chr) :
    pos = []
    neg = []
    dict = {}
    ## save strand status
    for a in open(str(type)+"_context_"+str(name).split(".")[0]+"_chr"+str(int(chr)+1)+".out", 'r') : # example: CpG_context_B73_ACAGTG_merged.sorted_chr10.out bsseq/methylKit/swDMR
        dict[a.split("\t")[3]] = a.split("\t")[1]
        
    outfile = open(str(type)+"_context_"+str(name).split(".")[0]+"_chr"+str(int(chr)+1)+".methylKit", 'w')
    for b in open(str(type)+"_context_"+str(name).split(".")[0]+"_chr"+str(int(chr)+1)+".bedGraph", 'r') :
        if  dict[str(int(b.split("\t")[1])+1)] == "+" :
            strand = "F"
        elif  dict[str(int(b.split("\t")[1])+1)] == "-" :
            strand = "R"
        #print str(b.split("\t")[0])+"."+str(b.split("\t")[1]), str(b.split("\t")[0]), str(b.split("\t")[1]), strand
        total = float(int(b.split("\t")[4])+int(b.split("\t")[5].strip("\n")))
        freqc = float(int(b.split("\t")[4]))
        freqt = float(int(b.split("\t")[5].strip("\n")))
        percc = freqc/total*100
        perct = freqt/total*100
        #print str(b.split("\t")[0])+"."+str(b.split("\t")[1]), str(b.split("\t")[0]), str(b.split("\t")[1]), strand, int(total), percc, perct
        outfile.write(str(b.split("\t")[0])+"."+str(b.split("\t")[1]))
        outfile.write("\t")
        outfile.write(str(b.split("\t")[0])+"\t"+str(b.split("\t")[1]))
        outfile.write("\t")
        outfile.write(str(strand)+"\t"+str(int(total))+"\t"+str(percc)+"\t"+str(perct)+"\n")
    outfile.close()

def rtest(type, name, chr) :
    os.system("Rscript /work/02114/wonaya/scripts/r_methylkit.R "+str(type)+"_context_B73_ACAGTG_merged_chr"+str(chr)+".methylKit "+str(type)+"_context_Mo17_GCCAAT_merged_chr"+str(chr)+".methylKit "+str(type)+" "+str(type)+"_context_chr"+str(chr)+".txt")
    
def post_methylKit(type, name, chr) :
    ## save DMR
    dmr_loc = {}
    dmr_loc_start_val = []
    for c in open(str(type)+"_context_chr"+str(int(chr)+1)+".txt", 'r') :
        if c.split("\t")[1].strip('"') == "chr"+str(int(chr)+1) :
            dmr_loc[int(c.split("\t")[2])] = int(c.split("\t")[3])
            dmr_loc_start_val.append(int(c.split("\t")[2]))     
    dmr_loc_start_val.sort()   
    
    ## calculate weighted methylation
    a = open(str(type)+"_context_"+str(name).split(".")[0]+"_chr"+str(int(chr)+1)+".out", 'r') 
    alines = a.readlines()
    outfile = open(str(type)+"_context_"+strain+"_merged.sorted_chr"+str(int(chr)+1)+"_calc_wMet.out", 'w')
    outfile.write("chrom\tstart\total C\ttotal Read\ttotal valid C\ttotal valid reads\tFrag\tMean\tWeighted\n") 
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
        outfile.write("chr"+str(int(chr)+1)+"\t")
        outfile.write(str(dmr)+"\t")
        outfile.write(str(dmr_loc[dmr])+"\t")
        outfile.write(str(total_c)+"\t")
        outfile.write(str(total_read)+"\t")
        outfile.write(str(val_met)+"\t")
        outfile.write(str(weighted_top)+"\t")
        outfile.write(str(frag)+"\t"+str(mean)+"\t"+str(weighted)+"\n")
    outfile.close()
    
def combine_chromosome(type) :
    outfile = open(str(type)+"_context_allchr.out", 'w')
    for x in range(1,11) :
        alines = []
        a = open(str(type)+"_context_chr"+str(x)+".out", 'r')
        alines = a.readlines()  
        if x == 1 :
            for aline in alines :
                outfile.write(aline)
        else :
            for aline in alines[1:] :
                outfile.write(aline)
    outfile.close()    

def main() :
    if sys.argv[1] == "split" : 
        make_split_input(sys.argv[2],sys.argv[3],sys.argv[4])
    elif sys.argv[1] == "all" :    
        combine_context(sys.argv[3],sys.argv[4])
    elif sys.argv[1] == "prep" :
        bis2bed(sys.argv[2],sys.argv[3],sys.argv[4])
        pre_methylKit(sys.argv[2],sys.argv[3],sys.argv[4])
        rtest(sys.argv[2],sys.argv[3],sys.argv[4])
    
    #post_methylKit(sys.argv[1],sys.argv[2],sys.argv[3])
    elif sys.argv[1] == "comb" :
        combine_chromosome(sys.argv[1])

if __name__ == "__main__":
    main()