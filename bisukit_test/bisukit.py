#!/usr/bin/env python
import os,sys
import resource
from optparse import OptionParser,OptionGroup
from datetime import datetime
try:
   import subprocess
except ImportError:
   print >> sys.stderr,"Could not import the subprocess module"
   raise
import multiprocessing
import readline, glob
import subprocess

def complete(text, state):
    return (glob.glob(text+'*')+[None])[state]
readline.set_completer_delims(' \t\n;')
readline.parse_and_bind("tab: complete")
readline.set_completer(complete)

class setup :
    @staticmethod
    def get_genome_size(genomefile) :
        chr_list = []
        count_list = []
        count = 0
        for genome in open(genomefile, 'r') :
            if genome[0] == ">" :
                if count > 0 :
                    count_list.append(count)
                chr_list.append(genome[1:].split(" ")[0].strip("\n"))
                count = 0
            else :
                count += len(genome)
        count_list.append(count)
        chr_dict = {}
        for chr in chr_list :
            chr_dict[chr] = count_list[chr_list.index(chr)]
        return chr_dict, chr_list
              
class bismark :
    @staticmethod
    def methylation_extractor(bismark, samfile, paired) :
        if not os.path.isfile(str(bismark)+"_methylation_extractor") :
            print "check the path to bismark_methylation_extractor"
        else :
            if paired == "yes" :
                os.system(str(bismark)+"_methylation_extractor -p --report --bedGraph "+str(samfile)) 
            elif paired == "no" :
                os.system(str(bismark)+"_methylation_extractor -s --report --bedGraph "+str(samfile)) 

class methylkit:
    @staticmethod
    def prefilter(chr, samname, context) :
        print "filtering by chr", chr, context, datetime.now()
        outfile = open(str(context)+"_OT_"+str(samname).split(".sam")[0]+"_chr"+str(chr)+".txt", 'w')
        for line1 in open(str(context)+"_OT_"+str(samname).split(".sam")[0]+".txt", 'r') :
            if len(line1.split("\t")) > 1 and line1.split("\t")[2] == "chr"+str(chr) :
                outfile.write(line1)
        outfile.close()
        outfile = open(str(context)+"_CTOT_"+str(samname).split(".sam")[0]+"_chr"+str(chr)+".txt", 'w')
        for line1 in open(str(context)+"_CTOT_"+str(samname).split(".sam")[0]+".txt", 'r') :
            if len(line1.split("\t")) > 1 and line1.split("\t")[2] == "chr"+str(chr) :
                outfile.write(line1)
        outfile.close()
        outfile = open(str(context)+"_OB_"+str(samname).split(".sam")[0]+"_chr"+str(chr)+".txt", 'w')
        for line1 in open(str(context)+"_OB_"+str(samname).split(".sam")[0]+".txt", 'r') :
            if len(line1.split("\t")) > 1 and line1.split("\t")[2] == "chr"+str(chr) :
                outfile.write(line1)
        outfile.close()
        outfile = open(str(context)+"_CTOB_"+str(samname).split(".sam")[0]+"_chr"+str(chr)+".txt", 'w')
        for line1 in open(str(context)+"_CTOB_"+str(samname).split(".sam")[0]+".txt", 'r') :
            if len(line1.split("\t")) > 1 and line1.split("\t")[2] == "chr"+str(chr) :
                outfile.write(line1)
        outfile.close()
        
    @staticmethod
    def prep_methylkit_nonb73(chr, write_chr, samname_ori, context, largemem) :
        import bisect
        if len(samname_ori.split("/")) > 1 :
            os.chdir(samname_ori.split("/")[0])
            samname = samname_ori.split("/")[1]
        else :
            samname = samname_ori
        genome_list = []
        whole_list = []
        chr_list = []
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
        if not os.path.isfile(str(context)+"_"+str(samname)+"_chr"+str(chr)+".methylKit2") :
            print "chr", chr, samname, context, "writing", str(context)+"_"+str(samname)+"_chr"+str(chr)+".methylKit"
            outfile = open(str(context)+"_"+str(samname)+"_chr"+str(chr)+".methylKit" ,'w')
            ## save strand status
            if largemem == "yes" :
                max_cycle = 4
                max_list = 100000000
            else :
                max_cycle = 40
                max_list = 10000000
            cg_loc_list = []
            if context == "CpG" :
                cg_loc_list = []
                loc = 1
                for base in whole_list[int(chr)-1] :
                    if whole_list[int(chr)-1][loc-1:loc+1] == "CG" :
                        cg_loc_list.append(loc)
                        cg_loc_list.append(loc+1)
                    loc += 1
            elif context == "CHG" :
                chg_list = ["CAG", "CTG", "CCG"]
                chg_rev = ["CTG", "CAG", "CGG"]
                cg_loc_list = []
                loc = 1 
                for base in whole_list[int(chr)-1] : 
                    if whole_list[int(chr)-1][loc-1:loc+2] in chg_list :
                        cg_loc_list.append(loc)
                    elif whole_list[int(chr)-1][loc-3:loc] in chg_rev :
                        cg_loc_list.append(loc)
                    loc += 1
            elif context == "CHH" :
                chh_list = ["CAA", "CAC", "CAT", "CCA", "CCC", "CCT", "CTA", "CTC", "CTT"]
                chh_rev = ["TTG", "GTG", "ATG", "TGG", "GGG", "AGG", "TAG", "GAG", "AAG"]
                cg_loc_list = []
                loc = 1 
                for base in whole_list[int(chr)-1] : 
                    if whole_list[int(chr)-1][loc-1:loc+2] in chh_list :
                        cg_loc_list.append(loc)
                    elif whole_list[int(chr)-1][loc-3:loc] in chh_rev :
                        cg_loc_list.append(loc)
                    loc += 1
            print len(cg_loc_list)
            for x in range(0, max_cycle) : ## max set to 4 0 to 399Mbp
                print "cycle", x, ":", "chr", chr, samname, context, x*max_list, (x+1)*max_list ,datetime.now()
                redo_cg_loc_list = []
                for cg_loc in cg_loc_list :
                    if cg_loc >= x*max_list and cg_loc < (x+1)*max_list : 
                        redo_cg_loc_list.append(cg_loc)
                if len(redo_cg_loc_list) == 0 :
                    break
                large_list = []
                for y in range(0,max_list+1): ### 56 percent of lonestar node, 2min to make
                    large_list.append([0]*4)
                print chr, context, samname, "reading OT",datetime.now()
                for a_f in open(str(context)+"_OT_"+str(samname).split(".sam")[0]+"_chr"+str(chr)+".txt", 'r') : # example: CpG_OT_11510-pool_GGTAGC_L004_R1_001_val_1.fq_bismark_bt2_pe.txt
                    if len(a_f.split("\t")) > 2 and len(a_f.split("\t")) != 0 and int(a_f.split("\t")[3]) >= x*max_list and int(a_f.split("\t")[3]) < (x+1)*max_list :
                        if a_f.split("\t")[1] == "+" :
                            large_list[int(float(a_f.split("\t")[3]))-(x*max_list)][0] += 1
                        elif a_f.split("\t")[1] == "-" :
                            large_list[int(float(a_f.split("\t")[3]))-(x*max_list)][1] += 1
                if os.path.isfile(str(context)+"_CTOT_"+str(samname).split(".sam")[0]+"_chr"+str(chr)+".txt"):
                    for a_f in open(str(context)+"_CTOT_"+str(samname).split(".sam")[0]+"_chr"+str(chr)+".txt", 'r') : # example: 
                        if len(a_f.split("\t")) > 2 and len(a_f.split("\t")) != 0 and int(a_f.split("\t")[3]) >= x*max_list and int(a_f.split("\t")[3]) < (x+1)*max_list :
                            if a_f.split("\t")[1] == "+" :
                                large_list[int(float(a_f.split("\t")[3]))-(x*max_list)][0] += 1
                            elif a_f.split("\t")[1] == "-" :
                                large_list[int(float(a_f.split("\t")[3]))-(x*max_list)][1] += 1
                print chr, context, samname, "reading OB",datetime.now()
                for a_r in open(str(context)+"_OB_"+str(samname).split(".sam")[0]+"_chr"+str(chr)+".txt", 'r') : 
                    if len(a_r.split("\t")) > 2 and len(a_r.split("\t")) != 0 and int(a_r.split("\t")[3]) >= x*max_list and int(a_r.split("\t")[3]) < (x+1)*max_list :
                        if a_r.split("\t")[1] == "+" :
                            large_list[int(float(a_r.split("\t")[3]))-(x*max_list)][2] += 1
                        elif a_r.split("\t")[1] == "-" :
                            large_list[int(float(a_r.split("\t")[3]))-(x*max_list)][3] += 1
                if os.path.isfile(str(context)+"_CTOB_"+str(samname).split(".sam")[0]+"_chr"+str(chr)+".txt"):
                    for a_r in open(str(context)+"_CTOB_"+str(samname).split(".sam")[0]+"_chr"+str(chr)+".txt", 'r') : 
                        if len(a_r.split("\t")) > 2 and len(a_r.split("\t")) != 0 and int(a_r.split("\t")[3]) >= x*max_list and int(a_r.split("\t")[3]) < (x+1)*max_list :
                            if a_r.split("\t")[1] == "+" :
                                large_list[int(float(a_r.split("\t")[3]))-(x*max_list)][2] += 1
                            elif a_r.split("\t")[1] == "-" :
                                large_list[int(float(a_r.split("\t")[3]))-(x*max_list)][3] += 1
                for cg_loc in redo_cg_loc_list :
                    if int(large_list[cg_loc-(x*max_list)][0])+int(large_list[cg_loc-(x*max_list)][1]) != 0 :
                        outfile.write(str(cg_loc)+"\t")
                        outfile.write(str(write_chr)+"\t")
                        outfile.write(str(cg_loc)+"\t")
                        outfile.write("1\t") ## change when later Forward and reverse strands are read in separately
                        outfile.write(str(int(large_list[cg_loc-(x*max_list)][0])+int(large_list[cg_loc-(x*max_list)][1]))+"\t")
                        outfile.write(str(float(large_list[cg_loc-(x*max_list)][0])/float(int(large_list[cg_loc-(x*max_list)][0])+int(large_list[cg_loc-(x*max_list)][1]))*100)+"\t")
                        outfile.write(str(float(large_list[cg_loc-(x*max_list)][1])/float(int(large_list[cg_loc-(x*max_list)][0])+int(large_list[cg_loc-(x*max_list)][1]))*100)+"\n")
                    elif int(large_list[cg_loc-(x*max_list)][2])+int(large_list[cg_loc-(x*max_list)][3]) != 0 :
                        outfile.write(str(cg_loc)+"\t")
                        outfile.write(str(write_chr)+"\t")
                        outfile.write(str(cg_loc)+"\t")
                        outfile.write("0\t") ## change when later Forward and reverse strands are read in separately
                        outfile.write(str(int(large_list[cg_loc-(x*max_list)][2])+int(large_list[cg_loc-(x*max_list)][3]))+"\t")
                        outfile.write(str(float(large_list[cg_loc-(x*max_list)][2])/float(int(large_list[cg_loc-(x*max_list)][2])+int(large_list[cg_loc-(x*max_list)][3]))*100)+"\t")
                        outfile.write(str(float(large_list[cg_loc-(x*max_list)][3])/float(int(large_list[cg_loc-(x*max_list)][2])+int(large_list[cg_loc-(x*max_list)][3]))*100)+"\n")
                del large_list, redo_cg_loc_list
            outfile.close()
        if len(samname_ori.split("/")) > 1 :
            os.chdir("..")
            
    @staticmethod
    def prep_methylkit(chr, write_chr, samname_ori, context, largemem) :
        if len(samname_ori.split("/")) > 1 :
            os.chdir(samname_ori.split("/")[0])
            samname = samname_ori.split("/")[1]
        else :
            samname = samname_ori
        if not os.path.isfile(str(context)+"_"+str(samname)+"_chr"+str(chr)+".methylKit") :
            print "writing", str(context)+"_"+str(samname)+"_chr"+str(chr)+".methylKit"
            outfile = open(str(context)+"_"+str(samname)+"_chr"+str(chr)+".methylKit" ,'w')
            ## save strand status
            if largemem == "yes" :
                max_cycle = 4
                max_list = 100000000
            else :
                max_cycle = 40
                max_list = 10000000
            for x in range(0, max_cycle) : ## max set to 4 0 to 399Mbp
                large_list = []
                for y in range(0,max_list): ### 56 percent of lonestar node, 2min to make
                    large_list.append([0]*4)
                print samname, "OT", datetime.now()
                for a_f in open(str(context)+"_OT_"+str(samname).split(".sam")[0]+".txt", 'r') : # example: CpG_OT_11510-pool_GGTAGC_L004_R1_001_val_1.fq_bismark_bt2_pe.txt
                    if len(a_f.split("\t")) > 2 and a_f.split("\t")[2] == "chr"+str(chr) and len(a_f.split("\t")) != 0 and int(a_f.split("\t")[3]) >= x*max_list and int(a_f.split("\t")[3]) < (x+1)*max_list :
                        if a_f.split("\t")[1] == "+" :
                            large_list[int(float(a_f.split("\t")[3]))-(x*max_list)][0] += 1
                        elif a_f.split("\t")[1] == "-" :
                            large_list[int(float(a_f.split("\t")[3]))-(x*max_list)][1] += 1
                print samname, "CTOT", datetime.now()
                if os.path.isfile(str(context)+"_CTOT_"+str(samname).split(".sam")[0]+".txt") :
                    for a_f in open(str(context)+"_CTOT_"+str(samname).split(".sam")[0]+".txt", 'r') : # example: 
                        if len(a_f.split("\t")) > 2 and a_f.split("\t")[2] == "chr"+str(chr) and len(a_f.split("\t")) != 0 and int(a_f.split("\t")[3]) >= x*max_list and int(a_f.split("\t")[3]) < (x+1)*max_list :
                            if a_f.split("\t")[1] == "+" :
                                large_list[int(float(a_f.split("\t")[3]))-(x*max_list)][0] += 1
                            elif a_f.split("\t")[1] == "-" :
                                large_list[int(float(a_f.split("\t")[3]))-(x*max_list)][1] += 1
                
                print samname, "OB", datetime.now()
                for a_r in open(str(context)+"_OB_"+str(samname).split(".sam")[0]+".txt", 'r') : 
                    if len(a_r.split("\t")) > 2 and a_r.split("\t")[2] == "chr"+str(chr) and len(a_r.split("\t")) != 0 and int(a_r.split("\t")[3]) >= x*max_list and int(a_r.split("\t")[3]) < (x+1)*max_list :
                        if a_r.split("\t")[1] == "+" :
                            large_list[int(float(a_r.split("\t")[3]))-(x*max_list)][2] += 1
                        elif a_r.split("\t")[1] == "-" :
                            large_list[int(float(a_r.split("\t")[3]))-(x*max_list)][3] += 1
                print samname, "CTOB", datetime.now()
                if os.path.isfile(str(context)+"_CTOB_"+str(samname).split(".sam")[0]+".txt") :
                    for a_r in open(str(context)+"_CTOB_"+str(samname).split(".sam")[0]+".txt", 'r') : 
                        if len(a_r.split("\t")) > 2 and a_r.split("\t")[2] == "chr"+str(chr) and len(a_r.split("\t")) != 0 and int(a_r.split("\t")[3]) >= x*max_list and int(a_r.split("\t")[3]) < (x+1)*max_list :
                            if a_r.split("\t")[1] == "+" :
                                large_list[int(float(a_r.split("\t")[3]))-(x*max_list)][2] += 1
                            elif a_r.split("\t")[1] == "-" :
                                large_list[int(float(a_r.split("\t")[3]))-(x*max_list)][3] += 1
                z = 0
                for list in large_list :
                    if int(list[0])+int(list[1]) != 0 :
                        outfile.write(str(int(z)+(int(x)*max_list))+"\t")
                        outfile.write(str(write_chr)+"\t")
                        outfile.write(str(int(z)+(int(x)*max_list))+"\t")
                        outfile.write("1\t") ## change when later Forward and reverse strands are read in separately
                        outfile.write(str(int(list[0])+int(list[1]))+"\t")
                        outfile.write(str(float(list[0])/float(int(list[0])+int(list[1]))*100)+"\t")
                        outfile.write(str(float(list[1])/float(int(list[0])+int(list[1]))*100)+"\n")
                    elif int(list[2])+int(list[3]) != 0 :
                        outfile.write(str(int(z)+(int(x)*max_list))+"\t")
                        outfile.write(str(write_chr)+"\t")
                        outfile.write(str(int(z)+(int(x)*max_list))+"\t")
                        outfile.write("0\t") ## change when later Forward and reverse strands are read in separately
                        outfile.write(str(int(list[2])+int(list[3]))+"\t")
                        outfile.write(str(float(list[2])/float(int(list[2])+int(list[3]))*100)+"\t")
                        outfile.write(str(float(list[3])/float(int(list[2])+int(list[3]))*100)+"\n")
                    z += 1
                del large_list
            outfile.close()
        if len(samname_ori.split("/")) > 1 :
            os.chdir("..")
        
    @staticmethod
    def methylkit_edmr(file1, file2, type, chr, specie) :
        try :
            from rpy2.robjects.packages import importr
        except :
            print "set R in $PATH" ; sys.exit()
        import rpy2.robjects as robjects
        
        print "generating eDMR for", file1, file2, "on", chr, type
        importr('GenomicRanges')
        importr('bigmemory')
        importr('data.table')
        robjects.r('''source('/work/02114/wonaya/scripts/base.r')''')
        filelist = [file1, file2]
        if len(file1.split("/")) > 1 : 
            file1name = file1.split("/")[0]+"/"+str(type)+"_"+file1.split("/")[1]+"_chr"+str(chr)+".methylKit"
            file1ID = file1.split("/")[1].split("_")[0]
        else :
            file1name = str(type)+"_"+file1+"_"+str(chr)+".methylKit"
            file1ID = file1.split("_")[0]
        if len(file2.split("/")) > 1 : 
            file2name = file2.split("/")[0]+"/"+str(type)+"_"+file2.split("/")[1]+"_chr"+str(chr)+".methylKit"
            file2ID = file2.split("/")[1].split("_")[0]
        else :
            file2name = str(type)+"_"+file2+"_chr"+str(chr)+".methylKit"
            file2ID = file2.split("_")[0]
        filelist = [file1name, file2name]
        robjects.r.assign('rfilelist',filelist)
        robjects.r.assign('type', type)
        robjects.r.assign('file1ID', file1ID)
        robjects.r.assign('file2ID', file2ID)
        robjects.r.assign('specie', specie)
        outfile = str(file1ID)+"_"+str(file2ID)+"_"+str(type)+"_dmr_chr"+str(chr)+".txt"
        robjects.r.assign('outfile', outfile)
        robjects.r('''suppressMessages(library(base))''')
        robjects.r('''myobj=read(rfilelist, sample.id=list(file1ID,file2ID), assembly=specie,treatment=c(0,1),context=type, resolution="base")''')
        robjects.r('''meth=unite(myobj, destrand=FALSE)''')
        robjects.r('''rm(valid.methylRawObj)''')
        robjects.r('''rm(rfilelist)''')
        robjects.r('''rm(read)''')
        robjects.r('''rm(myobj)''')
        robjects.r('''gc()''')
        robjects.r('''source("/work/02114/wonaya/scripts/diffMeth.R")''')
        robjects.r('''myDiff=calculateDiffMeth(meth)''')
        robjects.r('''rm(meth)''')
        robjects.r('''gc()''')
        robjects.r('''source("/work/02114/wonaya/scripts/eDMR_test.R")''')
        robjects.r('''myDMR=eDMR(myDiff,step=100,dist="none",DMC.qvalue=1,DMC.methdiff=25,num.DMCs=1,num.CpGs=3,DMR.methdiff=20,plot=FALSE,main="example",mode=1)''')
        robjects.r('''myDMR.sig=filter.dmr(myDMR, mean.meth.diff=20, num.CpGs=3, num.DMCs=1)''')
        robjects.r('''myDMR.sig.df <- as.data.frame(myDMR.sig)''')
        robjects.r('''write.table(myDMR.sig.df, file=outfile, sep="\t")''')
        
    @staticmethod
    def rewrite_dmr(file1, file2, type, chr, specie) :
        filelist = [file1, file2]
        if len(file1.split("/")) > 1 : 
            file1name = file1.split("/")[0]+"/"+str(type)+"_"+file1.split("/")[1]+"_chr"+str(chr)+".methylKit"
            file1ID = file1.split("/")[1].split("_")[0]
        else :
            file1name = str(type)+"_"+file1+"_"+str(chr)+".methylKit"
            file1ID = file1.split("_")[0]
        if len(file2.split("/")) > 1 : 
            file2name = file2.split("/")[0]+"/"+str(type)+"_"+file2.split("/")[1]+"_chr"+str(chr)+".methylKit"
            file2ID = file2.split("/")[1].split("_")[0]
        else :
            file2name = str(type)+"_"+file2+"_"+str(chr)+".methylKit"
            file2ID = file2.split("_")[0]
        infile = open(str(file1ID)+"_"+str(file2ID)+"_"+str(type)+"_dmr_chr"+str(chr)+".txt", 'r')
        inlines = infile.readlines()
        infile.close()
        outfile = open(str(file1ID)+"_"+str(file2ID)+"_"+str(type)+"_dmr_chr"+str(chr)+"_temp.txt", 'w')
        outfile.write("chr\tstart\tend\twidth\tstrand\tmean.meth.diff\tnum.CpG\tnum.DMC\tDMR.pvalue\tDMR.qvalue\n")
        for line in inlines :
            if line.split("\t")[0] != '"seqnames"' :
                for a in line.split("\t")[1:] :
                    if line.split("\t").index(a) == 1 :
                        outfile.write(str(chr)+"\t")
                    else :
                        outfile.write(a.strip('"\n')+"\t")
                outfile.write("\n")
        outfile.close()
        os.system("mv "+str(file1ID)+"_"+str(file2ID)+"_"+str(type)+"_dmr_chr"+str(chr)+"_temp.txt "+str(file1ID)+"_"+str(file2ID)+"_"+str(type)+"_dmr_chr"+str(chr)+".txt")
    
    @staticmethod
    def merge_dmr(file1, file2, type, specie) :
        filelist = [file1, file2]
        if len(file1.split("/")) > 1 : 
            file1ID = file1.split("/")[1].split("_")[0]
        else :
            file1ID = file1.split("_")[0]
        if len(file2.split("/")) > 1 : 
            file2ID = file2.split("/")[1].split("_")[0]
        else :
            file2ID = file2.split("_")[0]
        
        chr_list = {}
        for chr in setup.get_genome_size(options.genome)[0].keys() :
            if chr[:3] == "chr" : 
                chr_list[str(chr).split("chr")[1]] = []
            else :
                chr_list[str(chr)] = []
        chr_list_int = []
        chr_list_str = []
        for chr in chr_list :
            if chr.isdigit() == True :
                chr_list_int.append(int(chr))
            else :
                chr_list_str.append(str(chr))
        
        chr_list = []
        chr_list_int.sort()
        chr_list_str.sort()
        chr_list =chr_list_int+chr_list_str
        
        outfile = open(str(file1ID)+"_"+str(file2ID)+"_"+str(type)+"_dmr.txt", 'w')
        for chr in chr_list :
            if os.path.isfile(str(file1ID)+"_"+str(file2ID)+"_"+str(type)+"_dmr_chr"+str(chr)+".txt") :
                a = open(str(file1ID)+"_"+str(file2ID)+"_"+str(type)+"_dmr_chr"+str(chr)+".txt", 'r') 
                alines = a.readlines()
                if chr_list.index(chr) == 0 :
                    for aline in alines :
                        outfile.write(aline)
                else :
                    for aline in alines[1:] :
                        outfile.write(aline)
        outfile.close()
        
class tile_class :
    @staticmethod
    def extract_name(samname, genomefile) :
        dir = "."
        if len(samname.split("/")) > 1 :
            dir = samname.split("/")[0]
            filename = samname.split("/")[1]
        else :
            dir = "."
            samname = samname
        name = samname.split(".sam")[0]
        return name, dir
    
    @staticmethod
    def tile_context_sp_single(chr, samname, genomefile,context,window) :
        large_list = []
        for x in range(0,setup.get_genome_size(genomefile)[0][chr]/int(window)+1) :
            large_list.append([0]*4)
        if os.path.isfile(str(context)+"_OB_"+tile_class.extract_name(samname, genomefile)[0]+".txt") and os.path.isfile(str(context)+"_OT_"+tile_class.extract_name(samname, genomefile)[0]+".txt") :
            print "reading OB"
            for line1 in open(str(context)+"_OT_"+tile_class.extract_name(samname, genomefile)[0]+".txt", 'r') :
                if line1[0:7] != "Bismark" and line1.split("\t")[2] == "chr"+str(chr)  :
                    if line1.split("\t")[1] == "+" :
                        large_list[int(float(line1.split("\t")[3])/int(window))][0] += 1 
                        large_list[int(float(line1.split("\t")[3])/int(window))][1] += 1
                    elif line1.split("\t")[1] == "-" :
                        large_list[int(float(line1.split("\t")[3])/int(window))][1] += 1
            for line1 in open(str(context)+"_CTOT_"+tile_class.extract_name(samname, genomefile)[0]+".txt", 'r') :
                if line1[0:7] != "Bismark" and line1.split("\t")[2] == "chr"+str(chr)  :
                    if line1.split("\t")[1] == "+" :
                        large_list[int(float(line1.split("\t")[3])/int(window))][0] += 1 
                        large_list[int(float(line1.split("\t")[3])/int(window))][1] += 1
                    elif line1.split("\t")[1] == "-" :
                        large_list[int(float(line1.split("\t")[3])/int(window))][1] += 1
            print "reading OT"
            for line2 in open(str(context)+"_OB_"+tile_class.extract_name(samname, genomefile)[0]+".txt", 'r') :
                if line2[0:7] != "Bismark" and line2.split("\t")[2] == "chr"+str(chr)  :
                    if line2.split("\t")[1] == "+" :
                        large_list[int(float(line2.split("\t")[3])/int(window))][2] += 1 
                        large_list[int(float(line2.split("\t")[3])/int(window))][3] += 1
                    elif line2.split("\t")[1] == "-" :
                        large_list[int(float(line2.split("\t")[3])/int(window))][3] += 1
            for line2 in open(str(context)+"_CTOB_"+tile_class.extract_name(samname, genomefile)[0]+".txt", 'r') :
                if line2[0:7] != "Bismark" and line2.split("\t")[2] == "chr"+str(chr)  :
                    if line2.split("\t")[1] == "+" :
                        large_list[int(float(line2.split("\t")[3])/int(window))][2] += 1 
                        large_list[int(float(line2.split("\t")[3])/int(window))][3] += 1
                    elif line2.split("\t")[1] == "-" :
                        large_list[int(float(line2.split("\t")[3])/int(window))][3] += 1
        else : 
            print str(context)+"_OB_"+tile_class.extract_name(samname, genomefile)[0]+".txt missing, run bismark_methylation_extractor" ; sys.exit()
        return large_list
    
    @staticmethod
    def tile_context_sp_single_nonb73(chr, samname, genomefile,context,window) :
        genome_list = []
        whole_list = []
        chr_list = []
        print samname, context, chr, "reading in B73 genome"
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
        print samname, context, chr, "making cg_loc_list",datetime.now()
        cg_loc_list = []
        if context == "CpG" :
            cg_loc_list = []
            loc = 1
            for base in whole_list[int(chr)-1] :
                if whole_list[int(chr)-1][loc-1:loc+1] == "CG" :
                    cg_loc_list.append(loc)
                    cg_loc_list.append(loc+1)
                loc += 1
        elif context == "CHG" :
            chg_list = ["CAG", "CTG", "CCG"]
            chg_rev = ["CTG", "CAG", "CGG"]
            cg_loc_list = []
            loc = 1 
            for base in whole_list[int(chr)-1] : 
                if whole_list[int(chr)-1][loc-1:loc+2] in chg_list :
                    cg_loc_list.append(loc)
                elif whole_list[int(chr)-1][loc-3:loc] in chg_rev :
                    cg_loc_list.append(loc)
                loc += 1
        elif context == "CHH" :
            chh_list = ["CAA", "CAC", "CAT", "CCA", "CCC", "CCT", "CTA", "CTC", "CTT"]
            chh_rev = ["TTG", "GTG", "ATG", "TGG", "GGG", "AGG", "TAG", "GAG", "AAG"]
            cg_loc_list = []
            loc = 1 
            for base in whole_list[int(chr)-1] : 
                if whole_list[int(chr)-1][loc-1:loc+2] in chh_list :
                    cg_loc_list.append(loc)
                elif whole_list[int(chr)-1][loc-3:loc] in chh_rev :
                    cg_loc_list.append(loc)
                loc += 1
        
        print "filtering ", chr, "input by", len(cg_loc_list), context, datetime.now()
        if not os.path.isfile(str(context)+"_OT_"+tile_class.extract_name(samname, genomefile)[0]+"_chr"+str(chr)+".txt") :
            outfile = open(str(context)+"_OT_"+tile_class.extract_name(samname, genomefile)[0]+"_chr"+str(chr)+".txt", 'w')
            for line1 in open(str(context)+"_OT_"+tile_class.extract_name(samname, genomefile)[0]+".txt", 'r') :
                if len(line1.split("\t")) > 1 and line1.split("\t")[2] == "chr"+str(chr) :
                    outfile.write(line1)
            outfile.close()
        if not os.path.isfile(str(context)+"_CTOT_"+tile_class.extract_name(samname, genomefile)[0]+"_chr"+str(chr)+".txt") :
            outfile = open(str(context)+"_CTOT_"+tile_class.extract_name(samname, genomefile)[0]+"_chr"+str(chr)+".txt", 'w')
            for line1 in open(str(context)+"_CTOT_"+tile_class.extract_name(samname, genomefile)[0]+".txt", 'r') :
                if len(line1.split("\t")) > 1 and line1.split("\t")[2] == "chr"+str(chr) :
                    outfile.write(line1)
            outfile.close()
        if not os.path.isfile(str(context)+"_OB_"+tile_class.extract_name(samname, genomefile)[0]+"_chr"+str(chr)+".txt") :
            outfile = open(str(context)+"_OB_"+tile_class.extract_name(samname, genomefile)[0]+"_chr"+str(chr)+".txt", 'w')
            for line1 in open(str(context)+"_OB_"+tile_class.extract_name(samname, genomefile)[0]+".txt", 'r') :
                if len(line1.split("\t")) > 1 and line1.split("\t")[2] == "chr"+str(chr) :
                    outfile.write(line1)
            outfile.close()
        if not os.path.isfile(str(context)+"_CTOB_"+tile_class.extract_name(samname, genomefile)[0]+"_chr"+str(chr)+".txt") :
            outfile = open(str(context)+"_CTOB_"+tile_class.extract_name(samname, genomefile)[0]+"_chr"+str(chr)+".txt", 'w')
            for line1 in open(str(context)+"_CTOB_"+tile_class.extract_name(samname, genomefile)[0]+".txt", 'r') :
                if len(line1.split("\t")) > 1 and line1.split("\t")[2] == "chr"+str(chr) :
                    outfile.write(line1)
            outfile.close()
        
        large_list_1 = []
        print "making large_list", datetime.now()
        for x in range(0,len(whole_list[int(chr)-1])+1) :
            large_list_1.append([0]*4)
        print "adding to large_list", datetime.now()
        for line1 in open(str(context)+"_OT_"+tile_class.extract_name(samname, genomefile)[0]+"_chr"+str(chr)+".txt", 'r') :
            if int(line1.split("\t")[3]) <= len(whole_list[int(chr)-1]) :
                if line1.split("\t")[1] == "+" :
                    large_list_1[int(line1.split("\t")[3])][0] += 1
                    large_list_1[int(line1.split("\t")[3])][1] += 1
                elif line1.split("\t")[1] == "-" :
                    large_list_1[int(line1.split("\t")[3])][1] += 1
        for line1 in open(str(context)+"_CTOT_"+tile_class.extract_name(samname, genomefile)[0]+"_chr"+str(chr)+".txt", 'r') :
            if int(line1.split("\t")[3]) <= len(whole_list[int(chr)-1]) :
                if line1.split("\t")[1] == "+" :
                    large_list_1[int(line1.split("\t")[3])][0] += 1
                    large_list_1[int(line1.split("\t")[3])][1] += 1
                elif line1.split("\t")[1] == "-" :
                    large_list_1[int(line1.split("\t")[3])][1] += 1
        for line1 in open(str(context)+"_OB_"+tile_class.extract_name(samname, genomefile)[0]+"_chr"+str(chr)+".txt", 'r') :
            if int(line1.split("\t")[3]) <= len(whole_list[int(chr)-1]) :
                if line1.split("\t")[1] == "+" :
                    large_list_1[int(line1.split("\t")[3])][2] += 1
                    large_list_1[int(line1.split("\t")[3])][3] += 1
                elif line1.split("\t")[1] == "-" :
                    large_list_1[int(line1.split("\t")[3])][3] += 1
        for line1 in open(str(context)+"_CTOB_"+tile_class.extract_name(samname, genomefile)[0]+"_chr"+str(chr)+".txt", 'r') :
            if int(line1.split("\t")[3]) <= len(whole_list[int(chr)-1]) :
                if line1.split("\t")[1] == "+" :
                    large_list_1[int(line1.split("\t")[3])][2] += 1
                    large_list_1[int(line1.split("\t")[3])][3] += 1
                elif line1.split("\t")[1] == "-" :
                    large_list_1[int(line1.split("\t")[3])][3] += 1
        print "done making, make large_list", datetime.now()
        large_list = []
        for x in range(0,len(whole_list[int(chr)-1])/int(window)+1) :
            large_list.append([0]*4)
        print "done makinglarge_list", datetime.now()
        for cg_loc in cg_loc_list :
            large_list[int(float(cg_loc)/float(window))][0] += large_list_1[cg_loc][0]
            large_list[int(float(cg_loc)/float(window))][1] += large_list_1[cg_loc][1]
            large_list[int(float(cg_loc)/float(window))][2] += large_list_1[cg_loc][2]
            large_list[int(float(cg_loc)/float(window))][3] += large_list_1[cg_loc][3]
        return large_list
    
    @staticmethod
    def tile_context_sp_single_write(chr, samname, genomefile,context,window,nonb73) :
        if len(samname.split("/")) > 1 :
            os.chdir(samname.split("/")[0])
            samname = samname.split("/")[1]
        outfile = open(str(samname).strip(".sam")+"_tile_"+str(context)+"_"+str(window)+"bp_"+str(chr)+".txt", 'w')
        if nonb73 == "yes" :
            large_list = tile_class.tile_context_sp_single_nonb73(chr, samname, genomefile,context,window)  
        else :
            large_list = tile_class.tile_context_sp_single(chr, samname, genomefile,context,window)  
        x = 0
        for list in large_list :
            if list[1] != 0 or list[3] != 0 :
                outfile.write(str(x*int(window))+"\t")
                outfile.write(str((x+1)*int(window)-1)+"\t")
                outfile.write(str(chr)+"\t")
                list = map(str, list)
                outfile.write("\t".join(list))
                outfile.write("\n")
            x += 1
        outfile.close()
        if len(samname.split("/")) > 1 :
            os.chdir("..")
        return outfile                
    
    @staticmethod
    def tile_context_sp_single_merge(samname, chr_list ,context,window) :
        print "merging tile"
        outfile = open(str(samname).strip(".sam")+"_tile_"+str(context)+"_"+str(window)+"bp_new.txt", 'w')
        outfile.write("Start\tEnd\tChrom\tF CpG metC\tF CpG C\tR CpG metC\tR CpG C\n")
        for chr in chr_list :
            for line in open(str(samname).strip(".sam")+"_tile_"+str(context)+"_"+str(window)+"bp_"+str(chr)+".txt", 'r') :
                outfile.write(line)
            #os.remove(str(samname).strip(".sam")+"_tile_"+str(context)+"_"+str(window)+"bp_"+str(chr)+".txt")
        outfile.close()
    
    @staticmethod
    def make_bedGraph(samname, genomefile,context,window, stranded) :
        print "generating bedGraph file"
        outfile = open(str(samname).strip(".sam")+"_tile_"+str(context)+"_"+str(window)+"bp.bedGraph", 'w')
        if stranded == "no" :
            for a in open(str(samname).strip(".sam")+"_tile_"+str(context)+"_"+str(window)+"bp.txt", 'r') :
                if a.split("\t")[0] != "Start" :
                    outfile.write(a.split("\t")[2]+"\t"+a.split("\t")[0]+"\t"+a.split("\t")[1]+"\t"+str(((float(a.split("\t")[3])+float(a.split("\t")[5]))/((float(a.split("\t")[4])+float(a.split("\t")[6].strip("\n"))))))+"\n")
        outfile.close()
        
    @staticmethod
    def get_coverage(samname, context, window, samtools, multibamcov) :
        print "preparing sam to bam"
        #os.system(str(samtools)+" view -bS "+str(samname)+" > "+str(samname).split(".sam")[0]+".bam")
        os.system(str(samtools)+" sort "+str(samname).split(".sam")[0]+".bam "+str(samname).split(".sam")[0]+".sorted")
        os.system(str(samtools)+" index "+str(samname).split(".sam")[0]+".sorted.bam")
        print "running multiBamCov"
        os.system(str(multibamcov)+" -bams "+str(samname).split(".sam")[0]+".sorted.bam -bed "+str(samname).strip(".sam")+"_tile_"+str(context)+"_"+str(window)+"bp.bedGraph > "+str(samname).strip(".sam")+"_tile_"+str(context)+"_"+str(window)+"bp.cov")
    
    @staticmethod
    def circos(samname, context, window, genomefile, specie) :
        specie_dict = {'HG19':'hs', 'B73':'zm'}
        outfile = open(str(samname).strip(".sam")+"_tile_"+str(context)+"_"+str(window)+"bp_new.circos", 'w')
        for a in open(str(samname).strip(".sam")+"_tile_"+str(context)+"_"+str(window)+"bp_new.txt", 'r') :
            if a.split("\t")[0] != "Start" :
                outfile.write(specie_dict[specie]+str(a.split("\t")[2])+"\t")
                outfile.write(str(a.split("\t")[0]))
                outfile.write("\t")
                outfile.write(str(a.split("\t")[1]))
                outfile.write("\t")
                outfile.write(str(((float(a.split("\t")[3])+float(a.split("\t")[5]))/((float(a.split("\t")[4])+float(a.split("\t")[6].strip("\n")))))))
                outfile.write("\n")
        outfile.close()
        
parser = OptionParser()
allgroup = OptionGroup(parser, "Required for all function")
parser.add_option_group(allgroup)
allgroup.add_option("--run", dest="run", help="type of run : tile / methylation_extractor / methylkit / circos")
allgroup.add_option("--genome", dest="genome", help="name and directory of genome fasta file")
allgroup.add_option("--context", dest="context", help="CpG, CHH, CHG or all")
allgroup.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True, help="don't print status messages to stdout")
allgroup.add_option("--specie", dest="specie", help="Specie code, B73, MM9, HG19")
tilegroup = OptionGroup(parser, "Tile specific options","These are tile function specific options\n")
parser.add_option_group(tilegroup)
tilegroup.add_option("--window", dest="window", help="tiling window size", default=100)
tilegroup.add_option("--separate_strands", dest="stranded", help="merge strands=no, separate strands=yes", default="no")
tilegroup.add_option("--path_to_samtools", dest="samtools", help="full path to samtools")
tilegroup.add_option("--path_to_multibamcov", dest="multibamcov", help="full path to multiBamCov")
tilegroup.add_option("--sam", dest="samfile", help="name and directory of first sam file generated for tiling analysis, use unsorted")
tilegroup.add_option("--circos", dest="circos", help="Generate circos formatted file, default=no", default="no")
tilegroup.add_option("--bedGraph", dest="bedGraph", help="Generate bedGraph of tiles, default=no", default="no")
tilegroup.add_option("--coverage", dest="coverage", help="Generate read coverage for each tiles, default=no", default="no")
megroup = OptionGroup(parser, "Methylation extractor specific options","These are Bismark Methylation extractor options\n")
parser.add_option_group(megroup)
megroup.add_option("--paired", dest="paired", help="paired-end or single-end reads, yes or no")
megroup.add_option("--path_to_bismark", dest="bismark", help="full path to bismark")
mkgroup = OptionGroup(parser, "MethylKit specific options","These are MethylKit options\n")
parser.add_option_group(mkgroup)
mkgroup.add_option("--cores", dest="cores", help="no. of cores to use in running BisKit")
mkgroup.add_option("--sam1", dest="samfile1", help="name and directory of first sam file generated from bismark for DMR finding")
mkgroup.add_option("--sam2", dest="samfile2", help="name and directory of second sam file generated from bismark for DMR finding")
mkgroup.add_option("--largemem", dest="largemem", help="If using large memory notes = yes", default="no")
mkgroup.add_option("--nonb73", dest="nonb73", help="For maize genome, if looking at nonB73 genotype", default="no")

(options, args) = parser.parse_args()
    
if options.run == "tile" :
    cores = int(options.cores)
    #total_rounds = len(setup.get_genome_size(options.genome)[0])/cores
    #size_dict = setup.get_genome_size(options.genome)[0]
    if options.nonb73 == "yes" :
        chr_list = ['1','2','3','4','5','6','7','8','9','10']
    else :
        chr_list = setup.get_genome_size(options.genome)[1]
    total_rounds = (len(chr_list)/cores)+1
    for round in range(0, total_rounds) :
        jobs = []
        for chr in chr_list[round*cores:(round+1)*cores] :
            s1 = multiprocessing.Process(target=tile_class.tile_context_sp_single_write, args=(chr, options.samfile, options.genome, options.context, options.window, options.nonb73,))
            jobs.append(s1)
            s1.start()
        [x.join() for x in jobs]
    tile_class.tile_context_sp_single_merge(options.samfile, chr_list, options.context, options.window)
    if options.circos == "yes" :
        tile_class.circos(options.samfile, options.context, options.window, options.genome, options.specie)
    if options.bedGraph == "yes" :
        tile_class.make_bedGraph(options.samfile, options.genome, options.context, options.window, options.stranded)
    if options.coverage == "yes" :
        tile_class.get_coverage(options.samfile, options.context, options.window, options.samtools, options.multibamcov)
    print "complete" ; sys.exit()

elif options.run == "methylkation_extractor" :
    if options.samfile2 != None and options.cores > 1 :
        jobs = []
        s1 = multiprocessing.Process(target=bismark.methylation_extractor, args=(options.bismark, options.samfile1, options.paired, ))
        s2 = multiprocessing.Process(target=bismark.methylation_extractor, args=(options.bismark, options.samfile2, options.paired, ))
        jobs.append(s1)
        s1.start()
        jobs.append(s2)
        s2.start()
        [x.join() for x in jobs]
    elif options.samfile2 != None and options.core == 1 :
        bismark.methylation_extractor(options.bismark, options.samfile1, options.paired) 
        bismark.methylation_extractor(options.bismark, options.samfile2, options.paired) 
    else :
        bismark.methylation_extractor(options.bismark, options.samfile1, options.paired) 

elif options.run == "methylkit" :
    cores = int(options.cores)
    chr_list = {}
    for chr in setup.get_genome_size(options.genome)[0].keys() :
        if chr[:3] == "chr" : 
            chr_list[str(chr).split("chr")[1]] = []
        else :
            chr_list[str(chr)] = []
    rand_chr = 100
    for chr in chr_list.keys() :
        if chr.isdigit() == True :
            chr_list[chr].append(int(chr))
        else :
            chr_list[chr].append(int(rand_chr))
            rand_chr += 1
    if options.nonb73 == "yes" :
        chr_list = {'1':[1],'2':[2],'3':[3],'4':[4],'5':[5],'6':[6],'7':[7],'8':[8],'9':[9],'10':[10]} 
    total_rounds = (len(chr_list)/cores)+1
    """
    for round in range(0, total_rounds) :
        jobs = []
        for chr in chr_list.keys()[round*cores:(round+1)*cores] :
            print chr, chr_list[chr][0], options.samfile1, options.context, options.largemem, options.nonb73
            if options.nonb73 == "yes" :
                #s1 = multiprocessing.Process(target=methylkit.prep_methylkit_nonb73, args=(chr, chr_list[chr][0], options.samfile1, options.context, options.largemem, ))
                #jobs.append(s1)
                #s1.start()
                s2 = multiprocessing.Process(target=methylkit.prep_methylkit_nonb73, args=(chr, chr_list[chr][0], options.samfile2, options.context, options.largemem, ))
                jobs.append(s2)
                s2.start()
            else :
                s1 = multiprocessing.Process(target=methylkit.prep_methylkit, args=(chr, chr_list[chr][0], options.samfile1, options.context, options.largemem, ))
                jobs.append(s1)
                s1.start()
                s2 = multiprocessing.Process(target=methylkit.prep_methylkit, args=(chr, chr_list[chr][0], options.samfile2, options.context, options.largemem, ))
                jobs.append(s2)
                s2.start()
        [x.join() for x in jobs]
    print datetime.now()
    for round in range(0, total_rounds) :
        jobs = []
        for chr in chr_list.keys()[round*cores:(round+1)*cores] :
            s = multiprocessing.Process(target=methylkit.methylkit_edmr, args=(options.samfile1, options.samfile2, options.context, chr, options.specie,  ))
            jobs.append(s)
            s.start()
        [x.join() for x in jobs]
    print datetime.now()
    print "complete" ; sys.exit() 
    """
    for round in range(0, total_rounds) :
        jobs = []
        for chr in chr_list.keys()[round*cores:(round+1)*cores] :
            s = multiprocessing.Process(target=methylkit.rewrite_dmr, args=(options.samfile1, options.samfile2, options.context, chr, options.specie,  ))
            jobs.append(s)
            s.start()
        [x.join() for x in jobs]
    print datetime.now()
    methylkit.merge_dmr(options.samfile1, options.samfile2, options.context, options.specie)
    print "complete" ; sys.exit() 
    
elif options.run == "speed" :
    methylkit.speed_test(options.samfile1, options.samfile2, options.context, options.specie)
