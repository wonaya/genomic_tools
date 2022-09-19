#!/usr/bin/env python

### add option for normal mem

import os,sys
import resource
from optparse import OptionParser
from datetime import datetime
try:
   import subprocess
except ImportError:
   print >> sys.stderr,"Could not import the subprocess module"
   raise

from biskit import *
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
                chr_list.append(genome[1:].strip("\n"))
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
    def paired_alignment(outdir, fastq1, fastq2, paired) :
        os.system("src/bismark /work/02114/wonaya/genome/Zea_mays.AGPv2.14 -q -N 1 -p 3 --bowtie2 --path_to_bowtie=/work/02114/wonaya/software/bowtie2-2.0.2 -o "+str(outdir)+" -1 "+str(fastq1)+" -2 "+str(fastq2))
    @staticmethod
    def single_alignment(outdir, fastq1) :
        os.system("src/bismark /work/02114/wonaya/genome/Zea_mays.AGPv2.14 -q -N 1 -p 3 --bowtie2 --path_to_bowtie=/work/02114/wonaya/software/bowtie2-2.0.2 -o "+str(outdir)+" "+str(fastq1))
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
    def prep_methylkit(chr, write_chr, samname_ori, context, largemem) :
        if len(samname_ori.split("/")) > 1 :
            os.chdir(samname_ori.split("/")[0])
            samname = samname_ori.split("/")[1]
        else :
            samname = samname_ori
        if not os.path.isfile(str(context)+"_"+str(samname)+"_chr"+str(chr)+".methylKit2") :
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
                print samname, x
                large_list = []
                for y in range(0,max_list): ### 56 percent of lonestar node, 2min to make
                    large_list.append([0]*4)
                for a_f in open(str(context)+"_OT_"+str(samname).split(".sam")[0]+".txt", 'r') : # example: CpG_OT_11510-pool_GGTAGC_L004_R1_001_val_1.fq_bismark_bt2_pe.txt
                    if len(a_f.split("\t")) > 2 and a_f.split("\t")[2] == "chr"+str(chr) and len(a_f.split("\t")) != 0 and int(a_f.split("\t")[3]) >= x*max_list and int(a_f.split("\t")[3]) < (x+1)*max_list :
                        if a_f.split("\t")[1] == "+" :
                            large_list[int(float(a_f.split("\t")[3]))-(x*max_list)][0] += 1
                        elif a_f.split("\t")[1] == "-" :
                            large_list[int(float(a_f.split("\t")[3]))-(x*max_list)][1] += 1
                if os.path.isfile(str(context)+"_CTOT_"+str(samname).split(".sam")[0]+".txt") :
                    for a_f in open(str(context)+"_CTOT_"+str(samname).split(".sam")[0]+".txt", 'r') : # example: 
                        if len(a_f.split("\t")) > 2 and a_f.split("\t")[2] == "chr"+str(chr) and len(a_f.split("\t")) != 0 and int(a_f.split("\t")[3]) >= x*max_list and int(a_f.split("\t")[3]) < (x+1)*max_list :
                            if a_f.split("\t")[1] == "+" :
                                large_list[int(float(a_f.split("\t")[3]))-(x*max_list)][0] += 1
                            elif a_f.split("\t")[1] == "-" :
                                large_list[int(float(a_f.split("\t")[3]))-(x*max_list)][1] += 1
                
                for a_r in open(str(context)+"_OB_"+str(samname).split(".sam")[0]+".txt", 'r') : 
                    if len(a_r.split("\t")) > 2 and a_r.split("\t")[2] == "chr"+str(chr) and len(a_r.split("\t")) != 0 and int(a_r.split("\t")[3]) >= x*max_list and int(a_r.split("\t")[3]) < (x+1)*max_list :
                        if a_r.split("\t")[1] == "+" :
                            large_list[int(float(a_r.split("\t")[3]))-(x*max_list)][2] += 1
                        elif a_r.split("\t")[1] == "-" :
                            large_list[int(float(a_r.split("\t")[3]))-(x*max_list)][3] += 1
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
        #if len(samname_ori.split("/")) > 1 :
        #    os.chdir(samname_ori.split("/")[0])
        #    samname = samname_ori.split("/")[1]
        #else :
        #    samname = samname_ori
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
        #os.system("mv "+str(file1ID)+"_"+str(file2ID)+"_"+str(type)+"_dmr_chr"+str(chr)+"_temp.txt "+str(file1ID)+"_"+str(file2ID)+"_"+str(type)+"_dmr_chr"+str(chr)+".txt")
    
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
            for line1 in open(str(context)+"_OT_"+tile_class.extract_name(samname, genomefile)[0]+".txt", 'r') :
                if line1[0:7] != "Bismark" and line1.split("\t")[2] == str(chr) :
                    if line1.split("\t")[1] == "+" :
                        large_list[int(float(line1.split("\t")[3])/int(window))][0] += 1 
                        large_list[int(float(line1.split("\t")[3])/int(window))][1] += 1
                    elif line1.split("\t")[1] == "-" :
                        large_list[int(float(line1.split("\t")[3])/int(window))][1] += 1
            for line2 in open(str(context)+"_OB_"+tile_class.extract_name(samname, genomefile)[0]+".txt", 'r') :
                if line2[0:7] != "Bismark" and line2.split("\t")[2] == str(chr):
                    if line2.split("\t")[1] == "+" :
                        large_list[int(float(line2.split("\t")[3])/int(window))][2] += 1 
                        large_list[int(float(line2.split("\t")[3])/int(window))][3] += 1
                    elif line2.split("\t")[1] == "-" :
                        large_list[int(float(line2.split("\t")[3])/int(window))][3] += 1
        return large_list
    
    @staticmethod
    def tile_context_sp_single_write(chr, samname, genomefile,context,window) :
        outfile = open(str(samname)+"_tile_"+str(context)+"_"+str(window)+"bp_"+str(chr)+".txt", 'w')
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
        return outfile                
    
    @staticmethod
    def tile_context_sp_single_merge(samname, genomefile,context,window) :
        print "merge"
        outfile = open(str(samname)+"_tile_"+str(context)+"_"+str(window)+"bp.txt", 'w')
        outfile.write("Start\tEnd\tChrom\tF CpG metC\tF CpG C\tR CpG metC\tR CpG C\n")
        for chr in setup.get_genome_size(genomefile)[1]:
            for line in open(str(samname)+"_tile_"+str(context)+"_"+str(window)+"bp_"+str(chr)+".txt", 'r') :
                outfile.write(line)
            os.remove(str(samname)+"_tile_"+str(context)+"_"+str(window)+"bp_"+str(chr)+".txt")
        outfile.close()
    
    @staticmethod
    def make_bedGraph(samname, genomefile,context,window, stranded) :
        outfile = open(str(samname)+"_tile_"+str(context)+"_"+str(window)+"bp.bedGraph", 'w')
        if stranded == "no" :
            for a in open(str(samname)+"_tile_"+str(context)+"_"+str(window)+"bp.txt", 'r') :
                if a.split("\t")[0] != "Start" :
                    outfile.write(a.split("\t")[2]+"\t"+a.split("\t")[0]+"\t"+a.split("\t")[1]+"\t"+str(((float(a.split("\t")[3])+float(a.split("\t")[5]))/((float(a.split("\t")[4])+float(a.split("\t")[6].strip("\n"))))))+"\n")
        outfile.close()
        
    @staticmethod
    def get_coverage(samname, context, window, samtools, multibamcov) :
        os.system(str(samtools)+" view -bS "+str(samname)+" > "+str(samname).split(".sam")[0]+".bam")
        os.system(str(samtools)+" sort "+str(samname).split(".sam")[0]+".bam "+str(samname).split(".sam")[0]+".sorted")
        os.system(str(samtools)+" index "+str(samname).split(".sam")[0]+".sorted.bam")
        os.system(str(multibamcov)+" -bams "+str(samname).split(".sam")[0]+".sorted.bam -bed "+str(samname)+"_tile_"+str(context)+"_"+str(window)+"bp.bedGraph > "+str(samname)+"_tile_"+str(context)+"_"+str(window)+"bp.cov")

parser = OptionParser()
parser.add_option("--run", dest="run", help="type of run, tile / methylation_extractor / methylkit")
parser.add_option("--genome", dest="genome", help="name and directory of genome fasta file")
parser.add_option("--context", dest="context", help="CpG, CHH, CHG or all")
parser.add_option("--window", dest="window", help="tiling window size", default=100)
parser.add_option("--separate_strands", dest="stranded", help="merge strands=no, separate strands=yes", default="no")
parser.add_option("--path_to_samtools", dest="samtools", help="full path to samtools")
parser.add_option("--path_to_multibamcov", dest="multibamcov", help="full path to multiBamCov")
parser.add_option("--path_to_bismark", dest="bismark", help="full path to bismark")
parser.add_option("--paired", dest="paired", help="paired-end or single-end reads, yes or no")
parser.add_option("--cores", dest="cores", help="no. of cores to use in running BisKit")
parser.add_option("--sam1", dest="samfile1", help="name and directory of first sam file generated from bismark for DMR finding")
parser.add_option("--sam2", dest="samfile2", help="name and directory of second sam file generated from bismark for DMR finding")
parser.add_option("--specie", dest="specie", help="Specie code, B73, MM9, HG19")
parser.add_option("--largemem", dest="largemem", help="If using large memory notes = yes", default="no")

(options, args) = parser.parse_args()
if options.run == "tile" :
    cores = int(options.cores)
    total_rounds = len(setup.get_genome_size(options.genome)[0])/cores
    for round in range(0, total_rounds) :
        jobs = []
        for chr in setup.get_genome_size(options.genome)[0].keys()[round*cores:(round+1)*cores] :
            s1 = multiprocessing.Process(target=tile_class.tile_context_sp_single_write, args=(chr, options.samfile1, options.genome, options.context, options.window, ))
            jobs.append(s1)
            s1.start()
        [x.join() for x in jobs]
    tile_class.tile_context_sp_single_merge(options.samfile1, options.genome, options.context, options.window)
    tile_class.make_bedGraph(options.samfile1, options.genome, options.context, options.window, options.stranded)
    tile_class.get_coverage(options.samfile1, options.context, options.window, options.samtools, options.multibamcov)
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
    print datetime.now()
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
    total_rounds = len(chr_list)/cores
    print datetime.now()
    
    for round in range(0, total_rounds) :
        jobs = []
        #for chr in chr_list.keys()[round*cores:(round+1)*cores] :
        for chr in [1,2,3]:
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
        #for chr in setup.get_genome_size(options.genome)[0].keys()[round*cores:(round+1)*cores] :
        #for chr in chr_list.keys()[round*cores:(round+1)*cores] :
        for chr in [1,2,3]:
            s = multiprocessing.Process(target=methylkit.methylkit_edmr, args=(options.samfile1, options.samfile2, options.context, chr, options.specie,  ))
            jobs.append(s)
            s.start()
        [x.join() for x in jobs]
    print datetime.now()
    for round in range(0, total_rounds) :
        jobs = []
        #for chr in setup.get_genome_size(options.genome)[0].keys()[round*cores:(round+1)*cores] :
        #for chr in chr_list.keys()[round*cores:(round+1)*cores] :
        for chr in [1,2,3]:
            s = multiprocessing.Process(target=methylkit.rewrite_dmr, args=(options.samfile1, options.samfile2, options.context, chr, options.specie,  ))
            jobs.append(s)
            s.start()
        [x.join() for x in jobs]
    print datetime.now()
    print "complete" ; sys.exit() 
    
    
elif options.run == "merge" :
    methylkit.merge_dmr(options.samfile1, options.samfile2, options.context, options.specie)
    
    