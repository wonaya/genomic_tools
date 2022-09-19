#### written by Jawon Song
#### updated 12-11-2013
#### updated 12-17-2013

import os,sys
import datetime
import multiprocessing
import readline, glob
import subprocess
from optparse import OptionParser

def complete(text, state):
    return (glob.glob(text+'*')+[None])[state]

readline.set_completer_delims(' \t\n;')
readline.parse_and_bind("tab: complete")
readline.set_completer(complete)

class setup :
    def install_bismark() :
        print "copy over binary"
        
class bismark :
    @staticmethod
    def paired_alignment(outdir, fastq1, fastq2) :
        os.system("/opt/apps/bismark/0.10.0/bismark /work/02114/wonaya/genome/Zea_mays.AGPv2.14 -q -N 1 -p 3 --bowtie2 --path_to_bowtie=/work/02114/wonaya/software/bowtie2-2.0.2 -o "+str(outdir)+" -1 "+str(fastq1)+" -2 "+str(fastq2))
    @staticmethod
    def single_alignment(outdir, fastq1) :
        os.system("/opt/apps/bismark/0.10.0/bismark /work/02114/wonaya/genome/Zea_mays.AGPv2.14 -q -N 1 -p 3 --bowtie2 --path_to_bowtie=/work/02114/wonaya/software/bowtie2-2.0.2 -o "+str(outdir)+" "+str(fastq1))
    
    def methylation_extractor() :
        print "do methylation extractor"
          
class tile_class :
    @staticmethod
    def extract_name(filename) :
        ## given Bam file used for extraction get file ID for subsequent analysis
        dir = "."
        if len(filename.split("/")) > 1 :
            dir = filename.split("/")[0]
            filename = filename.split("/")[1]
        fieid = filename.split(".bam")[0]
        return fieid, dir
    
    @staticmethod
    def tile(name, chr, window) :
        ## define genotype and chromosome
        ## make array of 8 columns and put 0
        ## CpG/TotC/CHG/TotC/CHH/TotC/all/TotC forward (first 8 columns) and reverse
        ## making array
        large_list = []
        for x in range(0,500000000/int(window)):
            large_list.append([0]*16)

        ### CpG
        ### forward strand
        for CpG_OT in open("CpG_OT_"+str(name)+".txt", 'r') :
            if CpG_OT[0:7] != "Bismark" :
                if CpG_OT.split("\t")[2] == "chr"+str(chr) :
                    if CpG_OT.split("\t")[1] == "+" :
                        large_list[int(float(CpG_OT.split("\t")[3])/int(window))][0] += 1 
                        large_list[int(float(CpG_OT.split("\t")[3])/int(window))][1] += 1
                    elif CpG_OT.split("\t")[1] == "-" :
                        large_list[int(float(CpG_OT.split("\t")[3])/int(window))][1] += 1
        del CpG_OT
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
        for CpG_OB in open("CpG_OB_"+str(name)+".txt", 'r') :
            if CpG_OB[0:7] != "Bismark" :
                if CpG_OB.split("\t")[2] == "chr"+str(chr) :
                    if CpG_OB.split("\t")[1] == "+" :
                        large_list[int(float(CpG_OB.split("\t")[3])/int(window))][8] += 1 
                        large_list[int(float(CpG_OB.split("\t")[3])/int(window))][9] += 1
                    elif CpG_OB.split("\t")[1] == "-" :
                        large_list[int(float(CpG_OB.split("\t")[3])/int(window))][9] += 1
        del CpG_OB
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
        for CHG_OT in open("CHG_OT_"+str(name)+".txt", 'r') :
            if CHG_OT[0:7] != "Bismark" :
                if CHG_OT.split("\t")[2] == "chr"+str(chr) :
                    if CHG_OT.split("\t")[1] == "+" :
                        large_list[int(float(CHG_OT.split("\t")[3])/int(window))][2] += 1 
                        large_list[int(float(CHG_OT.split("\t")[3])/int(window))][3] += 1
                    elif CHG_OT.split("\t")[1] == "-" :
                        large_list[int(float(CHG_OT.split("\t")[3])/int(window))][3] += 1
        del CHG_OT
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
        for CHG_OB in open("CHG_OB_"+str(name)+".txt", 'r') :
            if CHG_OB[0:7] != "Bismark" :
                if CHG_OB.split("\t")[2] == "chr"+str(chr) :
                    if CHG_OB.split("\t")[1] == "+" :
                        large_list[int(float(CHG_OB.split("\t")[3])/int(window))][10] += 1 
                        large_list[int(float(CHG_OB.split("\t")[3])/int(window))][11] += 1
                    elif CHG_OB.split("\t")[1] == "-" :
                        large_list[int(float(CHG_OB.split("\t")[3])/int(window))][11] += 1
        del CHG_OB
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
        for CHH_OT in open("CHH_OT_"+str(name)+".txt", 'r') :
            if CHH_OT[0:7] != "Bismark" :
                if CHH_OT.split("\t")[2] == "chr"+str(chr) :
                    if CHH_OT.split("\t")[1] == "+" :
                        large_list[int(float(CHH_OT.split("\t")[3])/int(window))][4] += 1 
                        large_list[int(float(CHH_OT.split("\t")[3])/int(window))][5] += 1
                    elif CHH_OT.split("\t")[1] == "-" :
                        large_list[int(float(CHH_OT.split("\t")[3])/int(window))][5] += 1
        del CHH_OT
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
        for CHH_OB in open("CHH_OB_"+str(name)+".txt", 'r') :
            if CHH_OB[0:7] != "Bismark" :
                if CHH_OB.split("\t")[2] == "chr"+str(chr) :
                    if CHH_OB.split("\t")[1] == "+" :
                        large_list[int(float(CHH_OB.split("\t")[3])/int(window))][12] += 1 
                        large_list[int(float(CHH_OB.split("\t")[3])/int(window))][13] += 1
                    elif CHH_OB.split("\t")[1] == "-" :
                        large_list[int(float(CHH_OB.split("\t")[3])/int(window))][13] += 1
        del CHH_OB
        for CHH_CTOB in open("CHH_CTOB_"+str(name)+".txt", 'r') :
            if CHH_CTOB[0:7] != "Bismark" :
                if CHH_CTOB.split("\t")[2] == "chr"+str(chr) :
                    if CHH_CTOB.split("\t")[1] == "+" :
                        large_list[int(float(CHH_CTOB.split("\t")[3])/int(window))][12] += 1 
                        large_list[int(float(CHH_CTOB.split("\t")[3])/int(window))][13] += 1
                    elif CHH_CTOB.split("\t")[1] == "-" :
                        large_list[int(float(CHH_CTOB.split("\t")[3])/int(window))][13] += 1
        del CHH_CTOB

        ## all metC
        for list in large_list :
            list[6] += list[0]+list[2]+list[4]
            list[7] += list[1]+list[3]+list[5]
            list[14] += list[8]+list[10]+list[12]
            list[15] += list[9]+list[11]+list[13]
        #print "write to file", "chr"+str(chr),datetime.datetime.now()
        return large_list
    
    @staticmethod
    def write_to_file(chr, name, window) :
        print "chr", chr
        outfile = open(str(name)+"_tile_"+str(window)+"bp_chr"+str(chr)+".txt", 'w')
        outfile.write("Start\tEnd\tChrom\tF CpG metC\tF CpG C\tF CHG metC\tF CHG C\tF CHH metC\tF CHH C\tF all metC\tF all C\tR CpG metC\tR CpG C\tR CHG metC\tR CHG C\tR CHH metC\tR CHH C\tR all metC\tR all C\n")
        large_list = tile_class.tile(name, chr, window)   
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
        return outfile    
    
    @staticmethod
    def tile_merge(name, window) :
        outfile = open(str(name)+"_tile_merged_"+str(window)+"bp.txt", 'w')
        outfile.write("Start\tEnd\tChrom\tF CpG metC\tF CpG C\tF CHG metC\tF CHG C\tF CHH metC\tF CHH C\tF all metC\tF all C\tR CpG metC\tR CpG C\tR CHG metC\tR CHG C\tR CHH metC\tR CHH C\tR all metC\tR all C\n")
        for chr in range(1,11) :
            print "writing chr"+str(chr)
            for a in open(str(name)+"_tile_"+str(window)+"bp_chr"+str(chr)+".txt", 'r') :
                if a.split("\t")[0] != "Start" :
                    outfile.write(a)
        outfile.close()
        os.system("rm -Rf "+str(name)+"_tile_"+str(window)+"bp_chr*")
    
    @staticmethod
    def tile_flow(name, window) :
        name = tile_class.extract_name(sys.argv[1])[0]
        dir = tile_class.extract_name(sys.argv[1])[1]
        os.chdir(dir)

        jobs = []
        chrmax = 11
        for chr in range(1,chrmax) :
            s = multiprocessing.Process(target=tile_class.write_to_file, args=(chr, name, 1000, ))
            jobs.append(s)
        [x.start() for x in jobs]
        [x.join() for x in jobs]

        tile_class.tile_merge(name, 1000)

class methylkit:
    @staticmethod
    def prep_methylKit2(chr, type, name, dir) :
        os.chdir(dir)
        if not os.path.isfile(str(type)+"_"+str(name)+"_chr"+str(chr)+".methylKit") :
            outfile = open(str(type)+"_"+str(name)+"_chr"+str(chr)+".methylKit" ,'w')
            ## save strand status
            print chr, type, name, dir
            for x in range(0, 4) : ## max set to 4 0 to 399Mbp
                large_list = []
                for y in range(0,100000000): ### 56 percent of lonestar node, 2min to make
                    large_list.append([0]*4)
                print x*100000000, "to", (x+1)*100000000, datetime.datetime.now(), ": writing array"
                for a_f in open(str(type)+"_OT_"+str(name)+".txt", 'r') : # example: CpG_OT_11510-pool_GGTAGC_L004_R1_001_val_1.fq_bismark_bt2_pe.txt
                    if len(a_f.split("\t")) > 2 and a_f.split("\t")[2] == "chr"+str(chr) and len(a_f.split("\t")) != 0 and int(a_f.split("\t")[3]) >= x*100000000 and int(a_f.split("\t")[3]) < (x+1)*100000000 :
                        if a_f.split("\t")[1] == "+" :
                            large_list[int(float(a_f.split("\t")[3]))-(x*100000000)][0] += 1
                        elif a_f.split("\t")[1] == "-" :
                            large_list[int(float(a_f.split("\t")[3]))-(x*100000000)][1] += 1
                for a_f in open(str(type)+"_CTOT_"+str(name)+".txt", 'r') : # example: 
                    if len(a_f.split("\t")) > 2 and a_f.split("\t")[2] == "chr"+str(chr) and len(a_f.split("\t")) != 0 and int(a_f.split("\t")[3]) >= x*100000000 and int(a_f.split("\t")[3]) < (x+1)*100000000 :
                        if a_f.split("\t")[1] == "+" :
                            large_list[int(float(a_f.split("\t")[3]))-(x*100000000)][0] += 1
                        elif a_f.split("\t")[1] == "-" :
                            large_list[int(float(a_f.split("\t")[3]))-(x*100000000)][1] += 1
        
                for a_r in open(str(type)+"_OB_"+str(name)+".txt", 'r') : 
                    if len(a_r.split("\t")) > 2 and a_r.split("\t")[2] == "chr"+str(chr) and len(a_r.split("\t")) != 0 and int(a_r.split("\t")[3]) >= x*100000000 and int(a_r.split("\t")[3]) < (x+1)*100000000 :
                        if a_r.split("\t")[1] == "+" :
                            large_list[int(float(a_r.split("\t")[3]))-(x*100000000)][2] += 1
                        elif a_r.split("\t")[1] == "-" :
                            large_list[int(float(a_r.split("\t")[3]))-(x*100000000)][3] += 1
                for a_r in open(str(type)+"_CTOB_"+str(name)+".txt", 'r') : 
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
        os.chdir("..")
        
    @staticmethod
    def methylkit_rpy(file1, file2, type, chr) :
        try :
            from rpy2.robjects.packages import importr
        except :
            #print subprocess.Popen("module load R", stdout=subprocess.PIPE, shell=True, executable="/bin/bash").stdout.read()
            print "do module load R"
            sys.exit()
        import rpy2.robjects as robjects
        importr('GenomicRanges')
        importr('bigmemory')
        importr('data.table')
        robjects.r('''source('/work/02114/wonaya/scripts/base.r')''')
        filelist = [file1, file2]
        file1ID = filelist[0].split("_")[0]
        file2ID = filelist[1].split("_")[0]
        robjects.r.assign('rfilelist',filelist)
        robjects.r.assign('type', type)
        robjects.r.assign('file1ID', file1ID)
        robjects.r.assign('file2ID', file2ID)
        outfile = str(file1ID)+"_"+str(file2ID)+"_"+str(type)+"_dmr_chr"+str(chr)+".txt"
        robjects.r.assign('outfile', outfile)
        robjects.r('''myobj=read(rfilelist, sample.id=list(file1ID,file2ID), assembly="B73",treatment=c(0,1),context=type, resolution="base")''')
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
        robjects.r('''myDMR=eDMR(myDiff,step=100,dist="none",DMC.qvalue=0.01,DMC.methdiff=25,num.DMCs=1,num.CpGs=3,DMR.methdiff=20,plot=FALSE,main="example",mode=1)''')
        robjects.r('''myDMR.sig=filter.dmr(myDMR, DMR.qvalue=0.05, mean.meth.diff=20, num.CpGs=3, num.DMCs=1)''')
        robjects.r('''myDMR.sig.df <- as.data.frame(myDMR.sig)''')
        robjects.r('''write.table(myDMR.sig.df, file=outfile, sep="\t")''')
    
    @staticmethod
    def dummy(file1, file2, type, chr) :
        print file1, file2, type, chr
        
    @staticmethod
    def methylkit_flow(type) : ## give two bam file names as input
        #subprocess.call("module load R", shell=True, executable="/bin/bash")
        #subprocess.call("module list", shell=True, executable="/bin/bash")
        try :
            from rpy2.robjects.packages import importr
        except :
            print "do module load R"
            sys.exit()
        
        name1 = tile_class.extract_name(sys.argv[2])[0]
        name2 = tile_class.extract_name(sys.argv[3])[0]
        dir1 = tile_class.extract_name(sys.argv[2])[1]
        dir2 = tile_class.extract_name(sys.argv[3])[1]
        names = [name1, name2]
        dirs = [dir1,dir2]
        print names, dirs, datetime.datetime.now()
        chrmax = 11
        jobs = []
        for chr in range(1,chrmax) :
            for name in names :
                s = multiprocessing.Process(target=methylkit.prep_methylKit2, args=(chr, type, name, dirs[names.index(name)], ))
                jobs.append(s)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
        print "prep methylkit done", datetime.datetime.now()
        jobs = []
        for chr in range(1,11) :
            file1 = str(dir1)+"/"+str(type)+"_"+str(name1)+"_chr"+str(chr)+".methylKit"
            file2 = str(dir2)+"/"+str(type)+"_"+str(name2)+"_chr"+str(chr)+".methylKit"
            s = multiprocessing.Process(target=methylkit.methylkit_rpy, args=(file1, file2, type, chr, ))
            jobs.append(s)
        [x.start() for x in jobs]
        [x.join() for x in jobs]
        print "methylkit done", datetime.datetime.now()
        name1 = name1.split("_")[0]
        name2 = name2.split("_")[0]
        outfile = open(str(name1)+"_"+str(name2)+"_"+str(type)+"_dmr_merged.txt", 'w')
        for chr in range(1,chrmax) :
            for line in open(str(name1)+"_"+str(name2)+"_"+str(type)+"_dmr_chr"+str(chr)+".txt", 'r') :
                if line.split("\t")[0] != '"seqnames"' : 
                    outfile.write(line)
        outfile.close()
        print "merge done", datetime.datetime.now()
        
class analysis :
    @staticmethod
    def make_circos_file(bedgraph) :
        """
        name = str(bedgraph).split(".")[0]
        outfile = open(name+".circos", 'w')
        for a in open(bedgraph, 'r') :
            outfile.write("zm"+a.split("\t")[0]+"\t"+"\t".join(a.split("\t")[1:]))
        outfile.close()
        sys.exit()
        """
    def run_circos_file(bedgraph) :
        conf_file = open("/work/02114/wonaya/software/circos-0.64/circos_test/tmp.conf", 'w')
        conf_file.write("<<include etc/colors_fonts_patterns.conf>>\n")
        conf_file.write("<<include ideogram.conf>>\n")
        conf_file.write("<<include ticks.conf>>\n")
        conf_file.write("<image>\n")
        conf_file.write("<<include etc/image.conf>>\n")
        conf_file.write("</image>\n")
        conf_file.write("karyotype   = data/karyotype/karyotype.maize.txt\n")
        conf_file.write("chromosomes_units  = 1000000\n")
        conf_file.write("chromosomes        = zm10\n")
        conf_file.write("chromosomes_display_default = yes\n")
        conf_file.write("<plots>\n")
        conf_file.write("type    = heatmap\n")
        conf_file.write("color = set3-9-qua\n")
        conf_file.write("stroke_thickness = 1\n")
        conf_file.write("stroke_color     = black\n")
        conf_file.write("<plot>\n")
        conf_file.write("file = /scratch/02114/wonaya/UMN_Bismark/10-01-13_5genos/Mo17_CML322_array_dmr.circos\n")
        conf_file.write("r1 = 0.93r\n")
        conf_file.write("r0 = 0.98r\n")
        conf_file.write("color = rdbu-9-div\n")
        conf_file.write("stroke_thickness = 0\n")
        conf_file.write("min = -2\n")
        conf_file.write("max = 2\n")
        conf_file.write("</plot>\n")
        conf_file.write("<plot>\n")
        conf_file.write("file = /scratch/02114/wonaya/UMN_Bismark/10-01-13_5genos/Mo17_all3_CML322_all3_CpG_dmr_merged.circos\n")
        conf_file.write("r1 = 0.87r\n")
        conf_file.write("r0 = 0.92r\n")
        conf_file.write("color = rdbu-9-div\n")
        conf_file.write("stroke_thickness = 0\n")
        conf_file.write("min = -100\n")
        conf_file.write("max = 100\n")
        conf_file.write("</plot>\n")
        conf_file.write("<plot>\n")
        conf_file.write("file = /scratch/02114/wonaya/UMN_Bismark/10-01-13_5genos/Mo17_all3_CML322_all3_CHG_dmr_merged.circos\n")
        conf_file.write("r1 = 0.81r\n")
        conf_file.write("r0 = 0.86r\n")
        conf_file.write("color = rdbu-9-div\n")
        conf_file.write("stroke_thickness = 0\n")
        conf_file.write("min = -100\n")
        conf_file.write("max = 100\n")
        conf_file.write("</plot>\n")
        conf_file.write("<plot>\n")
        conf_file.write("file = /scratch/02114/wonaya/UMN_Bismark/10-01-13_5genos/Mo17_all3_CML322_all3_CHH_dmr_merged.circos\n")
        conf_file.write("r1 = 0.75r\n")
        conf_file.write("r0 = 0.80r\n")
        conf_file.write("color = rdbu-9-div\n")
        conf_file.write("stroke_thickness = 0\n")
        conf_file.write("min = -100\n")
        conf_file.write("max = 100\n")
        conf_file.write("</plot>\n")
        conf_file.write("</plots>\n")
        conf_file.write("<<include etc/housekeeping.conf>>\n")
        conf_file.close()
        os.system("/work/02114/wonaya/software/circos-0.64/bin/circos -conf /work/02114/wonaya/software/circos-0.64/circos_test/tmp.conf")
        os.system("mv circos.png Mo17_CML322.png")
        os.system("rm -Rf circos.svg")
        os.system("rm -Rf /work/02114/wonaya/software/circos-0.64/circos_test/tmp.conf")
        
    @staticmethod
    def make_bedGraph(dmrfile) :
        outfile = open(str(dmrfile).split(".")[0]+".bedGraph", 'w')
        for a in open(dmrfile, 'r') :
            if a.split("\t")[0] != "chrom" :
                if isinstance(a.split("\t")[0][0], int) == False :
                    #outfile.write(str(a.split("\t")[0][3:])+"\t")
                    outfile.write(str(a.split("\t")[1].strip('"'))+"\t")
                else :
                    print "update when chrom is written as integer only"
                    sys.exit()
                outfile.write(a.split("\t")[2]+"\t")
                outfile.write(a.split("\t")[3]+"\t")
                ### currently set to use score of mean. methylation difference
                outfile.write(a.split("\t")[6]+"\n")
        outfile.close()
    
    @staticmethod
    def commondmr(bedgraph1, bedgraph2) :
        ## output common and two different dmr bedgraphs 
        ## and percentage overlap by type
        no_of_overlap = 0
        common_file1 = open(str(bedgraph1).split(".")[0]+"_common.bedGraph", 'w') 
        #common_file2 = open(str(bedgraph2).split(".")[0]+"_common.bedGraph", 'w') 
        remain1 = open(str(bedgraph1).split(".")[0]+"_remainder.bedGraph", 'w') 
        #remain2 = open(str(bedgraph2).split(".")[0]+"_remainder.bedGraph", 'w') 
        for chr in range(1,11) :
            list_dmr_1 = []
            val_dmr_1  = []
            for a in open(bedgraph1, 'r') :
                if int(a.split("\t")[0]) == chr :
                    list_dmr_1.append(range(int(a.split("\t")[1]),int(a.split("\t")[2])+1))
                    val_dmr_1.append(a.split("\t")[3].strip("\n"))
            list_dmr_2 = []
            val_dmr_2  = []
            for b in open(bedgraph2, 'r') :
                if int(b.split("\t")[0]) == chr :    
                    list_dmr_2.append(range(int(b.split("\t")[1]),int(b.split("\t")[2])+1))
                    val_dmr_2.append(b.split("\t")[3].strip("\n"))
            common_dmr = []
            common_dmr_from_1 = []
            common_dmr_from_2 = []
            common_val = []
            for dmr_1 in list_dmr_1 :
                for dmr_2 in list_dmr_2 :
                    if len(set(dmr_1) & set(dmr_2)) > 0 :
                        commons = []
                        common_dmr_from_1.append(dmr_1)
                        common_dmr_from_2.append(dmr_2)
                        common_val.append([val_dmr_1[list_dmr_1.index(dmr_1)],val_dmr_2[list_dmr_2.index(dmr_2)]])
                        for dmr1 in dmr_1 :
                            if dmr1 not in commons :
                                commons.append(dmr1)
                        for dmr2 in dmr_2 :
                            if dmr2 not in commons :
                                commons.append(dmr2)
                        common_dmr.append(commons)
                        val_dmr_1.remove(val_dmr_1[list_dmr_1.index(dmr_1)])
                        list_dmr_1.remove(dmr_1)
                        val_dmr_2.remove(val_dmr_2[list_dmr_2.index(dmr_2)])
                        list_dmr_2.remove(dmr_2)
                        
            remain_dmr1 = []
            for dmr_1 in list_dmr_1 :
                remain_dmr1.append([min(dmr_1), max(dmr_1)])
            del list_dmr_1
            remain_dmr2 = []
            for dmr_2 in list_dmr_2 :
                remain_dmr2.append([min(dmr_2), max(dmr_2)])
            del list_dmr_2
            common_dmr_1 = []
            for common1 in common_dmr_from_1 :
                common_dmr_1.append([min(common1),max(common1)])
            common_dmr_2 = []
            for common2 in common_dmr_from_2 :
                common_dmr_2.append([min(common2),max(common2)])
                
            print "chr"+str(chr), len(remain_dmr1), len(remain_dmr2), len(common_dmr)
            for rd1 in remain_dmr1 :
                remain1.write(str(chr)+"\t"+str(rd1[0])+"\t"+str(rd1[1])+"\t"+str(val_dmr_1[remain_dmr1.index(rd1)])+"\n")
            #for rd2 in remain_dmr2 :
            #    remain2.write(str(chr)+"\t"+str(rd2[0])+"\t"+str(rd2[1])+"\t"+str(val_dmr_2[remain_dmr2.index(rd2)])+"\n")
            for cm1 in common_dmr_1 :
                common_file1.write(str(chr)+"\t"+str(cm1[0])+"\t"+str(cm1[1])+"\t"+str(common_val[common_dmr_1.index(cm1)][0])+"\n")
            #for cm2 in common_dmr_2 :
            #    common_file2.write(str(chr)+"\t"+str(cm2[0])+"\t"+str(cm2[1])+"\t"+str(common_val[common_dmr_2.index(cm2)][1])+"\n")
        common_file1.close()
        #common_file2.close()
        remain1.close()
        #remain2.close()
            
class unite :
    @staticmethod
    def pipeline(fastq1, fastq2) :
        print fastq1, fastq2
                  
#sys.exit()
#methylkit.methylkit_rpy(sys.argv[1], sys.argv[2], "CpG", 10)
if sys.argv[1] == "runCHH" :
    methylkit.methylkit_flow("CHH")
elif sys.argv[1] == "runCpG" :
    methylkit.methylkit_flow("CpG")
elif sys.argv[1] == "runCHG" :
    methylkit.methylkit_flow("CHG")
elif sys.argv[1] == "analysis" :
    analysis.make_bedGraph(sys.argv[2]) ## bedGraph of merged
elif sys.argv[1] == "circos" :
    analysis.make_circos_file(sys.argv[2]) ## bedGraph 
elif sys.argv[1] == "tile" :
    tile_class.tile_merge(name, 1000)
elif sys.argv[1] == "common" :
    analysis.commondmr(sys.argv[2], sys.argv[3]) ## bedgraph of DMR
elif sys.argv[1] == "test" :
    methylkit.prep_methylKit2(10, "CpG", "CML322_all3_bt202", "CML322_all3")
elif sys.argv[1] == "bismark" :
    bismark.single_alignment("test", "SRR641618.fastq")
#tile_class.tile_merge(name, 1000)
