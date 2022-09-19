#!/usr/bin/env python

import os,sys
from optparse import OptionParser

class tile_class :
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
    
    def tile_context_sp_single_merge(samname, genomefile,context,window) :
        outfile = open(str(samname)+"_tile_"+str(context)+"_"+str(window)+"bp.txt", 'w')
        outfile.write("Start\tEnd\tChrom\tF CpG metC\tF CpG C\tR CpG metC\tR CpG C\n")
        for chr in setup.get_genome_size(genomefile)[1]:
            for line in open(str(samname)+"_tile_"+str(context)+"_"+str(window)+"bp_"+str(chr)+".txt", 'r') :
                outfile.write(line)
            os.remove(str(samname)+"_tile_"+str(context)+"_"+str(window)+"bp_"+str(chr)+".txt")
        outfile.close()
    
    def make_bedGraph(samname, genomefile,context,window, stranded) :
        outfile = open(str(samname)+"_tile_"+str(context)+"_"+str(window)+"bp.bedGraph", 'w')
        if stranded == "no" :
            for a in open(str(samname)+"_tile_"+str(context)+"_"+str(window)+"bp.txt", 'r') :
                if a.split("\t")[0] != "Start" :
                    outfile.write(a.split("\t")[2]+"\t"+a.split("\t")[0]+"\t"+a.split("\t")[1]+"\t"+str(((float(a.split("\t")[3])+float(a.split("\t")[5]))/((float(a.split("\t")[4])+float(a.split("\t")[6].strip("\n"))))))+"\n")
        outfile.close()
        
    def get_coverage(samname, context, window, samtools, bedtools) :
        os.system(str(samtools)+" view -bS "+str(samname)+" > "+str(samname).split(".sam")[0]+".bam")
        os.system(str(samtools)+" sort "+str(samname).split(".sam")[0]+".bam "+str(samname).split(".sam")[0]+".sorted")
        os.system(str(samtools)+" index "+str(samname).split(".sam")[0]+".sorted.bam")
        os.system(str(bedtools)+" -bams "+str(samname).split(".sam")[0]+".sorted.bam -bed "+str(samname)+"_tile_"+str(context)+"_"+str(window)+"bp.bedGraph > "+str(samname)+"_tile_"+str(context)+"_"+str(window)+"bp.cov")

