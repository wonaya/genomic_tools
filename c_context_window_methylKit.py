### python.py CpG/CHH/CHG B73_ACAGTG/Mo17_GCCAAT chr1...10

import os,sys

## save DMR
dmr_loc = {}
dmr_loc_start_val = []
methyl_type = str(sys.argv[1])
strain = str(sys.argv[2])
for c in open(methyl_type+"_context_chr10.txt", 'r') :
    if c.split("\t")[1].strip('"') == "chr10" :
        dmr_loc[int(c.split("\t")[2])] = int(c.split("\t")[3])
        dmr_loc_start_val.append(int(c.split("\t")[2]))     
dmr_loc_start_val.sort()   

## calculate weighted methylation
a = open(methyl_type+"_context_"+strain+"_merged.sorted_chr10.out", 'r') 
alines = a.readlines()
outfile = open(methyl_type+"_context_"+strain+"_merged.sorted_chr10_calc_wMet.out", 'w')
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
    outfile.write("chr10\t")
    outfile.write(str(dmr)+"\t")
    outfile.write(str(dmr_loc[dmr])+"\t")
    outfile.write(str(total_c)+"\t")
    outfile.write(str(total_read)+"\t")
    outfile.write(str(val_met)+"\t")
    outfile.write(str(weighted_top)+"\t")
    outfile.write(str(frag)+"\t"+str(mean)+"\t"+str(weighted)+"\n")
outfile.close()

