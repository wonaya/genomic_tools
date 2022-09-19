### merge all rmdup bedgraph 

import os,sys

def get_genome_size() :
    chr_list = []
    count_list = []
    count = 0
    for genome in open("/work/02114/wonaya/genome/Zea_mays.AGPv3.18.dna.fa", 'r') :
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
"""
genos = ['B73','Mo17','Oh43','CML322','Tx303']
contexts = ['CHG','CHH','CpG']
outfile = open("rmdup_merge_summary.bed", 'w')
outfile.write("chr\tstart\tend\tB73 CHG\tB73 CHH\tB73 CG\tMo17 CHG\tMo17 CHH\tMo17 CG\tOh43 CHG\tOh43 CHH\tOh43 CG\tCML322 CHG\tCML322 CHH\tCML322 CG\tTx303 CHG\tTx303 CHH\tTx303 CG\n")
for chr in range(1,11) :
    print "chr", chr, get_genome_size()[0][str(chr)]
    large_list = []
    for x in range(0,get_genome_size()[0][str(chr)]/100+1) :
        large_list.append(['NA']*15)
    index = 0
    for geno in genos :
        for context in contexts :
            print index, geno, context
            if os.path.isfile(geno+"_all3/"+geno+"_all3_bt202_sorted_picard_rmdup_tile_"+context+"_100bp_merged_weighted.txt") :
                a = open(geno+"_all3/"+geno+"_all3_bt202_sorted_picard_rmdup_tile_"+context+"_100bp_merged_weighted.txt", 'r')
                alines = a.readlines()
                for aline in alines[1:] :
                    if str(aline.split("\t")[2]) == str(chr) : 
                        large_list[int(aline.split("\t")[0])/100][index] = (float(aline.split("\t")[3])+float(aline.split("\t")[5]))/(float(aline.split("\t")[4])+float(aline.split("\t")[6].strip("\n")))*100
            else : 
                print "somethings wrong",sys.exit()
            index += 1
    count = 0
    for x in large_list :
        if x != ['NA']*15 : 
            outfile.write("chr"+str(chr)+"\t")
            outfile.write(str(count*100)+"\t")
            outfile.write(str(count*100+99)+"\t")
            outfile.write("\t".join(map(str,x)))
            outfile.write("\n")
        count += 1

outfile.close()
sys.exit()
"""
import numpy 

## get correlation per chr
for chr in range(1,11) :
    coord_list = []
    val_list = []
    a = open("/corral-tacc/tacc/iplant/vaughn/springer_vaughn/eichten/jawon/05-29-14-100bp_tiles_redo/window100_5geno_data_redo.bed",'r') 
    alines = a.readlines()
    for aline in alines[1:] :
        if aline.split("\t")[0] == "chr"+str(chr) :
            coord_list.append(aline.split("\t")[1]) 
            val_list.append(aline.split("\t")[3])
    print len(alines)
    b = open("rmdup_merge_summary.bed", 'r') 
    blines = b.readlines()
    coord_list2 = []
    val_list2 = []
    for bline in blines[1:] :
        if bline.split("\t")[0] == "chr"+str(chr) :
            coord_list2.append(bline.split("\t")[1]) 
            val_list2.append(bline.split("\t")[3])
    corr_val1 = []
    corr_val2 = []
    print len(coord_list), len(coord_list2)
    for coord in coord_list :
        if coord in coord_list2 and val_list[coord_list.index(coord)] != "NA" and val_list2[coord_list2.index(coord)] != "NA" :
            #if float(val_list[coord_list.index(coord)]) > 0 and float(val_list2[coord_list2.index(coord)]) > 0 : 
            corr_val1.append(float(val_list[coord_list.index(coord)]))
            corr_val2.append(float(val_list2[coord_list2.index(coord)]))
    print chr, numpy.corrcoef(corr_val1,corr_val2)[0], len(corr_val1), len(corr_val2)
    sys.exit()
                
            
    
    