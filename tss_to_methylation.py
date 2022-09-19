import os,sys

### maize tss 

### extract start codon
outfile = open("b73_tss_10kb_metrate.txt", 'w')
for chr in range(10,11) :
    arb_list = []
    x = -10000 
    while x <= 10000 :
        arb_list.append(x)
        x += 100
    print arb_list
    sys.exit()
    outfile.write("Gene\t")
    outfile.write("\t".join(map(str,arb_list)))
    outfile.write("\n")
    gene_list = []
    tss_list = []
    for a in open("/work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23.gtf", 'r') :
        if a.split("\t")[0] == str(chr) and a.split("\t")[2] == "start_codon" : 
            tss_list.append(int(a.split("\t")[4]))
            gene_list.append(a.split("\t")[8].split(";")[0].split('"')[1])
    print chr, len(gene_list), len(tss_list)
    
    ### save methylation 
    
    tile_list = []
    metval_list = []
    for a in open("/scratch/02114/wonaya/UMN_Bismark/10-01-13_5genos/B73_all3/B73_all3_bt202_tile_CpG_100bp_merged_redone.txt", 'r'):
        if a.split("\t")[2] == str(chr) :
            tile_list.append(int(a.split("\t")[0]))
            metval_list.append(((float(a.split("\t")[3])+float(a.split("\t")[5]))/(float(a.split("\t")[4])+float(a.split("\t")[6].strip("\n"))))*100)
    ### get +/- 10kb of TSS
    index = 0
    metval_large = []
    
    for gene in list(set(gene_list)) :
        print index
        metval_small = [0]*len(arb_list)
        for tile in tile_list :
            if tile >= round(tss_list[index],-2)-10000 and tile <= round(tss_list[index],-2)+10000 :
                #print (round(tss_list[index],-2)-tile)*-1, arb_list.index((round(tss_list[index],-2)-tile)*-1)
                metval_small[arb_list.index((round(tss_list[index],-2)-tile)*-1)] = metval_list[tile_list.index(tile)]
                #metval_small.append(metval_list[tile_list.index(tile)])
        metval_small.insert(0, gene)
        metval_large.append(metval_small)
        index += 1
    
    for metval in metval_large :
        outfile.write("\t".join(map(str,metval)))
        outfile.write("\n")
outfile.close()