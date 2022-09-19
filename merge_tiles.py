import os,sys

genos = ['B73', 'Mo17', 'Oh43', 'CML322', 'Tx303']
contexts = ['CpG', 'CHG', 'CHH']

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

outfile = open("temp_merge.txt", 'w') 
for chr in range(1,11) :
    print chr, int(chr_dict[str(chr)])/100+1
    large_list = []
    for x in range(0,int(chr_dict[str(chr)])/100+2):
        large_list.append(['NA']*15)
    index = 0
    for geno in genos :
        for context in contexts :
            print geno, context, index
            for a in open(str(geno)+"_all3_bt202_sortedn_tile_"+str(context)+"_100bp_merged_redone.txt", 'r') :
                if a.split("\t")[2] == str(chr) :
                    metc = int(a.split("\t")[3])+int(a.split("\t")[5])
                    totc = int(a.split("\t")[4])+int(a.split("\t")[6].strip("\n"))
                    metperc = float(metc)/float(totc)*100
                    large_list[int(a.split("\t")[0])/100][index] = metperc
            index += 1
    
    count = 0
    for large in large_list :
        if large != ['NA']*15 :
            outfile.write(str(chr))
            outfile.write("\t")
            outfile.write(str(count*100))
            outfile.write("\t")
            outfile.write(str(count*100+99))
            outfile.write("\t")
            outfile.write("\t".join(map(str,large)))
            outfile.write("\n")
        count += 1
outfile.close()