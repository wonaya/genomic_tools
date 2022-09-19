import os,sys

genome_list = []
whole_list = []
chr_list = []
chr_meant_list = ['10', '9', '8', '7', '6', '5', '4', '3', '2', '1']
print "reading in genome"
for genome in open("/work/02114/wonaya/genome/Zea_mays.AGPv3.23/Zea_mays.AGPv3.23.dna.genome.fa", 'r') :
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

for chr in range(10,11) :
    print chr
    if not os.path.isfile("CpG_OT_unique_2_sortedn_chr"+str(chr)+".txt") :
        outfile = open("CpG_OT_unique_2_sortedn_chr"+str(chr)+".txt", 'w') 
        for a in open("CpG_OT_unique_2_sortedn.txt", 'r') :
            if len(a.split("\t")) > 1 and a.split("\t")[2] == str(chr) :
                outfile.write(a)
        outfile.close()
    if not os.path.isfile("CpG_OB_unique_2_sortedn_chr"+str(chr)+".txt") :
        outfile = open("CpG_OB_unique_2_sortedn_chr"+str(chr)+".txt", 'w') 
        for a in open("CpG_OB_unique_2_sortedn.txt", 'r') :
            if len(a.split("\t")) > 1 and a.split("\t")[2] == str(chr) :
                outfile.write(a)
        outfile.close()
    
    pos = []
    a = open("CpG_OT_unique_2_sortedn_chr"+str(chr)+".txt", 'r') 
    alines = a.readlines()
    for aline in alines[1:] :
        pos.append(int(aline.split("\t")[3]))
    
    a = open("CpG_OB_unique_2_sortedn_chr"+str(chr)+".txt", 'r') 
    alines = a.readlines()
    for aline in alines[1:] :
        pos.append(int(aline.split("\t")[3]))
    print len(set(pos))
    
    ## grab all the protein coding region
    coords = []
    for a in open("/work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23_protein_coding_exon_merged_transcripts.bed", 'r') :
        if a.split("\t")[0] == str(chr) and a.split("\t")[4].split("_")[1] == "T01\n" :
            coords.append(a.split("\t")[1]+"-"+(a.split("\t")[2]))
    large_metrate = []
    for coord in coords[:100] :
        metrate = []
        list_100 = []
        pos_coord = []
        for a in set(pos) :
            if a >= int(coord.split("-")[0])-10001 and a <= int(coord.split("-")[1])+10000 : 
                pos_coord.append(a)
        pos_met = [0]*len(pos_coord)
        pos_unmet = [0]*len(pos_coord)
        for a in open("CpG_OT_unique_2_sortedn_chr"+str(chr)+".txt", 'r') :
            if int(a.split("\t")[3]) in pos_coord :
                pos_coord.index(int(a.split("\t")[3]))
                if a.split("\t")[1] == "+" :
                    pos_met[pos_coord.index(int(a.split("\t")[3]))] += 1
                elif a.split("\t")[1] == "-" :
                    pos_unmet[pos_coord.index(int(a.split("\t")[3]))] += 1
        for a in open("CpG_OB_unique_2_sortedn_chr"+str(chr)+".txt", 'r') :
            if int(a.split("\t")[3]) in pos_coord :
                pos_coord.index(int(a.split("\t")[3]))
                if a.split("\t")[1] == "+" :
                    pos_met[pos_coord.index(int(a.split("\t")[3]))] += 1
                elif a.split("\t")[1] == "-" :
                    pos_unmet[pos_coord.index(int(a.split("\t")[3]))] += 1
        
        ## -10kb to 0 
        for x in range(0,100) :
            total_met = 0
            total_unmet = 0
            for loc in range(int(coord.split("-")[0])-(10001-(x*100)),int(coord.split("-")[0])-(10001-((x+1)*100))) :
                if loc in pos_coord :
                    total_met += pos_met[pos_coord.index(loc)]
                    total_unmet += pos_unmet[pos_coord.index(loc)]
            if total_met+total_unmet == 0 :
                metrate.append(0)
            else: 
                metrate.append((float(total_met)/float(total_met+total_unmet)*100))
        
            
        """
        ## gene body
        body_length = coord.split("-")[1]-coord.split("-")[0]
        for nucl in whole_list[chr_list.index(str(chr))][int(coord.split("-")[0]):int(coord.split("-")[1])+1] :
            ## normalize length
            if len(list_100) == int(round(float(body_length)/100,0)) :
                gc_content.append(list_100.count("G")+list_100.count("C"))
                list_100 = []
                list_100.append(nucl)
            else :
                list_100.append(nucl)
        
        """
        ## 0 to +10kb
        for x in range(0,100) :
            total_met = 0
            total_unmet = 0
            for loc in range(int(coord.split("-")[1])+(x*100),int(coord.split("-")[1])+((x+1)*100)) :
                if loc in pos_coord :
                    total_met += pos_met[pos_coord.index(loc)]
                    total_unmet += pos_unmet[pos_coord.index(loc)]
            if total_met+total_unmet == 0 :
                metrate.append(0)
            else: 
                metrate.append((float(total_met)/float(total_met+total_unmet)*100))
        if len(metrate) == 200 :
            large_metrate.append(metrate)
                
    ### average out by number of gene coordinates and draw a graph
    total_list = []
    for x in range(0,200) :
        total = 0
        for metrate in large_metrate : 
            total += metrate[x]
        total_list.append(float(total)/float(len(large_metrate)))
    
    import numpy as np
    import matplotlib.pyplot as plt
    y = total_list
    x = range(0,200)
    colors = np.random.rand(200)
    plt.scatter(x, y, c=colors, alpha=0.5)
    plt.show()
