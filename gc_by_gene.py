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
## grab all the protein coding region
for chr in range(1,11) :
    print chr
    coords = []
    for a in open("/work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23_protein_coding_exon_merged_transcripts.bed", 'r') :
        if a.split("\t")[0] == str(chr) and a.split("\t")[4].split("_")[1] == "T01\n" :
            coords.append(a.split("\t")[1]+"-"+(a.split("\t")[2]))
    large_gc_content = []
    for coord in coords :
        gc_content = []
        list_100 = []
        ## -10kb to 0 
        for nucl in whole_list[chr_list.index(str(chr))][int(coord.split("-")[0])-10001:int(coord.split("-")[0])] :
            if len(list_100) == 100 :
                gc_content.append(list_100.count("G")+list_100.count("C"))
                list_100 = []
                list_100.append(nucl)
            else :
                list_100.append(nucl)
        ## gene body
        body_length = len( whole_list[chr_list.index(str(chr))][int(coord.split("-")[0]):int(coord.split("-")[1])+1])
        for nucl in whole_list[chr_list.index(str(chr))][int(coord.split("-")[0]):int(coord.split("-")[1])+1] :
            ## normalize length
            if len(list_100) == int(round(float(body_length)/100,0)) :
                gc_content.append(list_100.count("G")+list_100.count("C"))
                list_100 = []
                list_100.append(nucl)
            else :
                list_100.append(nucl)
        ## 0 to +10kb
        for nucl in whole_list[chr_list.index(str(chr))][int(coord.split("-")[1]):int(coord.split("-")[1])+10001] :
            if len(list_100) == 100 :
                gc_content.append(list_100.count("G")+list_100.count("C"))
                list_100 = []
                list_100.append(nucl)
            else :
                list_100.append(nucl)
        if len(gc_content) == 300 :
            large_gc_content.append(gc_content)
    ### average out by number of gene coordinates and draw a graph
    total_list = []
    for x in range(0,300) :
        total = 0
        for gc in large_gc_content : 
            total += gc[x]
        total_list.append(float(total)/float(len(large_gc_content)))
    
    import numpy as np
    import matplotlib.pyplot as plt
    y = total_list
    x = range(0,300)
    colors = np.random.rand(300)
    plt.scatter(x, y, c=colors, alpha=0.5)
    plt.show()
