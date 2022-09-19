import os,sys
### annotated

gene_list = []
for a in open(sys.argv[1], 'r') :
    gene_list.append(a.split("\t")[4].split("_")[0].strip("\n"))

mapman_list = []
for a in open("/work/02114/wonaya/genome/annotation/MapMan_B73_annot.txt", 'r') :
    mapman_list.append(a.split("\t")[1].strip("\n"))
    
mapman_dict = {}
for mapman in list(set(mapman_list)) :
    mapman_dict[mapman] = 0

for gene in list(set(gene_list)) :
    for a in open("/work/02114/wonaya/genome/annotation/MapMan_B73_annot.txt", 'r') : 
        if gene == a.split("\t")[0] :
            mapman_dict[a.split("\t")[1].strip("\n")] += 1

outfile = open(sys.argv[2], 'w') 
for mapman in mapman_dict.keys(): 
    if mapman_dict[mapman] > 0 : 
        outfile.write(mapman+"\t"+str(mapman_dict[mapman])+"\n")
outfile.close() 
