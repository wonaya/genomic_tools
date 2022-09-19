import os,sys
import numpy as np

logFC_list = []
for x in open("edgeR-05-vs-01.csv", 'r') :
    if x.split(",")[3] != '"PValue"' : 
        logFC_list.append(float(x.split(",")[1]))

a = np.array(logFC_list)
p95 = np.percentile(a, 95)
p05 = np.percentile(a, 5)

## grab into genes
enriched_01 = [] 
enriched_05 = [] 
for x in open("edgeR-05-vs-01.csv", 'r') :
    if x.split(",")[3] != '"PValue"' : 
        if float(x.split(",")[1]) >= p95 : 
            enriched_01.append(x.split(",")[0][1:-1])
        elif float(x.split(",")[1]) <= p05 : 
            enriched_05.append(x.split(",")[0][1:-1])

outfile1 = open("0-1mm_enriched_genes.txt", 'w')
for x in enriched_01 :
    outfile1.write(x)
    outfile1.write("\n")
outfile1.close()

outfile2 = open("0-5mm_enriched_genes.txt", 'w')
for x in enriched_05 :
    outfile2.write(x)
    outfile2.write("\n")
outfile2.close()
