import os, sys
import numpy as np

### 
#genos = ['B73', 'Mo17', 'Oh43', 'Tx303', 'CML322']
geno = sys.argv[1]

chrom_length = {}
for a in open("/work/03076/gzynda/public/AGPv4.33/Zea_mays.AGPv4.dna.chrom.fa.fai", 'r') : 
    #chrom_length[a.split("\t")[0]] = int(round(int(a.split("\t")[1]),-2))
    chrom_length[a.split("\t")[0]] = np.zeros(((int(a.split("\t")[1])/100)+2,3))
#for chrom in chrom_length.keys():
#    print chrom, len(chrom_length[chrom])
#sys.exit()    
#os.chdir(geno)
#os.system('gzip -d '+geno+'_methratio.txt.gz')
with open(geno+"_methratio.txt") as f : 
    next(f)
    for line in f :  
        if line.split("\t")[3] == "CG" :
            chrom_length[line.split("\t")[0]][int(line.split("\t")[1])/100][0] += 1
        elif line.split("\t")[3] == "CHG" :
            chrom_length[line.split("\t")[0]][int(line.split("\t")[1])/100][1] += 1
        elif line.split("\t")[3] == "CHH" :
            chrom_length[line.split("\t")[0]][int(line.split("\t")[1])/100][2] += 1
 
outfile = open(geno+"_summary.bed", 'w') 
outfile.write("chr\tstart\tend\t#CG\t#CHG\t#CHH\n")
for chrom in [1,2,3,4,5,6,7,8,9,10,"Pt","Mt"] :
    chrom = str(chrom)
    for x in range(0,len(chrom_length[chrom])) :
        if chrom_length[chrom][x].tolist() != [0,0,0] :
            #print chrom, x, len(chrom_length[chrom]), chrom_length[chrom][x]
            outfile.write(chrom+"\t")
            outfile.write(str(x*100+1)+"\t")
            outfile.write(str((x+1)*100)+"\t")
            listval = chrom_length[chrom][x].tolist()
            outfile.write(str(int(listval[0]))+"\t")
            outfile.write(str(int(listval[1]))+"\t")
            outfile.write(str(int(listval[2]))+"\n")
            #print chrom_length[chrom][x]
            #outfile.write(str("\t".join(chrom_length[chrom][x])+"\n"))
outfile.close()

