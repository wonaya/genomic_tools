import os,sys
"""
annot = open("/work/02114/wonaya/genome/annotation/TAIR10_GFF3_genes.gff", 'r')
annotlines = annot.readlines()
outannot = open("TAIR10_polyphemus.txt", 'w')
for annotline in annotlines :
    if annotline.split("\t")[2] == "gene" :
        outannot.write("chr"+str(annotline.split("\t")[0].strip("Chr")))
        outannot.write("\t")
        outannot.write(str(annotline.split("\t")[3])+"\t")
        outannot.write(str(annotline.split("\t")[4])+"\t")
        outannot.write(str(annotline.split("\t")[8].split(";")[0].strip("ID="))+"\t")
        outannot.write(str(annotline.split("\t")[6])+"\t")
        outannot.write(str(annotline.split("\t")[8].split(";")[2].strip("Name=")))
sys.exit()
"""
a = open(sys.argv[1], 'r')
alines = a.readlines()
#chr_list = ['1','2','3','4','5','Mt','Pt']
chr_list = ['1']
for chr in chr_list :
    print chr
    os.chdir("pooled")
    outfile = open(sys.argv[1].split(".xls")[0]+"_chr"+str(chr)+".txt", 'w')
    for aline in alines[1:] :
        if aline[0] is not '#' and len(aline.split("\t")) > 3 and aline.split("\t")[0] == str(chr) :
            outfile.write(aline.split("\t")[0]+"\t")
            outfile.write(str(int(aline.split("\t")[1])+int(aline.split("\t")[4]))+"\t")
            outfile.write(aline.split("\t")[5]+"\t")
            outfile.write(str(10**(-1*float(aline.split("\t")[6]))))
            outfile.write("\n")
    outfile.close()
    os.chdir("..")
