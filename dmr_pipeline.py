import os,sys

if not os.path.isfile("/scratch/02114/wonaya/UMN_Bismark/B73_Class_array_dmr.bedGraph") :
    a = open("/scratch/02114/wonaya/UMN_Bismark/misc/9899_DMRs_across_51_genos.seg", 'r') 
    alines = a.readlines()
    class_interest = ['B73_Class', 'CML322_Class', 'MO17_Class', 'OH43_Class', 'TX303_Class']
    for classi in class_interest :
        print classi
        outfile = open(classi+"_array_dmr.bedGraph", 'w')
        for aline in alines[1:] :
            if aline.split("\t")[0] == classi :
                outfile.write(str(aline.split("\t")[1])+"\t")
                outfile.write(str(aline.split("\t")[2])+"\t")
                outfile.write(str(aline.split("\t")[3])+"\t")
                outfile.write(str(aline.split("\t")[5]))
        outfile.close()

b = open(sys.argv[1], 'r')
blines = b.readlines()
outfile = open(str(sys.argv[1]).split(".txt")[0]+".bedGraph", 'w')
for bline in blines[1:] :
    outfile.write(str(bline.split("\t")[0].split("chr")[1])+"\t")
    outfile.write(str(bline.split("\t")[1])+"\t")
    outfile.write(str(bline.split("\t")[2])+"\t")
    outfile.write(str(bline.split("\t")[4])+"\n")
outfile.close()