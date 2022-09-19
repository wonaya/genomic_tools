import os,sys

## smooth bedGraphs
bedGraph = sys.argv[1]
level = sys.argv[2]

os.mkdir("tmp")
os.system("rm -Rf "+sys.argv[1].split(".bedGraph")[0]+str(level)+".smooth.wig")
os.system("sort -k1,1 -k2,2g -o tmp/"+sys.argv[1].split(".bedGraph")[0]+".wig "+sys.argv[1])
for chr in range(1,11) :
    print chr
    outfile = open("tmp/"+str(chr)+".out", 'w')
    for a in open("tmp/"+sys.argv[1].split(".bedGraph")[0]+".wig" ,'r') :
        if a.split("\t")[0] == chr :
            outfile.write(a)
    outfile.close()
    os.system("cut -f4 tmp/"+str(chr)+".out > tmp/"+str(chr)+".val")
    os.system("/work/02114/wonaya/software/hotspot-distr-v4/hotspot-deploy/bin/wavelets --level "+str(level)+" --to-stdout --boundary reflected --filter Haar tmp/"+str(chr)+".val > tmp/"+str(chr)+".smooth")
    os.system("paste tmp/"+str(chr)+".out tmp/"+str(chr)+".smooth | cut -f 1,2,3,5 >> "+sys.argv[1].split(".bedGraph")[0]+"_logFC_"+str(level)+".smooth.wig")
os.system("rm -Rf tmp")        
    
