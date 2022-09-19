import os,sys

def smooth(level,list_of_bams):
    chr_list = range(1,11)
    chr_list2 = []
    for chr in chr_list :
        chr_list2.append("Chr"+str(chr))
    chr_list = chr_list2
    bams = []
    for a in open(list_of_bams, 'r' ): 
        bams.append(a.strip("\n"))
    for bam in bams[1:] :
        os.system("rm -Rf tmp") 
        os.mkdir("tmp")
        os.system("rm -Rf "+bam.strip(".bam")+"_logFC."+str(level)+".smooth.wig")
        print bam.strip(".bam")+"_logFC.bedGraph"
        os.system("sort -k1,1 -k2,2g -o tmp/"+bam.strip(".bam")+"_logFC.wig "+bam.strip(".bam")+"_logFC.bedGraph")
        for chr in chr_list :
            print chr
            outfile = open("tmp/"+str(chr)+".out", 'w')
            for a in open("tmp/"+bam.strip(".bam")+"_logFC.wig" ,'r') :
                if a.split("\t")[0] == chr :
                    outfile.write(a)
            outfile.close()
            os.system("cut -f4 tmp/"+str(chr)+".out > tmp/"+str(chr)+".val")
            os.system("/work/02114/wonaya/software/hotspot-distr-v4/hotspot-deploy/bin/wavelets --level "+str(level)+" --to-stdout --boundary reflected --filter Haar tmp/"+str(chr)+".val > tmp/"+str(chr)+".smooth")
            os.system("paste tmp/"+str(chr)+".out tmp/"+str(chr)+".smooth | cut -f 1,2,3,5 >> "+bam.strip(".bam")+"_logFC."+str(level)+".smooth.wig")
        os.system("rm -Rf tmp") 
        
smooth(1, "list_bam.txt")
smooth(2, "list_bam.txt")
smooth(3, "list_bam.txt")
smooth(4, "list_bam.txt")
smooth(5, "list_bam.txt")
