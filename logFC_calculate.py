import os,sys
import math

def roundup(x, window_size):
    return int(math.ceil(x / int(window_size))) * int(window_size)
       
def make_window_bed(window_size, genome) :
    outfile = open("window_"+str(window_size)+".bed", 'w')
    for chr in range(1,11): ## for the case of maize
        for x in range(0, (int(roundup(get_genome_size()[0][str(chr)], window_size))/window_size)+1) :
            outfile.write(str(chr)+"\t")
            outfile.write(str((x*window_size*1)+1)+"\t")
            outfile.write(str((x+1)*window_size)+"\n")
    outfile.close()

def run_logFC(list_of_bams, window_size) :
    bams = []
    for a in open(list_of_bams, 'r' ): 
        bams.append(a.strip("\n"))
        if not os.path.isfile(a.strip("\n")+".bai") :
            os.system("/opt/apps/samtools/0.1.19/samtools index "+a.strip("\n")) 
    print " ".join(bams)
    os.system("/opt/apps/bedtools/2.19.0/bin/multiBamCov -bams "+" ".join(bams)+" -bed window_"+str(window_size)+".bed > "+list_of_bams+".test.count")
    
    bams = []
    for a in open(list_of_bams, 'r' ): 
        bams.append(a.strip("\n"))
    print len(bams), bams
    sum_list = [0]*len(bams)
    for a in open(list_of_bams+".test.count", 'r') :
        count = 0
        for val in a.split("\t")[3:] :
            sum_list[count] += int(val.strip("\n"))
            count += 1
    print sum_list
    norm_factor = [] 
    for sum in sum_list :
        norm_factor.append(float(sum)/float(min(sum_list)))
    print norm_factor
    outfile = open(list_of_bams+".test_norm.count", 'w')
    for a in open(list_of_bams+".test.count", 'r') :  
        outfile.write("\t".join(a.split("\t")[:3]))
        count = 0
        for val in a.split("\t")[3:] : 
            outfile.write("\t")
            outfile.write(str(float(val)/norm_factor[count]))
            count += 1
        outfile.write("\n")
    outfile.close()

    #Log2 (ip*k / cntr+1)  k= read counts of cntrl/read counts of IP
    bams = []
    for a in open(list_of_bams, 'r' ): 
        bams.append(a.strip("\n"))
    for bam in bams :
        print bam, bams.index(bam)
        outfile = open(bam.strip(".bam")+"_logFC.bedGraph", 'w')
        for a in open(list_of_bams+".test_norm.count", 'r') :
            outfile.write("\t".join(a.split("\t")[:3]))
            outfile.write("\t")
            if float(a.split("\t")[3+bams.index(bam)]) == 0 : 
                outfile.write("0")
            else :
                ### change to whichever input used
                outfile.write(str(math.log(float(a.split("\t")[3+bams.index(bam)])/(float(a.split("\t")[3])+1),2)))
            outfile.write("\n")
        outfile.close()

def smooth(level,list_of_bams):
    chr_list = range(1,11)
    chr_list2 = []
    for chr in chr_list :
        chr_list2.append(str(chr))
    chr_list = chr_list2
    print chr_list
    bams = []
    for a in open(list_of_bams, 'r' ): 
        bams.append(a.strip("\n"))
    for bam in bams[1:] :
        print bam
        os.mkdir("tmp")
        os.system("rm -Rf "+bam.strip(".bam")+"_logFC."+str(level)+".smooth.wig")
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
            
make_window_bed(1000, "/work/02114/wonaya/genome/Zea_mays.AGPv3.18.dna.fa")
run_logFC(sys.argv[1], 1000)        
smooth(1,sys.argv[1])
smooth(2,sys.argv[1])
smooth(3,sys.argv[1])
smooth(4,sys.argv[1])
smooth(5,sys.argv[1])
