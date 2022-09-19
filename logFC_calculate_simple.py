import os,sys
import math

def roundup(x, window_size):
    return int(math.ceil(x / int(window_size))) * int(window_size)

def different_samples() :
    a = open("g1.test.count", 'r') 
    alines = a.readlines()
    b = open("brdu.merged.test.count", 'r') 
    blines = b.readlines()
    c = open("g1_brdu.test.count",'w')
    for aline in alines :
        c.write("\t".join(aline.strip("\n").split("\t")))
        c.write("\t")
        c.write("\t".join(blines[alines.index(aline)].split("\t")[3:]))
#different_samples()

def get_genome_size() :
    chr_list = []
    count_list = []
    count = 0
    for genome in open("/work/02114/wonaya/genome/Zea_mays.AGPv3.18/Zea_mays.AGPv3.18.dna.fa", 'r') :
        if genome[0] == ">" :
            if count > 0 :
                count_list.append(count)
            chr_list.append(genome[1:].split(" dna")[0])
            count = 0
        else :
            count += len(genome)
    count_list.append(count)
    chr_dict = {}
    for chr in chr_list :
        chr_dict[chr] = count_list[chr_list.index(chr)]
    return chr_dict, chr_list
       
def make_window_bed(window_size, genome) :
    outfile = open("window_"+str(window_size)+".bed", 'w')
    for chr in range(1,11): ## for the case of maize
        for x in range(0, (int(roundup(get_genome_size()[0][str(chr)], window_size))/window_size)+1) :
            #outfile.write("chr"+str(chr)+"\t")
            outfile.write(str(chr)+"\t")
            outfile.write(str((x*window_size*1)+1)+"\t")
            outfile.write(str((x+1)*window_size)+"\n")
    outfile.close()

def run_logFC(list_of_bams, window_size) :
    ### always put input at the top row
    bams = []
    """
    for a in open(list_of_bams, 'r' ): 
        bams.append(a.strip("\n"))
        if not os.path.isfile(a.strip("\n")+".bai") :
            os.system("/opt/apps/samtools/0.1.19/samtools index "+a.strip("\n")) 
    print " ".join(bams)
    os.system("/opt/apps/bedtools/2.19.0/bin/multiBamCov -bams "+" ".join(bams)+" -bed window_"+str(window_size)+".bed > "+list_of_bams+".test.count")
    """
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
    
    sum_list_norm = [0]*len(bams)
    for a in open(list_of_bams+".test_norm.count", 'r') :
        count = 0
        for val in a.split("\t")[3:] :
            sum_list_norm[count] += float(val.strip("\n"))
            count += 1
    print sum_list_norm
    
    #Log2 (ip*k / cntr+1)  k= read counts of cntrl/read counts of IP
    bams = []
    for a in open(list_of_bams, 'r' ): 
        bams.append(a.strip("\n"))
    for bam in bams[1:] :
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
    
def replicon(level):
    bams = []
    for a in open(list_of_bams, 'r' ): 
        bams.append(a.strip("\n"))
    for bam in bams[1:] :
        print bam
        outfile = open(bam.strip(".bam")+"_logFC."+str(level)+".replicon.bedGraph", 'w')
        for chr in range(1,11) :
            coord_list = []
            value_list = []
            start = []
            end = []
            first_loc_list = []
            second_loc_list = []
            val_list = []
            for a in open(bam.strip(".bam")+"_logFC."+str(level)+".smooth.wig", 'r') :
                if int(a.split("\t")[0]) == chr :
                    first_loc_list.append(int(a.split("\t")[1]))
                    second_loc_list.append(int(a.split("\t")[2]))
                    val_list.append(float(a.split("\t")[3].strip("\n")))
            if len(first_loc_list) == len(second_loc_list) == len(val_list) :
                for x in range(0, len(first_loc_list)) :
                    if x == 0 :
                        if float(val_list[x]) >= 0 :
                            start.append(1)
                    else :
                        if float(val_list[x]) >= 0 :
                            if len(start) == len(end) :
                                start.append(first_loc_list[x]-5000)
                        elif float(val_list[x]) < 0 :
                            if len(start) > len(end) :
                                end.append(first_loc_list[x]-5000)
            if len(start) > len(end) :
                end.append(second_loc_list[-1])
                
            print chr, len(start), len(end) ## should be the same?
            for no in range(0, len(start)) :
                if start[no] >= end[no] : 
                    print no, start[no], end[no]
                    sys.exit()
                outfile.write(str(chr)+"\t"+str(int(start[no]))+"\t"+str(int(end[no]))+"\t1\n")
        outfile.close()
        
def spreading():
    # spreading replicons    
    for chr in range(1,11) :
        outfile = open("ES-MS-spreading.6.replicon.bedGraph", 'w')
        for a in open("ES-merge_sorted_logFC.6.replicon.bedGraph", 'r') :
            for b in open("MS-merge_sorted_logFC.6.replicon.bedGraph", 'r') :
                if int(b.split("\t")[1]) <= int(a.split("\t")[1]) and int(b.split("\t")[2]) >= int(a.split("\t")[2]) :
                    outfile.write(b)
    outfile.close()
    
def overlap() : 
    # gene association 
    ## get GFF
    outfile = open("EdU-overlay.bedGraph",'w')
    chr_dict = get_genome_size()[0]
    for chr in range(1,11):
        print chr, chr_dict[str(chr)]
        large_list = []
        for coord in range(0, chr_dict[str(chr)]) :
             large_list.append([])
        for a in open("ES-merge_sorted_logFC.1.replicon.bedGraph", 'r') :
            if a.split("\t")[0] == str(chr) :
                for coord in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
                    coord.append('E')
        for a in open("MS-merge_sorted_logFC.1.replicon.bedGraph", 'r') :
            if a.split("\t")[0] == str(chr) :
                for coord in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
                    coord.append('M')
        for a in open("LS-merge_sorted_logFC.1.replicon.bedGraph", 'r') :
            if a.split("\t")[0] == str(chr) :
                for coord in large_list[int(a.split("\t")[1]):int(a.split("\t")[2])+1] :
                    coord.append('L')
        
        test = []
        score_dict = {'E':-3,'M':-2,'L':-1,'EM':0,'EL':1,'ML':2,'EML':3, '':4}
        for coord in range(1, chr_dict[str(chr)]) :
            if len(test) == 0 :
                id = "".join(large_list[coord])
                test.append(coord)
            else :
                if "".join(large_list[coord]) == "".join(large_list[coord-1]) :
                    test.append(coord)
                else :
                    outfile.write(str(chr)+"\t")
                    outfile.write(str(test[0])+"\t")
                    outfile.write(str(test[-1])+"\t")
                    outfile.write(str(score_dict[id])+"\n")
                    test = []
                    id = "".join(large_list[coord])
                    test.append(coord)
        outfile.close()
        sys.exit()
            
#make_window_bed(1000, "/work/02114/wonaya/genome/Zea_mays.AGPv3.18/Zea_mays.AGPv3.18.dna.fa")
run_logFC(sys.argv[1], 1000)        
#smooth(1,sys.argv[1])
#smooth(2,sys.argv[1])
#smooth(3,sys.argv[1])
#smooth(4,sys.argv[1])
#smooth(5,sys.argv[1])
#replicon(level)
#spreading()
#overlap()