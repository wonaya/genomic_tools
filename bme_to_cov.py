import os,sys
import multiprocessing

def get_genome_size() :
    chr_list = []
    count_list = []
    count = 0
    for genome in open("/work/02114/wonaya/genome/Zea_mays.AGPv3.18.dna.fa", 'r') :
        if genome[0] == ">" :
            if count > 0 :
                count_list.append(count)
            chr_list.append(genome[1:].split(" ")[0].strip("\n"))
            count = 0
        else :
            count += len(genome)
    count_list.append(count)
    chr_dict = {}
    for chr in chr_list :
        chr_dict[chr] = count_list[chr_list.index(chr)]
    return chr_dict, chr_list

def extract_name(samname) :
    dir = "."
    if len(samname.split("/")) > 1 :
        dir = samname.split("/")[0]
        filename = samname.split("/")[1]
    else :
        dir = "."
        samname = samname
    name = samname.split(".sam")[0]
    return name, dir

def bme_to_cov(chr, samname, context) :
    
    large_list_1 = []
    for x in range(0,get_genome_size()[0][chr]+1) :
        large_list_1.append([0]*4)
    for line1 in open(str(context)+"_OT_"+extract_name(samname)[0]+"_"+str(chr)+".txt", 'r') :
        if line1.split("\t")[1] == "+" and int(line1.split("\t")[3]) <= len(large_list_1) :
            large_list_1[int(line1.split("\t")[3])][0] += 1
            large_list_1[int(line1.split("\t")[3])][1] += 1
        elif line1.split("\t")[1] == "-" and int(line1.split("\t")[3]) <= len(large_list_1):
            large_list_1[int(line1.split("\t")[3])][1] += 1
    for line1 in open(str(context)+"_CTOT_"+extract_name(samname)[0]+"_"+str(chr)+".txt", 'r') :
        if line1.split("\t")[1] == "+" and int(line1.split("\t")[3]) <= len(large_list_1):
            large_list_1[int(line1.split("\t")[3])][0] += 1
            large_list_1[int(line1.split("\t")[3])][1] += 1
        elif line1.split("\t")[1] == "-" and int(line1.split("\t")[3]) <= len(large_list_1):
            large_list_1[int(line1.split("\t")[3])][1] += 1
    for line1 in open(str(context)+"_OB_"+extract_name(samname)[0]+"_"+str(chr)+".txt", 'r') :
        if line1.split("\t")[1] == "+" and int(line1.split("\t")[3]) <= len(large_list_1):
            large_list_1[int(line1.split("\t")[3])][2] += 1
            large_list_1[int(line1.split("\t")[3])][3] += 1
        elif line1.split("\t")[1] == "-" and int(line1.split("\t")[3]) <= len(large_list_1):
            large_list_1[int(line1.split("\t")[3])][3] += 1
    for line1 in open(str(context)+"_CTOB_"+extract_name(samname)[0]+"_"+str(chr)+".txt", 'r') :
        if line1.split("\t")[1] == "+" and int(line1.split("\t")[3]) <= len(large_list_1):
            large_list_1[int(line1.split("\t")[3])][2] += 1
            large_list_1[int(line1.split("\t")[3])][3] += 1
        elif line1.split("\t")[1] == "-" and int(line1.split("\t")[3]) <= len(large_list_1) :
            large_list_1[int(line1.split("\t")[3])][3] += 1
    outfile = open(str(context)+"_"+extract_name(samname)[0]+"_"+str(chr)+".cov", 'w')
    coord = 0
    for large_list in large_list_1 :
        if large_list[1] > 0 or large_list[3] > 0 :
            outfile.write("chr"+str(chr)+"\t")
            outfile.write(str(coord)+"\t"+str(coord)+"\t")
            metperc = (float(large_list[0]+large_list[2])/float(large_list[1]+large_list[3]))*100
            outfile.write(str(metperc)+"\t")
            outfile.write(str(int(large_list[0]+large_list[2]))+"\t")
            outfile.write(str(int(large_list[1]+large_list[3])-int(large_list[0]+large_list[2]))+"\n")
        coord += 1 
    outfile.close()

context = sys.argv[1]
geno = sys.argv[2]
jobs = []
for chr in range(1,11) :
    print chr
    if geno == "B73" :
        s = multiprocessing.Process(target=bme_to_cov, args=(str(chr), str(geno)+"_all3_R1_bt202_sortedn.sam", str(context), ))
    else :
        s = multiprocessing.Process(target=bme_to_cov, args=(str(chr), str(geno)+"_all3_bt202_sortedn.sam", str(context), ))
    jobs.append(s)
    s.start()
[x.join() for x in jobs]

sys.exit()
for chr in range(1,11) :
    print chr
    if geno == "B73" :
        os.system("cat "+str(context)+"_"+str(geno)+"_all3_R1_bt202_sortedn_"+str(chr)+".cov >> "+str(context)+"_"+str(geno)+"_all3_bt202_sortedn.cov")
        #os.system("rm "+str(context)+"_"+str(geno)+"_all3_R1_bt202_sortedn_"+str(chr)+".cov")
    else :
        os.system("cat "+str(context)+"_"+str(geno)+"_all3_bt202_sortedn_"+str(chr)+".cov >> "+str(context)+"_"+str(geno)+"_all3_bt202_sortedn.cov")
        #os.system("rm "+str(context)+"_"+str(geno)+"_all3_bt202_sortedn_"+str(chr)+".cov")
    
