import os,sys

for chr in range(1,11) :
    print chr
    afile =  open("B73_Oh43_CpG_dmr_50.txt", 'r') 
    alines = afile.readlines()
    a_coord_list = []
    for a in alines[1:] :
        if a.split("\t")[0] == str(chr) :
            a_coord_list.append(range(int(a.split("\t")[1]), int(a.split("\t")[2])+1))
    
    bfile = open("B73_Oh43_CHG_dmr_50.txt", 'r')
    blines = bfile.readlines()
    b_coord_list = []
    for b in blines[1:] :
        if b.split("\t")[0] == str(chr) :
            b_coord_list.append(range(int(b.split("\t")[1]), int(b.split("\t")[2])+1))
    overlap_coord_list = []      
    for a_coord in a_coord_list :
        for b_coord in b_coord_list :
            if len(set(a_coord)&set(b_coord)) > 0 :  
                if range(min(min(a_coord), max(a_coord), min(b_coord), max(b_coord)), max(min(a_coord), max(a_coord), min(b_coord), max(b_coord))+1) not in overlap_coord_list :
                    overlap_coord_list.append(range(min(min(a_coord), max(a_coord), min(b_coord), max(b_coord)), max(min(a_coord), max(a_coord), min(b_coord), max(b_coord))+1))
    legit_a_coord_list  =[]
    for a_coord in a_coord_list :
        count = 0
        for overlap_coord in overlap_coord_list :
            if len(set(a_coord)&set(overlap_coord)) > 0 :
                count += 1
        if count == 0 :
            legit_a_coord_list.append(a_coord)
    
    legit_b_coord_list = []
    for b_coord in b_coord_list :
        count = 0
        for overlap_coord in overlap_coord_list :
            if len(set(b_coord)&set(overlap_coord)) > 0 :
                count += 1
        if count == 0 :
            legit_b_coord_list.append(b_coord)
    print len(legit_a_coord_list), len(legit_b_coord_list),len(overlap_coord_list)
    outfile = open("B73_Oh43_CpG_CHG_dmr_50_merged_temp_chr"+str(chr)+".txt", 'w') 
    for a_coord in legit_a_coord_list :
        outfile.write("chr"+str(chr)+"\t"+str(min(a_coord))+"\t"+str(max(a_coord))+"\t1\n")
    for b_coord in legit_b_coord_list :
        outfile.write("chr"+str(chr)+"\t"+str(min(b_coord))+"\t"+str(max(b_coord))+"\t2\n")
    for overlap_coord in overlap_coord_list :
        outfile.write("chr"+str(chr)+"\t"+str(min(overlap_coord))+"\t"+str(max(overlap_coord))+"\t3\n")
    outfile.close()
    os.system("sort -nk2 B73_Oh43_CpG_CHG_dmr_50_merged_temp_chr"+str(chr)+".txt > B73_Oh43_CpG_CHG_dmr_50_merged_chr"+str(chr)+".txt")
    os.system("rm -Rf B73_Oh43_CpG_CHG_dmr_50_merged_temp_chr"+str(chr)+".txt")
for chr in range(1,11) :
    os.system("cat B73_Oh43_CpG_CHG_dmr_50_merged_chr"+str(chr)+".txt >> B73_Oh43_CpG_CHG_dmr_50_merged.bedGraph")
    os.system("rm -Rf B73_Oh43_CpG_CHG_dmr_50_merged_chr"+str(chr)+".txt")