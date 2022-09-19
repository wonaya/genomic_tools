import os,sys
import datetime

for chr in range(10,11) :
    if not os.path.isfile("S1_methratio_CG_chr"+str(chr)+".txt") :
        outfile = open("S1_methratio_CG_chr"+str(chr)+".txt", 'w') 
        for a in open("S1_methratio.txt", 'r') :
            if a.split("\t")[0] == str(chr) and a.split("\t")[3] == "CG" :
                outfile.write(a)
        outfile.close()
    
    if not os.path.isfile("S2_methratio_CG_chr"+str(chr)+".txt") :
        outfile = open("S2_methratio_CG_chr"+str(chr)+".txt", 'w') 
        for a in open("S2_methratio.txt", 'r') :
            if a.split("\t")[0] == str(chr) and a.split("\t")[3] == "CG" :
                outfile.write(a)
        outfile.close()
    
    if not os.path.isfile("S3_methratio_CG_chr"+str(chr)+".txt") :
        outfile = open("S3_methratio_CG_chr"+str(chr)+".txt", 'w') 
        for a in open("S3_methratio.txt", 'r') :
            if a.split("\t")[0] == str(chr) and a.split("\t")[3] == "CG" :
                outfile.write(a)
        outfile.close()
    print "start"
    pos_list = []
    metc_list = []
    totc_list = []
    for a in open("S1_methratio_CG_chr"+str(chr)+".txt", 'r') :
        pos_list.append(int(a.split("\t")[1]))
    for a in open("S2_methratio_CG_chr"+str(chr)+".txt", 'r') :
        pos_list.append(int(a.split("\t")[1]))
    for a in open("S3_methratio_CG_chr"+str(chr)+".txt", 'r') :
        pos_list.append(int(a.split("\t")[1]))
    print "done"
    set_pos_list = []
    for pos in set(pos_list) :
        set_pos_list.append(pos)
    set_pos_list.sort()
    print len(set(pos_list)), datetime.datetime.now()
    
    metc_list = [0]*len(set_pos_list)
    totc_list = [0]*len(set_pos_list)
    
    for a in open("S1_methratio_CG_chr"+str(chr)+".txt", 'r') :
        metc_list[set_pos_list.index(int(a.split("\t")[1]))] += int(a.split("\t")[6])
        totc_list[set_pos_list.index(int(a.split("\t")[1]))] += int(a.split("\t")[7])
    for a in open("S2_methratio_CG_chr"+str(chr)+".txt", 'r') :
        metc_list[set_pos_list.index(int(a.split("\t")[1]))] += int(a.split("\t")[6])
        totc_list[set_pos_list.index(int(a.split("\t")[1]))] += int(a.split("\t")[7])
    for a in open("S3_methratio_CG_chr"+str(chr)+".txt", 'r') :
        metc_list[set_pos_list.index(int(a.split("\t")[1]))] += int(a.split("\t")[6])
        totc_list[set_pos_list.index(int(a.split("\t")[1]))] += int(a.split("\t")[7])
    
    for pos in set_pos_list : 
        print pos, metc_list[set_pos_list.index(pos)], totc_list[set_pos_list.index(pos)]
        sys.exit()
    