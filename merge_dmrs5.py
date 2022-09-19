import os,sys
import matplotlib.pyplot as plt
import math
## density plot first
cg_distance = []
chg_distance = []
    
for chr in range(1,11) :
    cg_dmr_list = []
    afile = open("test_chr1_distance43.txt", 'r') 
    alines = afile.readlines()
    for a in alines[1:]:
        ## save all locations
        if float(a.split("\t")[9]) <= 0.05 and a.split("\t")[1] == '"'+str(chr)+'"' :
            cg_dmr_list.append(range(int(a.split("\t")[2]),int(a.split("\t")[3])+1))
    chg_dmr_list = []
    afile = open("B73_Mo17_CHG_dmr_1000.txt", 'r') 
    alines = afile.readlines()
    for a in alines[1:]:
        ## save all locations
        if float(a.split("\t")[8]) <= 0.05 and a.split("\t")[0] == str(chr) :
            chg_dmr_list.append(range(int(a.split("\t")[1]),int(a.split("\t")[2])+1))
    print len(cg_dmr_list),len(chg_dmr_list)
    x = 0
    for cg_dmr in cg_dmr_list : 
        if x >= 1 :
            cg_distance.append(math.log(cg_dmr[0]-cg_dmr_list[x-1][-1],10))
            #cg_distance.append(cg_dmr[0]-cg_dmr_list[x-1][-1])
        x += 1
    x = 0
    for chg_dmr in chg_dmr_list : 
        if x >= 1 :
            chg_distance.append(math.log(chg_dmr[0]-chg_dmr_list[x-1][-1],10))
        x += 1

fig, ax = plt.subplots()
#ax.set_yscale('log',basey=2)
n, bins, patches = ax.hist(cg_distance, 50, facecolor='green', alpha=0.75)
plt.xlabel('Log 10 Distance')
plt.ylabel('Frequency')
#plt.xlim([0,10])
#plt.yscale('log')
#plt.show()
plt.savefig("CpG_B73_Mo17_dist43_eDMR_DMR.png")
sys.exit()
for chr in range(1,11) :
    cg_dmr_list = []
    afile = open("B73_all3_Mo17_all3_CpG_dmr_merged.txt", 'r') 
    alines = afile.readlines()
    for a in alines[1:]:
        ## save all locations
        if float(a.split("\t")[7]) <= 0.05 and a.split("\t")[0] == "chr"+str(chr) :
            cg_dmr_list.append(range(int(a.split("\t")[1]),int(a.split("\t")[2])+1))
    chg_dmr_list = []
    afile = open("B73_all3_Mo17_all3_CHG_dmr_merged.txt", 'r') 
    alines = afile.readlines()
    for a in alines[1:]:
        ## save all locations
        if float(a.split("\t")[7]) <= 0.05 and a.split("\t")[0] == "chr"+str(chr) :
            chg_dmr_list.append(range(int(a.split("\t")[1]),int(a.split("\t")[2])+1))
    print len(cg_dmr_list),len(chg_dmr_list)
    x = 0
    
