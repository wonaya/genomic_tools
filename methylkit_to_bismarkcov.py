import os,sys

### chr10

"""
outfile = open("CpG_context_Tx303_all3_bt202_chr10.bismark.cov", 'w')
for a in open("CpG_Tx303_all3_bt202.sam_chr10.methylKit", 'r') :
    outfile.write("chr"+a.split("\t")[1]+"\t")
    outfile.write(a.split("\t")[0]+"\t"+a.split("\t")[0]+"\t")
    outfile.write(a.split("\t")[5]+"\t")
    metc = float(a.split("\t")[4])*(float(a.split("\t")[5])/100)
    unmetc = float(a.split("\t")[4])*(float(a.split("\t")[6].strip("\n"))/100)
    outfile.write(str(int(metc))+"\t"+str(int(unmetc))+"\n")
outfile.close()
"""

large_list = []
for x in range(0,151000000):
    large_list.append([0]*2)
    
for a in open("CpG_OB_Tx303_all3_bt202_10.txt", 'r') :
    if a.split("\t")[1] == "+" :
        large_list[int(a.split("\t")[3])][0] += 1
    elif a.split("\t")[1] == "-" :
        large_list[int(a.split("\t")[3])][1] += 1

for a in open("CpG_OT_Tx303_all3_bt202_10.txt", 'r') :
    if a.split("\t")[1] == "+" :
        large_list[int(a.split("\t")[3])][0] += 1
    elif a.split("\t")[1] == "-" :
        large_list[int(a.split("\t")[3])][1] += 1

"""
for a in open("CpG_CTOB_Tx303_all3_bt202_10.txt", 'r') :
    if a.split("\t")[1] == "+" :
        large_list[int(a.split("\t")[3])][0] += 1
    elif a.split("\t")[1] == "-" :
        large_list[int(a.split("\t")[3])][1] += 1

for a in open("CpG_CTOT_Tx303_all3_bt202_10.txt", 'r') :
    if a.split("\t")[1] == "+" :
        large_list[int(a.split("\t")[3])][0] += 1
    elif a.split("\t")[1] == "-" :
        large_list[int(a.split("\t")[3])][1] += 1
"""
x = 0
for val in large_list :
    if val != [0,0] :
        print x, val
        sys.exit()
    x += 1