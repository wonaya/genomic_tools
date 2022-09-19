import os,sys


for chr in range(1,2) :
    coord_list = []
    value_list = []
    
    print "chr"+str(chr)
    start = []
    end = []
    for a in open("LS-merge_sorted_logFC.6.smooth.wig", 'r') :
        if int(a.split("\t")[0]) == chr :
            coord_list.append(int(a.split("\t")[1]))
            value_list.append(float(a.split("\t")[3].strip("\n")))
            if len(value_list) < 2 :
                if value_list[-1] >= 0 : 
                    start.append(1)
            if len(value_list) >= 2:
                if value_list[-2] < 0 and value_list[-1] > 0 : 
                    start.append(coord_list[-2]+((float(coord_list[-1])-float(coord_list[-2]))/2)-1)
                elif value_list[-2] > 0 and value_list[-1] < 0 : 
                    end.append(coord_list[-1]+((float(coord_list[-1])-float(coord_list[-2]))/2)-1)
    print len(start), len(end)
    if len(start) > len(end) :
        start = start[:-1]
        print len(start), len(end)
    outfile = open("LS_chr1.6.bedGraph", 'w')
    for no in range(0, len(start)) :
        outfile.write(str(chr)+"\t"+str(int(start[no]))+"\t"+str(int(end[no]))+"\t1\n")
    outfile.close()
    sys.exit()
"""
# spreading replicons    
for chr in range(1,11) :
    outfile = open("MS_LS_spreading.bedGraph", 'w')
    for a in open("MS_chr1.bedGraph", 'r') :
        for b in open("LS_chr1.bedGraph", 'r') :
            if int(b.split("\t")[1]) <= int(a.split("\t")[1]) and int(b.split("\t")[2]) >= int(a.split("\t")[2]) :
                outfile.write(b)
                
    outfile.close()
    sys.exit()
"""    
    