import os,sys

window_size = 500000
chr = 1
max_count = (249250621/window_size)+1

def wmet() :
    large_list = []
    for x in range(0, max_count) :
        large_list.append([0]*2)
    all_list = []
    for a in open("SRR641628/CpG_OB_SRR641628_unsorted.coverage.bismark.cov", 'r') :
        if a.split("\t")[0] == "chr1" : 
            all_list.append(a)
    for a in open("SRR641628/CpG_OT_SRR641628_unsorted.coverage.bismark.cov", 'r') :
        if a.split("\t")[0] == "chr1" : 
            all_list.append(a)
    for a in all_list :    
        list_loc = int(a.split("\t")[1])/window_size
        if float(a.split("\t")[4])+float(a.split("\t")[5].strip("\n")) > 0 : 
            if float(a.split("\t")[4])/(float(a.split("\t")[4])+float(a.split("\t")[5].strip("\n"))) >= 0.126 : 
                large_list[list_loc][0] += int(a.split("\t")[4])
            large_list[list_loc][1] += float(a.split("\t")[4])+float(a.split("\t")[5].strip("\n"))
    outfile = open("SRR641628_chr"+str(chr)+"_wmet_"+str(window_size)+".bedGraph", 'w')
    for x in range(0, max_count) :
        if float(large_list[x][1]) > 0 :
            outfile.write(str(chr))
            outfile.write("\t")
            outfile.write(str(x*window_size+1))
            outfile.write("\t")
            outfile.write(str((x+1)*window_size))
            outfile.write("\t")
            outfile.write(str((float(large_list[x][0])/float(large_list[x][1]))*100))
            outfile.write("\n")
    outfile.close()
    outfile = open("SRR641628_chr"+str(chr)+"_wmet_"+str(window_size)+".circos", 'w')
    for x in range(0, max_count) :
        if float(large_list[x][1]) > 0 :
            outfile.write("hs"+str(chr))
            outfile.write("\t")
            outfile.write(str(x*window_size+1))
            outfile.write("\t")
            outfile.write(str((x+1)*window_size))
            outfile.write("\t")
            outfile.write(str((float(large_list[x][0])/float(large_list[x][1]))*100))
            outfile.write("\n")
    outfile.close()
print "wmet"
wmet()

