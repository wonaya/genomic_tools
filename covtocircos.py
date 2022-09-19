import os,sys

window_size = 100000
chr = 1
max_count = (301354095/window_size)+1

def average() :
    large_list = []
    for x in range(0, max_count) :
        large_list.append([0]*2)
    for a in open("CML322_all3/CpG_context_CML322_all3_bt202_sortedn_chr"+str(chr)+".bismark.cov", 'r') :
        list_loc = int(a.split("\t")[1])/window_size
        large_list[list_loc][0] += int(a.split("\t")[4])
        large_list[list_loc][1] += int(a.split("\t")[5].strip("\n"))

    outfile = open("CpG_context_CML322_all3_bt202_chr"+str(chr)+"_average_"+str(window_size)+".bedGraph", 'w')
    for x in range(0, max_count) :
        if (float(large_list[x][0])+float(large_list[x][1]))*100 > 0 :
            outfile.write(str(chr))
            outfile.write("\t")
            outfile.write(str(x*window_size+1))
            outfile.write("\t")
            outfile.write(str((x+1)*window_size))
            outfile.write("\t")
            outfile.write(str(float(large_list[x][0])/(float(large_list[x][0])+float(large_list[x][1]))*100))
            outfile.write("\n")
    outfile.close()
    
    outfile = open("CpG_context_CML322_all3_bt202_chr"+str(chr)+"_average_"+str(window_size)+".circos", 'w')
    for x in range(0, max_count) :
        if (float(large_list[x][0])+float(large_list[x][1]))*100 > 0 :
            outfile.write("zm"+str(chr))
            outfile.write("\t")
            outfile.write(str(x*window_size+1))
            outfile.write("\t")
            outfile.write(str((x+1)*window_size))
            outfile.write("\t")
            outfile.write(str(float(large_list[x][0])/(float(large_list[x][0])+float(large_list[x][1]))*100))
            outfile.write("\n")
    outfile.close()
    

def fraction():
    ## fraction
    large_list = []
    for x in range(0, max_count) :
        large_list.append([0]*2)
    for a in open("CML322_all3/CpG_context_CML322_all3_bt202_sortedn_chr"+str(chr)+".bismark.cov", 'r') :
        list_loc = int(a.split("\t")[1])/window_size
        if float(a.split("\t")[4])/(float(a.split("\t")[4])+float(a.split("\t")[5].strip("\n"))) >= 0.126 :
            large_list[list_loc][0] += 1
        large_list[list_loc][1] += 1
    outfile = open("CpG_context_CML322_all3_bt202_chr"+str(chr)+"_frac_"+str(window_size)+".bedGraph", 'w')
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
    outfile = open("CpG_context_CML322_all3_bt202_chr"+str(chr)+"_frac_"+str(window_size)+".circos", 'w')
    for x in range(0, max_count) :
        if float(large_list[x][1]) > 0 :
            outfile.write("zm"+str(chr))
            outfile.write("\t")
            outfile.write(str(x*window_size+1))
            outfile.write("\t")
            outfile.write(str((x+1)*window_size))
            outfile.write("\t")
            outfile.write(str((float(large_list[x][0])/float(large_list[x][1]))*100))
            outfile.write("\n")
    outfile.close()


def mean():
    large_list = []
    for x in range(0, max_count) :
        large_list.append([0]*2)
    for a in open("CML322_all3/CpG_context_CML322_all3_bt202_sortedn_chr"+str(chr)+".bismark.cov", 'r') :
        list_loc = int(a.split("\t")[1])/window_size
        if float(a.split("\t")[4])+float(a.split("\t")[5].strip("\n")) > 0 : 
            if float(a.split("\t")[4])/(float(a.split("\t")[4])+float(a.split("\t")[5].strip("\n"))) >= 0.126 :
                large_list[list_loc][0] += float(a.split("\t")[4])/(float(a.split("\t")[4])+float(a.split("\t")[5].strip("\n")))
            large_list[list_loc][1] += 1
    outfile = open("CpG_context_CML322_all3_bt202_chr"+str(chr)+"_mean_"+str(window_size)+".bedGraph", 'w')
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
    outfile = open("CpG_context_CML322_all3_bt202_chr"+str(chr)+"_mean_"+str(window_size)+".circos", 'w')
    for x in range(0, max_count) :
        if float(large_list[x][1]) > 0 :
            outfile.write("zm"+str(chr))
            outfile.write("\t")
            outfile.write(str(x*window_size+1))
            outfile.write("\t")
            outfile.write(str((x+1)*window_size))
            outfile.write("\t")
            outfile.write(str((float(large_list[x][0])/float(large_list[x][1]))*100))
            outfile.write("\n")
    outfile.close()
    
def wmet() :
    large_list = []
    for x in range(0, max_count) :
        large_list.append([0]*2)
    for a in open("CML322_all3/CpG_context_CML322_all3_bt202_sortedn_chr"+str(chr)+".bismark.cov", 'r') :
        list_loc = int(a.split("\t")[1])/window_size
        if float(a.split("\t")[4])+float(a.split("\t")[5].strip("\n")) > 0 : 
            if float(a.split("\t")[4])/(float(a.split("\t")[4])+float(a.split("\t")[5].strip("\n"))) >= 0.126 : 
                large_list[list_loc][0] += int(a.split("\t")[4])
            large_list[list_loc][1] += float(a.split("\t")[4])+float(a.split("\t")[5].strip("\n"))
    outfile = open("CpG_context_CML322_all3_bt202_chr"+str(chr)+"_wmet_"+str(window_size)+".bedGraph", 'w')
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
    outfile = open("CpG_context_CML322_all3_bt202_chr"+str(chr)+"_wmet_"+str(window_size)+".circos", 'w')
    for x in range(0, max_count) :
        if float(large_list[x][1]) > 0 :
            outfile.write("zm"+str(chr))
            outfile.write("\t")
            outfile.write(str(x*window_size+1))
            outfile.write("\t")
            outfile.write(str((x+1)*window_size))
            outfile.write("\t")
            outfile.write(str((float(large_list[x][0])/float(large_list[x][1]))*100))
            outfile.write("\n")
    outfile.close()
print "average"
average()
print "fraction"
fraction()
print "mean"
mean()
print "wmet"
wmet()

