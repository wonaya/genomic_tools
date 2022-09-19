import os,sys
import collections

## use macs_peak_compare.py control.peak.bed file1 file2 ... fileN

## save all gene coordinates
### get into one file 
dict = {}
chr_list = ['Mt', 'Pt','5', '4', '3', '2', '1']
outfile = open("macs_output.txt", 'w')
list_dict = []
for chr in chr_list :
    print chr
    dict = {}
    for a in open("/work/02114/wonaya/genome/annotation/Arabidopsis_thaliana.TAIR10.19.gtf", 'r') :
        if a.split("\t")[0] == chr and a.split("\t")[1] == "protein_coding" and a.split("\t")[2] == "exon" :
            dict[a.split("\t")[8].split(";")[0].split('gene_id "')[1].strip('"')] = []
    for a in open("/work/02114/wonaya/genome/annotation/Arabidopsis_thaliana.TAIR10.19.gtf", 'r') :
        if a.split("\t")[0] == chr and a.split("\t")[1] == "protein_coding" and a.split("\t")[2] == "exon" :
            dict[a.split("\t")[8].split(";")[0].split('gene_id "')[1].strip('"')].append(int(a.split("\t")[3]))
            dict[a.split("\t")[8].split(";")[0].split('gene_id "')[1].strip('"')].append(int(a.split("\t")[4]))
    final_dict = {}
    for b in dict.keys() :
        if len(dict[b]) > 2 : 
            final_dict[b] = range(min(dict[b]), max(dict[b])+1)
        else :
            final_dict[b] = range(dict[b][0], dict[b][1]+1)
        
    ## get bed from chromosome
    result_dict1 = {}
    for gene in final_dict.keys() :
        score = 0
        for c in open(sys.argv[1], 'r') :
            if str(c.split("\t")[0]) == str(chr) :  
                if len(set(range(int(c.split("\t")[1]),int(c.split("\t")[2])+1)) & set(final_dict[gene])) > 0 :
                    score += float(c.split("\t")[4].strip("\n"))
        result_dict1[gene] = score

    result_dict2 = {}
    for gene in final_dict.keys() :
        score = 0
        for c in open(sys.argv[2], 'r') :
            if str(c.split("\t")[0]) == str(chr) :  
                if len(set(range(int(c.split("\t")[1]),int(c.split("\t")[2])+1)) & set(final_dict[gene])) > 0 :
                    score += float(c.split("\t")[4].strip("\n"))
        result_dict2[gene] = score

    result_dict3 = {}
    for gene in final_dict.keys() :
        score = 0
        for c in open(sys.argv[3], 'r') :
            if str(c.split("\t")[0]) == str(chr) :  
                if len(set(range(int(c.split("\t")[1]),int(c.split("\t")[2])+1)) & set(final_dict[gene])) > 0 :
                    score += float(c.split("\t")[4].strip("\n"))
        result_dict3[gene] = score
    
    result_dict4 = {}
    for gene in final_dict.keys() :
        score = 0
        for c in open(sys.argv[4], 'r') :
            if str(c.split("\t")[0]) == str(chr) :  
                if len(set(range(int(c.split("\t")[1]),int(c.split("\t")[2])+1)) & set(final_dict[gene])) > 0 :
                    score += float(c.split("\t")[4].strip("\n"))
        result_dict4[gene] = score

    result_dict5 = {}
    for gene in final_dict.keys() :
        score = 0
        for c in open(sys.argv[5], 'r') :
            if str(c.split("\t")[0]) == str(chr) :  
                if len(set(range(int(c.split("\t")[1]),int(c.split("\t")[2])+1)) & set(final_dict[gene])) > 0 :
                    score += float(c.split("\t")[4].strip("\n"))
        result_dict5[gene] = score

    result_dict6 = {}
    for gene in final_dict.keys() :
        score = 0
        for c in open(sys.argv[6], 'r') :
            if str(c.split("\t")[0]) == str(chr) :  
                if len(set(range(int(c.split("\t")[1]),int(c.split("\t")[2])+1)) & set(final_dict[gene])) > 0 :
                    score += float(c.split("\t")[4].strip("\n"))
        result_dict6[gene] = score

    result_dict7 = {}
    for gene in final_dict.keys() :
        score = 0
        for c in open(sys.argv[7], 'r') :
            if str(c.split("\t")[0]) == str(chr) :  
                if len(set(range(int(c.split("\t")[1]),int(c.split("\t")[2])+1)) & set(final_dict[gene])) > 0 :
                    score += float(c.split("\t")[4].strip("\n"))
        result_dict7[gene] = score

    result_dict8 = {}
    for gene in final_dict.keys() :
        score = 0
        for c in open(sys.argv[8], 'r') :
            if str(c.split("\t")[0]) == str(chr) :  
                if len(set(range(int(c.split("\t")[1]),int(c.split("\t")[2])+1)) & set(final_dict[gene])) > 0 :
                    score += float(c.split("\t")[4].strip("\n"))
        result_dict8[gene] = score

    for results in result_dict1 :
        print results, int(result_dict1[results]), int(result_dict2[results]), int(result_dict3[results]), int(result_dict4[results]), int(result_dict5[results]), int(result_dict6[results]), int(result_dict7[results]), int(result_dict8[results])
    
        outfile.write(str(results))
        outfile.write("\t")
        outfile.write(str(int(result_dict1[results])))
        outfile.write("\t")
        outfile.write(str(int(result_dict2[results])))
        outfile.write("\t")
        outfile.write(str(int(result_dict3[results])))
        outfile.write("\t")
        outfile.write(str(int(result_dict4[results])))
        outfile.write("\t")
        outfile.write(str(int(result_dict5[results])))
        outfile.write("\t")
        outfile.write(str(int(result_dict6[results])))
        outfile.write("\t")
        outfile.write(str(int(result_dict7[results])))
        outfile.write("\t")
        outfile.write(str(int(result_dict8[results])))
        outfile.write("\n")
outfile.close()
sys.exit()
"""
nooffiles = len(sys.argv)

a0 = open(sys.argv[1],'r')
a0lines = a0.readlines()

chr_list = []
for a0line in a0lines :
    if a0line.split("\t")[0] not in chr_list :
        chr_list.append(a0line.split("\t")[0])

peaks_large = []
for chr in chr_list :
    peaks = []
    for a0line in a0lines :
        if a0line.split("\t")[0] == chr :
            peaks.append(a0line.split("\t")[1]+"-"+a0line.split("\t")[2])
    peaks_large.append(peaks)
    
a1 = open(sys.argv[2],'r')
a1lines = a1.readlines()

test_peaks_large = []
for chr in chr_list :
    test_peaks = []
    for a1line in a1lines :
        if a1line.split("\t")[0] == chr :
            test_peaks.append(a1line.split("\t")[1]+"-"+a1line.split("\t")[2])
    test_peaks_large.append(test_peaks)

for x in range(0,13) :
    count = 0
    for peak in peaks_large[x] :
        for test_peak in test_peaks_large[x] :
            if int(test_peak.split("-")[0]) < int(peak.split("-")[0]) :
                if int(test_peak.split("-")[1]) > int(peak.split("-")[0]) and int(test_peak.split("-")[1]) < int(peak.split("-")[1]) :
                    count += 1
                elif int(test_peak.split("-")[1]) > int(peak.split("-")[1]) :
                    count += 1
            elif int(test_peak.split("-")[0]) > int(peak.split("-")[0]) :
                if int(test_peak.split("-")[0]) < int(peak.split("-")[1]) and int(test_peak.split("-")[1]) > int(peak.split("-")[1]) :
                    count += 1
                if int(test_peak.split("-")[1]) < int(peak.split("-")[1]) :
                    count += 1
    print "chr:",chr_list[x], "control:",len(peaks_large[x]), "test:",len(test_peaks_large[x]), "overlap count:",count
"""

def overlap_by_chr():
    ## use macs_peak_compare.py control.bedfile, control.bedgrah file1,file2...fileN
    if len(sys.argv) < 6 :
        print "usage python /work/02114/wonaya/scripts/macs_peak_compare.py 1_S1.unamp_peaks.bed 1_S1.unamp.sorted.bedGraph 1_S1.ds.bedGraph 2_S2.ds.bedGraph 3_S3.ds.bedGraph 4_S4.ds.bedGraph (0 to 12)"
        sys.exit()
    a0 = open(sys.argv[1],'r')
    a0lines = a0.readlines()

    ## peak regions from control
    chr_list = [] # list of chromosomes
    for a0line in a0lines :
        if a0line.split("\t")[0] not in chr_list :
            chr_list.append(a0line.split("\t")[0])
    chr_list.sort()

    peaks_starts = []
    peaks_ends   = []
    for a0line in a0lines :
        if a0line.split("\t")[0] == chr_list[int(sys.argv[7])] :
            peaks_ends.append(int(a0line.split("\t")[2]))
            peaks_starts.append(int(a0line.split("\t")[1]))
    x = 0
    peaks_starts.sort()
    peaks_dict = {}
    while x <= (int(max(peaks_ends))/10000000)+1 :
    #while x <= 0 :
        ## peak values from control bedgraph 
        peaks_start = []
        peaks_end   = []
        for a0line in a0lines :
            if a0line.split("\t")[0] == chr_list[int(sys.argv[7])] and int(a0line.split("\t")[1]) >= x*10000000 and int(a0line.split("\t")[1]) < (x+1)*10000000:
                peaks_start.append(int(a0line.split("\t")[1]))
                peaks_end.append(int(a0line.split("\t")[2]))
                peaks_dict[int(a0line.split("\t")[1])] = [int(a0line.split("\t")[1]),int(a0line.split("\t")[2])]
        print "total peaks of chr", chr_list[int(sys.argv[7])], "between", x*10000000, "and", (x+1)*10000000, "are", len(peaks_end)
        #if len(peaks_end) == 0 :
        #    print "no more peaks between", x*10000000, "and", (x+1)*10000000, ": break"
        #    break
        #del a0line, a0lines
        print "reading control file:", sys.argv[2]
        bed_large1 = []
        a1 = open(sys.argv[2], 'r') 
        a1lines = a1.readlines()
        for a1line in a1lines :
            if a1line.split("\t")[0] == chr_list[int(sys.argv[7])] and int(a1line.split("\t")[1]) >= x*10000000 and int(a1line.split("\t")[1]) < (x+1)*10000000 :  
                bed_large1.append(a1line)
        #del a1line, a1lines
        #print "counting", sys.argv[2]
        for peak in peaks_start :
            total = 0
            for a1line in bed_large1 :
                a_multiset = range(int(a1line.split("\t")[1]),int(a1line.split("\t")[2]))
                b_multiset = range(peak,peaks_end[peaks_start.index(peak)])
                overlap = set(a_multiset).intersection(set(b_multiset))
                if len(overlap) > 0 :
                    total += int(a1line.split("\t")[3].strip("\n"))
            peaks_dict[peak].append(total)
        del bed_large1
        print "reading sample1 file:", sys.argv[3]    
        bed_large2 = []
        a2 = open(sys.argv[3], 'r') 
        a2lines = a2.readlines()
        for a2line in a2lines :
            if a2line.split("\t")[0] == chr_list[int(sys.argv[7])] and int(a2line.split("\t")[1]) >= x*10000000 and int(a2line.split("\t")[1]) < (x+1)*10000000 :  
                bed_large2.append(a2line)
        #del a2line, a2lines
        #print "counting", sys.argv[3]
        for peak in peaks_start :
            total = 0
            for a2line in bed_large2 :
                a_multiset = range(int(a2line.split("\t")[1]),int(a2line.split("\t")[2]))
                b_multiset = range(peak,peaks_end[peaks_start.index(peak)])
                overlap = set(a_multiset).intersection(set(b_multiset))
                if len(overlap) > 0 :
                    total += int(a2line.split("\t")[3].strip("\n"))
            peaks_dict[peak].append(total)
        del bed_large2
        print "reading sample2 file:", sys.argv[4]    
        bed_large3 = []
        a3 = open(sys.argv[4], 'r') 
        a3lines = a3.readlines()
        for a3line in a3lines :
            if a3line.split("\t")[0] == chr_list[int(sys.argv[7])] and int(a3line.split("\t")[1]) >= x*10000000 and int(a3line.split("\t")[1]) < (x+1)*10000000 :  
                bed_large3.append(a3line)
        #del a3line, a3lines
        #print "counting", sys.argv[4]
        for peak in peaks_start :
            total = 0
            for a3line in bed_large3 :
                a_multiset = range(int(a3line.split("\t")[1]),int(a3line.split("\t")[2]))
                b_multiset = range(peak,peaks_end[peaks_start.index(peak)])
                overlap = set(a_multiset).intersection(set(b_multiset))
                if len(overlap) > 0 :
                    total += int(a3line.split("\t")[3].strip("\n"))
            peaks_dict[peak].append(total)
        del bed_large3
    
        print "reading sample3 file:", sys.argv[5]    
        bed_large4 = []
        a4 = open(sys.argv[5], 'r') 
        a4lines = a4.readlines()
        for a4line in a4lines :
            if a4line.split("\t")[0] == chr_list[int(sys.argv[7])] and int(a4line.split("\t")[1]) >= x*10000000 and int(a4line.split("\t")[1]) < (x+1)*10000000 :  
                bed_large4.append(a4line)
        #del a4line, a4lines
        #print "counting", sys.argv[5]
        for peak in peaks_start :
            total = 0
            for a4line in bed_large4 :
                a_multiset = range(int(a4line.split("\t")[1]),int(a4line.split("\t")[2]))
                b_multiset = range(peak,peaks_end[peaks_start.index(peak)])
                overlap = set(a_multiset).intersection(set(b_multiset))
                if len(overlap) > 0 :
                    total += int(a4line.split("\t")[3].strip("\n"))
            peaks_dict[peak].append(total)
        del bed_large4
    
        print "reading sample4 file:", sys.argv[6]    
        bed_large5 = []
        a5 = open(sys.argv[6], 'r') 
        a5lines = a5.readlines()
        for a5line in a5lines :
            if a5line.split("\t")[0] == chr_list[int(sys.argv[7])] and int(a5line.split("\t")[1]) >= x*10000000 and int(a5line.split("\t")[1]) < (x+1)*10000000 :  
                bed_large5.append(a5line)
        #del a5line, a5lines
        #print "counting", sys.argv[6]
        for peak in peaks_start :
            total = 0
            for a5line in bed_large5 :
                a_multiset = range(int(a5line.split("\t")[1]),int(a5line.split("\t")[2]))
                b_multiset = range(peak,peaks_end[peaks_start.index(peak)])
                overlap = set(a_multiset).intersection(set(b_multiset))
                if len(overlap) > 0 :
                    total += int(a5line.split("\t")[3].strip("\n"))
            peaks_dict[peak].append(total)
        del bed_large5
        x += 1
    print "writing file"
    outfile = open(str(sys.argv[1]).split(".bed")[0]+"_"+str(chr_list[int(sys.argv[7])])+".overlap", 'w')
    #outfile.write("Peak start\tPeak end\tControl count\tSample1\tSample2\tSample3\tSample4\n")
    for peak in peaks_starts:
        #print peaks_dict[peak]
        outfile.write(chr_list[int(sys.argv[7])])
        outfile.write("\t")
        outfile.write("\t".join([str(x) for x in peaks_dict[peak]]))
        outfile.write("\n")
    outfile.close()     

def merge_all():
    chr_list = [1,2,3,4,5,6,7,8,9,10,'chloroplast','mitochondrion','UNKNOWN']
    outfile = open(str(sys.argv[1]).split(".bed")[0]+".overlap", 'w')
    outfile.write("Chr\tPeak start\tPeak end\tControl count\tSample1\tSample2\tSample3\tSample4\n")
    for chr in chr_list :
        for line in open("1_S1.unamp-bwa-pe-macs1M_peaks_"+str(chr)+".overlap", 'r') :
            outfile.write(line)
    outfile.close()
def main():
    #overlap_by_chr()
    merge_all()
    
if __name__ == '__main__':
    main()    