import os,sys

## make negative bed from negative_peaks.bed
outfile = open(sys.argv[1].split("_peaks.bed")[0]+"_negative_peaks.bed", 'w')
for a in open(sys.argv[1].split("_peaks.bed")[0]+"_negative_peaks.xls", 'r') :
    count = 1
    if a.split("\t")[0] != "chr" :
        outfile.write(a.split("\t")[0])
        outfile.write("\t")
        outfile.write(a.split("\t")[1])
        outfile.write("\t")
        outfile.write(a.split("\t")[2])
        outfile.write("\t")
        outfile.write("negative_peak"+str(count))
        outfile.write("\t")
        outfile.write(a.split("\t")[6])
        outfile.write("\n")
        count += 1
outfile.close()

## remove negative peaks from peaks_bed
chr_list = ['1','2','3','4','5','Mt','Pt']
## save locations by chromosome
outfile_2 = open(sys.argv[1].split("_peaks.bed")[0]+"_filtered_peaks.bed", 'w')
for chr in chr_list :
    print chr
    ## save negative peak location
    large_list = []
    for a in open(sys.argv[1].split("_peaks.bed")[0]+"_negative_peaks.bed",'r') :
        if str(a.split("\t")[0]) == chr :
            large_list.append(range(int(a.split("\t")[1]),int(a.split("\t")[2])))
    for b in open(sys.argv[1], 'r') :
        if str(b.split("\t")[0]) == chr :
            bed_list = range(int(b.split("\t")[1]),int(b.split("\t")[2]))
            count_overlap = 0
            for large in large_list :
                if len(set(bed_list) & set(large)) > 0 :
                    count_overlap += 1
            if count_overlap == 0 :
                outfile_2.write(b)
outfile.close()