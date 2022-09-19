import os,sys

## usage python c_context_window.py CHH_context.out chromosme_index(0to9) "split = if split required"
# .py infile (methylation extractor output), window size, minimum count, "label"

if sys.argv[3] == "split" :
    ## divide by chr names
    chr_list = ['chr1','chr2','chr3','chr4','chr5','chr6', 'chr7', 'chr8','chr9','chr10']
    outfile = open(str(sys.argv[1]).split(".")[0]+"_"+chr_list[int(sys.argv[2])]+".out", 'w')
    for a in open(sys.argv[1], 'r') :
        if len(a.split("\t")) > 2 and str(a.split("\t")[2]) == chr_list[int(sys.argv[2])] :  
            outfile.write(a) 
    outfile.close()
    sys.exit()    

sys.exit()
## get list of chromosome names
#chr_list = []
#for a0 in open(sys.argv[1], 'r') :
#    if a0[:7] != "Bismark" and a0.split("\t")[2] not in chr_list :
#        chr_list.append(a0.split("\t")[2])
#chr_list.sort()

chr = chr_list[int(sys.argv[2])]
# run by chromosomes
outfile = open(str(sys.argv[1])+"_wMet_calc_"+chr_list[int(sys.argv[2])]+".out", 'w')
outfile.write("chrom\tstart\tend\tNo.Me+reads in region\tNo.Me+reads>threshold\tNo.Me+reads/All reads\tCount of methylated C over threshold\tTotal no. of reads\tTotal C's\tFractional\tMean\tWeighted\n") 

positions = []
print chr
for a0 in open(sys.argv[1], 'r') :
    if a0[:7] != "Bismark" and a0.split("\t")[2] == chr :
        positions.append(int(a0.split("\t")[3])) # get all the positions in the file for that chromosome

x = 0
#while x <= (int(max(positions))/10000000)+2 : ## cut positions into 10M portions
while x <= 10: ## max should be 32, 
    print str(x)+"th iteration"
    locations = []
    for a0 in open(sys.argv[1], 'r') :
        if a0[:7] != "Bismark" and a0.split("\t")[2] == chr :
            if int(a0.split("\t")[3]) >= x*10000000 and int(a0.split("\t")[3]) < (x+1)*10000000 and len(a0.split("\t")) == 5 : # 10M at a time
                locations.append(a0.strip("\n"))
    print len(locations)
    if len(locations) == 0 :
        print "end of file, closing file"
        break
    huge_list = []
    for y in range(0, 10000000) :
        huge_list.append([])
    for a in locations :
        huge_list[int(a.split("\t")[3])-(x*10000000)].append(a.split("\t")[1])
            
    ### analysis part
    for count in range(0,100000) :
        y = 0
        total = 0
        total_count = 0
        met = 0
        valid_met = 0
        valid_met_val = 0
        valid_met_count = 0
        for b in huge_list[count*100:(count+1)*100]: ## window size
            if len(b) > 0 :
                total += len(b)
                total_count += 1
                met += b.count("+")
                if float(b.count("+"))/float(len(b)) >= 0.126 :
                    valid_met += int(b.count("+"))
                    valid_met_val += float(b.count("+"))/float(len(b))
                    valid_met_count += 1
            y += 1
        
        #print count*100,(count+1)*100, met, valid_met, valid_met_val, valid_met_count, total, total_count
        if total > 0 :
            #print count*100,(count+1)*100, met, valid_met, valid_met_val, valid_met_count, total, total_count, float(valid_met_count)/float(total_count), float(valid_met_val)/float(total_count), float(valid_met)/float(total) 
            outfile.write(str(chr))
            outfile.write("\t")
            outfile.write(str((count*100)+(x*10000000))+"\t"+str(((count+1)*100)+(x*10000000)))
            outfile.write("\t")
            outfile.write(str(met))
            outfile.write("\t")
            outfile.write(str(valid_met))
            outfile.write("\t")
            outfile.write(str(valid_met_val))
            outfile.write("\t")
            outfile.write(str(valid_met_count))
            outfile.write("\t")
            outfile.write(str(total))
            outfile.write("\t")
            outfile.write(str(total_count))
            outfile.write("\t")
            outfile.write(str(float(valid_met_count)/float(total_count)))
            outfile.write("\t")
            outfile.write(str(float(valid_met_val)/float(total_count)))
            outfile.write("\t")
            outfile.write(str(float(valid_met)/float(total)))
            outfile.write("\n")
        count += 1
    x += 1
outfile.close()    
   