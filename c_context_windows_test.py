import os,sys

# .py infile (methylation extractor output), window size, minimum count, "label"
def print_chr_name(data) :
    # get list of chromosome names
    chr_list = []
    for a0 in open(sys.argv[1], 'r') :
        if a0[:7] != "Bismark" and a0.split("\t")[2] not in chr_list :
            chr_list.append(a0.split("\t")[2])
    chr_list.sort()
    print chr_list
    return chr_list

def separate_by_chr(data, chr) :
    # run by chromosomes
    positions = []
    outfile = open(str(data).strip(".txt")+"_"+str(chr)+".out", 'w')
    for a0 in open(sys.argv[1], 'r') :
        if a0[:7] != "Bismark" and a0.split("\t")[2] == str(chr) :
            outfile.write(a0)
            #positions.append(int(a0.split("\t")[3])) # get all the positions in the file for that chromosome
    outfile.close()
    #return positions

def calculation(data, chr) :
    #while x <= (int(max(positions))/10000000)+2 : ## cut positions into 10M portions
    #while x <= 32: ## max should be 32, 
    if not os.path.isfile(str(data).strip(".txt")+"_"+str(chr)+".out") :
        print "chr-specific file not there, running separation script"
        separate_by_chr(data, chr)
    x = 0
    result_list = []
    for x in range(0, 1000000): ## 1M iteration in total
        locations = []
        for a0 in open(str(data).strip(".txt")+"_"+str(chr)+".out", 'r')  :
            if int(a0.split("\t")[3]) >= x*10000000 and int(a0.split("\t")[3]) < (x+1)*10000000 and len(a0.split("\t")) == 5 : # 10M at a time
                locations.append(a0.strip("\n"))
        print str(x)+"th iteration, no.of C is",len(locations), "between", x*10000000, "to", (x+1)*10000000
        if len(locations) == 0 :
            print "end of file, closing file"
            break
        #huge_list = []
        #for y in range(0, 10000000) : # making list for each C site
        #    huge_list.append([])
        huge_list = [[] for _ in range(10000000)]
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
                result_list.append(str(chr)+"\t"+str((count*100)+(x*10000000))+"\t"+str(((count+1)*100)+(x*10000000))+"\t"+str(met)+"\t"+str(valid_met)+"\t"+str(valid_met_val)+"\t"+str(valid_met_count)+"\t"+str(total)+"\t"+str(total_count)+"\t"+str(float(valid_met_count)/float(total_count))+"\t"+str(float(valid_met_val)/float(total_count))+"\t"+str(float(valid_met)/float(total))+"\n")
            count += 1
        
        x += 1
    outfile = open(str(data).strip(".txt")+"_wMet_calc_"+str(chr)+".txt", 'w')
    #outfile.write("chrom\tstart\tend\tNo.Me+reads in region\tNo.Me+reads>threshold\tNo.Me+reads/All reads\tCount of methylated C over threshold\tTotal no. of reads\tTotal C's\tFractional\tMean\tWeighted\n") 
    for result in result_list :
        outfile.write(result)

    outfile.close()    
def merge(data):
    file_list = []
    for files in os.listdir("."):
        if files.startswith(str(data).strip(".txt")+"_wMet_calc_chr") and files.endswith(".txt"):
        #if files.startswith(str(data)+"_wMet_calc_chr") and files.endswith(".txt"):
           file_list.append(files)
    outfile = open(str(data).strip(".txt")+"_wMet_calc_combined.txt", 'w')
    outfile.write("chrom\tstart\tend\tNo.Me+reads in region\tNo.Me+reads>threshold\tNo.Me+reads/All reads\tCount of methylated C over threshold\tTotal no. of reads\tTotal C's\tFractional\tMean\tWeighted\n") 
    outfile.close()
    chr_order = [1,2,3,4,5,6,7,8,9,10,'Mt']
    for chr in chr_order :
        os.system("cat "+str(data).strip(".txt")+"_wMet_calc_chr"+str(chr)+".txt >> "+str(data).strip(".txt")+"_wMet_calc_combined.txt")
        #os.system("cat "+str(data)+"_wMet_calc_chr"+str(chr)+".txt >> "+str(data).strip(".txt")+"_wMet_calc_combined.txt")
    os.system("rm -Rf "+str(data).strip(".txt")+"_wMet_calc_chr*")
def main():
    #print_chr_name(sys.argv[1])
    #separate_by_chr(sys.argv[1], sys.argv[2])
    calculation(sys.argv[1], sys.argv[2])
    #merge(sys.argv[1])
if __name__ == "__main__":
    main()
  