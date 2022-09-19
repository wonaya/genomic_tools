import os,sys

# Use python DMR_calling.py SAM_ID
# Need to extract info. from methylation extractor output
# merge by each chromosome
# order them with position/chr

def separate_chr(name, type):
    if type == "CpG" :
        methylation_type = ["CpG"]
    
    for chrom in range(1,2):
        outfile = open(str(name)+"_chr"+str(chrom)+".txt", 'w')
        list = []
        for a1line in open(str(type)+"_OT_"+name+".sorted.sam.txt", 'r') :
            if a1line[:7] != "Bismark" and int(a1line.split("\t")[2].split("chr")[1]) == chrom:
                list.append(a1line)
        for lis in list :
            outfile.write(lis)
        del list
        list = []
        for a2line in open(str(type)+"_OB_"+name+".sorted.sam.txt", 'r') :
            if a1line[:7] != "Bismark" and int(a1line.split("\t")[2].split("chr")[1]) == chrom:
                list.append(a2line)
        for lis in list :
            outfile.write(lis)
        del list
        outfile.close()
    
def merge_sort_chr(name, type,chrom):
    import glob
    outfile=name+"_merge_"+type+"_chr"+chrom+".out"
    if not os.path.isfile(name+"_merge_"+type+"_chr"+chrom+".out") :
        for files in glob.glob(name+"*chr"+chrom+".txt"):
            os.system("cat "+str(files)+" >> "+str(outfile))
    if not os.path.isfile(name+"_merge_"+type+"_chr"+chrom+"_sort.out") :
        sortfile=name+"_merge_"+type+"_chr"+chrom+"_sort.out"
        os.system("sort -nk4 -T $SCRATCH "+str(outfile)+" > "+str(sortfile))
        os.system("rm "+str(outfile))
def swDMR_input(name,type,chrom) :
    base = []
    stype = []
    strand = []
    outfile = open(name+"_swDMR_"+str(type)+"_chr"+str(chrom)+".txt", 'w')
    a = open(name+"_merge_"+str(type)+"_chr"+str(chrom)+"_sort.out", 'r')
    alines = a.readlines()
    start = int(alines[0].split("\t")[3])
    end  = int(alines[-1].split("\t")[3])
    for a in alines :
        if int(a.split("\t")[3]) == start:
            base.append(a.split("\t")[1])
        else :
            if len(base) > 1 : # take points where two more more methylations per point
                if len(set(strand)) > 1 :
                    print "some error in strand", strand
                    sys.exit()
                outfile.write("chr"+str(chrom))
                outfile.write("\t")
                outfile.write(str(start))
                outfile.write("\t")
                outfile.write("+")
                outfile.write("\t")
                outfile.write(str(type))
                outfile.write("\t")
                outfile.write(str(base.count("+")))
                outfile.write("\t")
                outfile.write(str(base.count("-")))
                outfile.write("\n")
                
            base = []
            stype = []
            strand = []
            base.append(a.split("\t")[1])
            start = int(a.split("\t")[3])
    outfile.close()   

def findDMR(name,type,chrom) :
    print "hello"
         
def main():
    #os.system("methylation_extractor -p --no_overlap "+sys.argv[1]+"_bismark_bt2_pe.sam")
    #separate_chr(sys.argv[1], sys.argv[2])
    #merge_sort_chr(sys.argv[1], sys.argv[2], sys.argv[3])
    #swDMR_input(sys.argv[1], sys.argv[2], sys.argv[3])
    
if __name__ == "__main__":
    main()

        