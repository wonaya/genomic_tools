import os,sys


## merge htseq-count outfiles
dict_count = {}
a00 = open("1_S1-bwa-pe.1200k.sorted.htseq.out")
a00lines = a00.readlines()
for a00line in a00lines[:-5] :
    dict_count[a00line.split("\t")[0]] = []
order = [1,6,2,3,4,5,7,8,9,10]
total = ["Total"]
for x in order :
    Total = 0
    a0 = open(str(x)+"_S"+str(x)+"-bwa-pe.1200k.sorted.htseq.out")
    a0lines = a0.readlines()
    for a0line in a0lines[:-5] :
        Total += int(a0line.split("\t")[1].strip("\n"))
        dict_count[a0line.split("\t")[0]].append(a0line.split("\t")[1].strip("\n"))
    total.append(str(Total))
outfile = open("htseq_control_vs_light_miseq.csv", 'w')
#outfile.write("gene,c1,c2,l1,l2,l3,l4,d1,d2,d3,d4\n")
outfile.write("gene,c1,c2,l1,l2,l3,l4,d1,d2,d3,d4\n")
for gene in dict_count.keys():
    outfile.write(str(gene)+","+",".join(dict_count[gene]))
    outfile.write("\n")
outfile.write(",".join(total))
outfile.write("\n")
outfile.close()
sys.exit()

"""
chrlist = []
for x in range(1,11) :
    a0 = open(str(x)+"_S"+str(x)+"-bwa-pe.1200k.wig", 'r') 
    a0lines = a0.readlines()
    for a0line_chr in a0lines :
        if a0line_chr[:3] == "var" and a0line_chr.split("chrom=")[1].split(" span=")[0] not in chrlist:
            chrlist.append(a0line_chr.split("chrom=")[1].split(" span=")[0])
print chrlist
max_result = []
for chr in chrlist :
    max_by_chr = []
    for x in range(1,11) :
        a0 = open(str(x)+"_S"+str(x)+"-bwa-pe.1200k.wig", 'r') 
        a0lines = a0.readlines()
        ends = []
        for a0line_test in a0lines :
            if chr != chrlist[-1] :
                if a0line_test[:3] == "var" and a0line_test.split("chrom=")[1].split(" span=")[0] == chrlist[int(chrlist.index(chr)+1)] :
                    ends.append(int(a0lines.index(a0line_test)))
            elif chr == chrlist[-1] :
                ends.append(int(len(a0lines)))
                end = ends[0]
        if chr != chrlist[-1] :
            end =  ends[0]-1
        elif chr == chrlist[-1] :
            end = ends[0]
        max_by_chr.append(int(a0lines[int(end)-1].split("\t")[0]))
    max_result.append(max(max_by_chr))
print max_result

x = sys.argv[1]
a0 = open(str(x)+"_S"+str(x)+"-bwa-pe.1200k.wig", 'r') 
a0lines = a0.readlines()
chrlist = []
for a0line_chr in a0lines :
    if a0line_chr[:3] == "var" and a0line_chr.split("chrom=")[1].split(" span=")[0] not in chrlist:
        chrlist.append(a0line_chr.split("chrom=")[1].split(" span=")[0])
result = []
for chr in chrlist :
    result_by_chr = {}
    #print "chr", chr
    ends = []
    for a0line_test in a0lines :
        if chr != chrlist[-1] :
            if a0line_test[:3] == "var" and a0line_test.split("chrom=")[1].split(" span=")[0] == chrlist[int(chrlist.index(chr)+1)] :
                ends.append(int(a0lines.index(a0line_test)))
        elif chr == chrlist[-1] :
            ends.append(int(len(a0lines)))
            end = ends[0]
    if chr != chrlist[-1] :
        end =  ends[0]-1
    elif chr == chrlist[-1] :
        end = ends[0]
    starts =[]
    for a0line_test in a0lines :
        if a0line_test[:3] == "var" and a0line_test.split("chrom=")[1].split(" span=")[0] == chrlist[int(chrlist.index(chr))] :
            starts.append(int(a0lines.index(a0line_test)))
    start =starts[0]+1
    for result_line in a0lines[start:end] :
        if not result_line[0] == "#" :
            if not result_line[:3] == "var" : 
                result_by_chr[int(result_line.split("\t")[0])] = result_line.split("\t")[1].strip("\n")
    result.append(result_by_chr)

outfile = open(str(x)+"_wig_sum.txt", 'w')
for z in range(0, len(chrlist)) :
    y = 1
    while y <= max_result[0] :
        if y in (result[0]).keys() :
            outfile.write(str(chrlist[z])+"\t"+str(y)+"\t"+str(result[z][y])+"\n") ### first file, chr
        else :
            outfile.write(str(chrlist[z])+"\t"+str(y)+"\t0\n")
        y += 100
outfile.close()

large_list = []
for x in range(0, len(chrlist)) :
    y = 1
    while y <= max_result[x] :
        small_list = []
        small_list.append(y)
        for whole in whole_result :
            test_list = []
            for location in whole[x] :
                if int(location.split("\t")[0]) == y :
                    #small_list.append("file"+whole_result.index(whole)+1+":"+float(location.split("\t")[1].strip("\n")))
                    test_list.append(float(location.split("\t")[1].strip("\n")))
            if len(test_list) == 0 :
                small_list.append(0)
            else :
                small_list.append(test_list[0])
        print small_list
        y += 100
    sys.exit() 
"""   
sig_gene = []
a = open(sys.argv[1],'r') 
alines = a.readlines()
for aline in alines[1:] :
    if float(aline.split(",")[3]) < 0.05 and aline.split(",")[0].strip('"') not in sig_gene and float(aline.split(",")[1]) > 0 : 
        sig_gene.append(aline.split(",")[0].strip('"'))

annotation_dict = {}
for gene in sig_gene :
    annotation_dict[gene] = []
evid_code = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP']
process_code = ['P']
annotation = []
annotation_desc = {}
for b in open("/work/02114/wonaya/genome/annotation/ATH_GO_GOSLIM.txt", 'r') :
    if b.split("\t")[0] in annotation_dict.keys() and b.split("\t")[7] in process_code and b.split("\t")[9] in evid_code :
        annotation_dict[b.split("\t")[0]].append(",".join([b.split("\t")[5], b.split("\t")[3], b.split("\t")[4]]))
        annotation.append(b.split("\t")[5])
        annotation_desc[b.split("\t")[5]] = str(b.split("\t")[3])+"\t"+str(b.split("\t")[4])
        
annotation_count = {}
annotation_count_order = []
for annotate in annotation :
    annotation_count[annotation.count(annotate)] = [] 
    if int(annotation.count(annotate)) not in annotation_count_order :
        annotation_count_order.append(int(annotation.count(annotate)))
annotation_count_order.sort()
annotation_count_order.reverse()
for annotate in annotation :
    if annotate not in annotation_count[annotation.count(annotate)] :
        annotation_count[annotation.count(annotate)].append(annotate)

outfile = open(str(sys.argv[1]).strip(".csv")+"_GOterms_p005_P_upreginD.out", 'w')
for annot in annotation_count_order: 
    for goterm in annotation_count[annot] :
        outfile.write(str(goterm)+"\t"+str(annot)+"\t"+str(annotation_desc[goterm])+"\n")
outfile.close()
    