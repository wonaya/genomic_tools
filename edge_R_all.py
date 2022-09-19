import csv
import collections
import os,sys

def gtf_to_bed(gtf) :
    chr = ['1','2','3','4','5','6','7','8','9','10','Mt','Pt','UNKNOWN']
    genes = []
    genes_dict = {}
    
    for a in open(gtf, 'r'): 
        if a.split("\t")[8].split(";")[0].split(' gene_id "')[1].strip('"') not in genes  and str(a.split("\t")[0]) in chr :
            genes.append(a.split("\t")[8].split(";")[0].split(' gene_id "')[1].strip('"'))
            genes_dict[a.split("\t")[8].split(";")[0].split(' gene_id "')[1].strip('"')] = []
    for a in open(gtf, 'r'): 
        if str(a.split("\t")[0]) in chr :
            genes_dict[a.split("\t")[8].split(";")[0].split(' gene_id "')[1].strip('"')].append(int(a.split("\t")[3]))
            genes_dict[a.split("\t")[8].split(";")[0].split(' gene_id "')[1].strip('"')].append(int(a.split("\t")[4]))
    outfile = open("Zea_mays.AGPv3.18.bed", 'w')
    for a in open(gtf, 'r'): 
        if a.split("\t")[8].split(";")[0].split(' gene_id "')[1].strip('"') in genes :
            outfile.write(a.split("\t")[0])
            outfile.write("\t")
            outfile.write(str(min(genes_dict[a.split("\t")[8].split(";")[0].split(' gene_id "')[1].strip('"')])))
            outfile.write("\t")
            outfile.write(str(max(genes_dict[a.split("\t")[8].split(";")[0].split(' gene_id "')[1].strip('"')])))
            outfile.write("\t")
            outfile.write(a.split("\t")[8].split(";")[0].split(' gene_id "')[1].strip('"'))
            outfile.write("\t")
            outfile.write(a.split("\t")[1])
            outfile.write("\t")
            outfile.write(a.split("\t")[2])
            outfile.write("\n")
    outfile.close() 

def remove_duplicate() :
    outlist = []
    outfile = open("Zea_mays.AGPv3.18_rmdup.bed", 'w')
    for a in open("Zea_mays.AGPv3.18.bed", 'r') :
        if " ".join(a.split("\t")[:4]) not in outlist :
            outlist.append(" ".join(a.split("\t")[:4]))
    outfile = open("Zea_mays.AGPv3.18_rmdup.bed", 'w')
    for out in outlist :
        outfile.write("\t".join(out.split(" ")))
        outfile.write("\n")
    outfile.close()
            
def convert_multibam_to_csv(list) :
    bams = []
    for a in open(list, 'r') :
        bams.append(a.strip("\n"))
        if a.strip("\n") != "Input_CGATGT.F.bam" : 
            print "2"
            os.system("/opt/apps/samtools/0.1.18/samtools reheader header "+a.strip("\n")+" > "+a.strip("\n")+".fixed")
            print "3"
            os.system("mv "+a.strip("\n")+".fixed "+a.strip("\n"))
            print "4"
            os.system("/opt/apps/samtools/0.1.18/samtools index "+a.strip("\n"))
        
    os.system("/opt/apps/bedtools/2.19.0/bin/multiBamCov -bams "+" ".join(bams)+" -bed Zea_mays.AGPv3.18_rmdup.bed > "+str(list).split(".txt")[0]+"_count.txt")

def test2(list):    
    genes = []
    for a in open("Zea_mays.AGPv3.18_rmdup.bed", 'r') :
        genes.append(a.split("\t")[3].strip("\n"))
    write_list = {}
    for a in open(str(list).split(".txt")[0]+"_count.txt", 'r')  :
        write_list[a.split("\t")[3]] = a.split("\t")[4:]
    outfile = open(str(list).split(".txt")[0]+"_edgeR_inputds.txt", 'w')
    outfile.write("Gene\tInput\tAB\tAB\tAB\n")
    for gene in genes :
        if gene in write_list.keys() :
            outfile.write(gene)
            outfile.write("\t")
            outfile.write("\t".join(write_list[gene]))
    outfile.close()

def edger_prep(list) :
    os.system("Rscript edgeR.R "+str(list).split(".txt")[0]+"_edgeR_inputds.txt "+str(list).split(".txt")[0]+"_edgeR_inputds.txt.txt")
    b = open(str(list).split(".txt")[0]+"_edgeR_inputds_result.txt" , 'w')
    a = open(str(list).split(".txt")[0]+"_edgeR_inputds.txt.txt" , 'r') 
    alines = a.readlines()
    count = 0
    sig_gene = 0
    b.write("Gene\tlogFC\tlogCPM\tLR\tp.value\tFDR\n")
    for aline in alines :
        if count != 0 :
            if float(aline.split("\t")[1]) > 1 and float(aline.split("\t")[4]) <= 0.05 :
                b.write(aline.split("\t")[0].strip('"'))
                b.write("\t")
                b.write("\t".join(aline.split("\t")[1:]))
                sig_gene += 1
        count += 1
    b.close()
    print list, sig_gene

def test(infile) :
    a = open(infile, 'r')
    alines = a.readlines()
    count = 0
    sig_gene = 0
    for aline in alines :
        if count != 0 :
            if float(aline.split("\t")[1]) > 1 and float(aline.split("\t")[4]) <= 0.05 :
                sig_gene += 1
        count += 1
    print sig_gene


def edger_post(list) :
    dict_test = {}
    a = open(str(list).split(".txt")[0]+"_edgeR_result.txt", 'r')
    alines = a.readlines()
    count = 0
    for aline in alines :
        if count != 0 :
            dict_test[aline.split("\t")[0]] = [aline.split("\t")[1],aline.split("\t")[4]]
        count += 1
    
    outfile = open("Top-98-pct-genes_edgeR_ab1.txt", 'w')
    for genes in open("Top-98-pct_rpkm.txt", 'r') :
        if genes.split("\t")[0] in dict_test.keys() :
            outfile.write(genes.split("\t")[0])
            outfile.write("\t")
            outfile.write(genes.split("\t")[1].strip("\n"))
            outfile.write("\t")
            outfile.write("\t".join(dict_test[genes.split("\t")[0]]))
            outfile.write("\n")
        else :
            outfile.write(genes.split("\t")[0])
            outfile.write("\t")
            outfile.write(genes.split("\t")[1].strip("\n"))
            outfile.write("\tNA\tNA\n")
    outfile.close()
def annotate_genes(list) :
    sig_gene = []
    for a in open(str(list).split(".txt")[0]+"_edgeR_result.txt", 'r') :
        if a.split("\t")[0] != "Gene" :
            sig_gene.append(a.split("\t")[0])
    print len(sig_gene)
    mapman_terms = []
    mapman_terms_dict = {}
    for b in open("/work/02114/wonaya/genome/annotation/MapMan_B73_annot.txt", 'r') :
        if b.split("\t")[0] in sig_gene :
            mapman_terms.append(b.split("\t")[1].strip("\n"))
            mapman_terms_dict[b.split("\t")[1].strip("\n")] = 0
    print len(mapman_terms)
    print len(mapman_terms_dict)
    for terms in mapman_terms :
        mapman_terms_dict[terms] += 1
    outfile = open(str(list).split(".txt")[0]+"_edgeR_mapman.txt", 'w')
    for terms in mapman_terms_dict :
        outfile.write(terms+"\t"+str(mapman_terms_dict[terms])+"\n")
    outfile.close()
        
annotate_genes(sys.argv[1])
#test2(sys.argv[1])
#convert_multibam_to_csv(sys.argv[1])
#test(sys.argv[1])
#print "edger_prep"    
#edger_prep(sys.argv[1])
#print "edger_post"
#edger_post(sys.argv[1])
