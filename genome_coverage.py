import os,sys

### length of genome (Zea mays 3.23)
chr_len_dict = {"Pt":140384,"Mt":569630,10:149632204,9:157038028,6:169407836,8:175377492,7:176826311,5:217959525,3:232245527,2:237917468,4:242062272,1:301476924}

### how many are associated with gene
# /work/02114/wonaya/genome/Zea_mays.AGPv3.23/annotation/Zea_mays.AGPv3.23_protein_coding_exon_merged_transcripts.bed

### remove redundant lines
a = open(sys.argv[1],'r')
alines = a.readlines()
print "no. of peaks", len(set(alines))

b = open(sys.argv[2], 'r') 
blines = b.readlines()
print "no. of genes", len(set(blines))

os.system("intersectBed -wa -a "+sys.argv[1]+" -b "+sys.argv[2]+" > "+sys.argv[1].split(".bedGraph")[0]+"_"+sys.argv[2].split("annotation/")[1].split(".bed")[0]+"_whole_peak.bedGraph")

c = open(sys.argv[1].split(".bedGraph")[0]+"_"+sys.argv[2].split("annotation/")[1].split(".bed")[0]+"_whole_peak.bedGraph", 'r')
clines = c.readlines()
print "no. of peaks overlapping with genes", len(set(clines))


  

### gene coverage in genome
