import os,sys

outfile = open(sys.argv[1].split(".gtf")[0]+"_protein_coding_exon.bed", 'w')
for a in open(sys.argv[1], 'r') :
    if a.split("\t")[1] == "protein_coding" and a.split("\t")[2] == "exon" :
        outfile.write(a.split("\t")[0]+"\t")
        outfile.write(a.split("\t")[3]+"\t")
        outfile.write(a.split("\t")[4]+"\t")
        outfile.write(a.split("\t")[8].split(";")[1].split(' transcript_id "')[1].strip('"'))
        outfile.write("\n")
outfile.close()

### merge same gene 
transcripts = []
for a in open(sys.argv[1].split(".gtf")[0]+"_protein_coding_exon.bed", 'r') : 
    if a.split("\t")[3].strip("\n") not in transcripts :
        transcripts.append(a.split("\t")[3].strip("\n"))

outfile = open(sys.argv[1].split(".gtf")[0]+"_protein_coding_exon_merged_transcripts.bed", 'w')
for transcript in transcripts :
    loc_list = []
    for a in open(sys.argv[1].split(".gtf")[0]+"_protein_coding_exon.bed", 'r') : 
        if a.split("\t")[3].strip("\n") == transcript :
            loc_list.append(int(a.split("\t")[1]))
            loc_list.append(int(a.split("\t")[2]))
            chrom = a.split("\t")[0]
    outfile.write(chrom)
    outfile.write("\t")
    outfile.write(str(min(loc_list)))
    outfile.write("\t")
    outfile.write(str(max(loc_list)))
    outfile.write("\t")
    outfile.write(".\t")
    outfile.write(transcript)
    outfile.write("\n")
    transcript, min(loc_list), max(loc_list)
outfile.close() 
sys.exit()

## sample 2, 4, 7
large_count = []
outfile = open("2_7_12_count.txt", 'w')
outfile.write("area\t2C BR1\t2C BR2\t2C BR3\n")
area = 0
for a in open(sys.argv[1], 'r') :
    if a.split("\t")[0] == "10" :
        large_count.append([int(a.split("\t")[5]),int(a.split("\t")[10]),int(a.split("\t")[15])])
        if [int(a.split("\t")[5]),int(a.split("\t")[10]),int(a.split("\t")[15])].count(0) != 3 :
            area += 1
            outfile.write(str(area)+"\t"+a.split("\t")[5]+"\t"+a.split("\t")[10]+"\t"+a.split("\t")[15]+"\n")
outfile.close()
