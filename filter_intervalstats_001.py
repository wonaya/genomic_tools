import os,sys

file_q = sys.argv[1]
file_r = sys.argv[2]
outfile = sys.argv[3]

os.system("/work/02114/wonaya/software/IntervalStats/IntervalStats -q "+str(file_q)+" -r "+str(file_r)+" -d /work/02114/wonaya/genome/annotation/Zea_mays.AGPv2.14_protein_coding_exon.bed -o "+str(outfile))

"""
1. Query interval
2. Closest reference interval
3. Length of query
4. Distance
5. Numerator
6. Denominator
7. p-value (quotient of above two)
"""
outfile2 = open(str(sys.argv)[3]+"_p001.txt", 'w')
count_total_q = 0
count_sig_peak = 0
for a in open(outfile, 'r') :
    count_total_q += 1
    if float(a.split("\t")[6].strip("\n")) <= 0.01 :
        outfile2.write(a)
        count_sig_peak += 1

outfile2.close()
print "total peak in GFF:", count_total_q, "overlap peak:", count_sig_peak
    