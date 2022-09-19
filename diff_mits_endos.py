import os,sys

### overlap between mitotic and endocycle 

### early S

os.system("intersectBed -a /scratch/02114/wonaya/NCSU_HiSeq/7-21-14_RelicationTiming_Hiseq/LS_sig_region.bedGraph -b /scratch/02114/wonaya/NCSU_HiSeq/09-02-14_Maize_endoreplicationTiming_Hiseq/Endo_LS_sig_region.2.bedGraph > LS_mitS_endoS_overlap.bedGraph")

### find length
mits = 0
for a in open("/scratch/02114/wonaya/NCSU_HiSeq/7-21-14_RelicationTiming_Hiseq/LS_sig_region.bedGraph", 'r') :
    mits += int(a.split("\t")[2])-int(a.split("\t")[1])

endos = 0
for a in open("/scratch/02114/wonaya/NCSU_HiSeq/09-02-14_Maize_endoreplicationTiming_Hiseq/Endo_LS_sig_region.2.bedGraph", 'r') :
    endos += int(a.split("\t")[2])-int(a.split("\t")[1])

overlaps = 0
for a in open("LS_mitS_endoS_overlap.bedGraph", 'r') :
    overlaps += int(a.split("\t")[2])-int(a.split("\t")[1])

print overlaps, mits, endos
sys.exit()
#print float(overlaps)/float(mits), float(overlaps)/float(endos)
