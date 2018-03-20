# sortBed -i Input_8C_allmerge_RPGCnorm_1kb.bedgraph > Input_8C_allmerge_RPGCnorm_1kb_sort.bedgraph
# intersectBed -wa -wb -a Input_8C_allmerge_RPGCnorm_1kb_sort.bedgraph -b maize_4.33_w1000.bed > Input_8C_allmerge_RPGCnorm_1kb_merge.bedgraph
# python merge_two_beds.py Input_2C_allmerge_RPGCnorm_1kb_merge.bedgraph Input_2C_allmerge_RPGCnorm_1kb_fix.bedgraph
# python similarity.py Input_2C_allmerge_RPGCnorm_1kb_fix.bedgraph Input_8C_allmerge_RPGCnorm_1kb_fix.bedgraph
import os,sys
import difflib
import numpy as np
list1 = []
for a in open(sys.argv[1], 'r'):
    list1.append(float(a.split("\t")[-1].strip("\n")))
list2 = []
for b in open(sys.argv[2], 'r'):
    list2.append(float(b.split("\t")[-1].strip("\n")))

print len(list1), len(list2)
#sm=difflib.SequenceMatcher(None, list1, list2)
#print sm.ratio()
print np.corrcoef(list1, list2)


