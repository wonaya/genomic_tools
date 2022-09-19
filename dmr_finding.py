import os,sys

## save all coordinates by chromosome

## if more than x count of dmcs within 10bp, same degree

large_list = []
for a in open("B73_all3/CpG_B73_all3_bt202.sam_chr10.methylKit", 'r') :
    if 