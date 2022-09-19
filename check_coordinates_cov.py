import os,sys

for a in open(sys.argv[1], 'r') :
    if int(a.split("\t")[1]) >= int(sys.argv[2]) and int(a.split("\t")[1]) <= int(sys.argv[3]) :
        print a
