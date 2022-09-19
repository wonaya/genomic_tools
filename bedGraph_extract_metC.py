import os,sys

name = str(sys.argv[1]).split(".")[0]
fromint = str(sys.argv[2])
toint   = str(sys.argv[3])
outfile = name+"_"+fromint+"_"+toint+".bedGraph"
#outfile = name+"_whole.bedGraph"

outfile = open(outfile,'w')
for a in open(sys.argv[1], 'r') :
    if int(a.split("\t")[0]) >= int(sys.argv[2]) and int(a.split("\t")[0]) <= int(sys.argv[3]) :
    #if len(a.split("\t")) > 0 :
        outfile.write(a.split("\t")[1])
        outfile.write("\t")
        outfile.write(a.split("\t")[0])
        outfile.write("\t")
        outfile.write(a.split("\t")[0])
        outfile.write("\t")
        outfile.write(a.split("\t")[5])
        outfile.write("\n")
outfile.close()