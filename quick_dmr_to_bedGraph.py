import os,sys
## native form of DMR output
"""
a = open(sys.argv[1], 'r')
b = open(sys.argv[1]+".bedGraph", 'w')
alines = a.readlines()
for aline in alines[1:] :
    if int(aline.split("\t")[8]) >= 3 :
        b.write(aline.split("\t")[1].strip('"'))
        b.write("\t")
        b.write(aline.split("\t")[2])
        b.write("\t")
        b.write(aline.split("\t")[3])
        b.write("\t")
        b.write(aline.split("\t")[6])
        b.write("\n")
b.close()  
"""
### bisukit processed

a = open(sys.argv[1], 'r')
b = open(sys.argv[1]+".bedGraph", 'w')
alines = a.readlines()
for aline in alines[1:] :
    if int(aline.split("\t")[7]) >= 3 :
        b.write(aline.split("\t")[0])
        b.write("\t")
        b.write(aline.split("\t")[1])
        b.write("\t")
        b.write(aline.split("\t")[2])
        b.write("\t")
        b.write(aline.split("\t")[5])
        b.write("\n")
b.close()  
