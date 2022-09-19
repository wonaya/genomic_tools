import os, sys
import scipy
import scikits.bootstrap as bootstrap

lines = []
for a in open(sys.argv[1], 'r') :
    lines.append(float(a.strip("\n")))
print len(lines)

print bootstrap.ci(data=lines, statfunction=scipy.mean) 
