import os,sys
import numpy
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages

a = open("1_S1.unamp-bwa-pe-macs1M_peaks.overlap", 'r') 
alines = a.readlines()
x1 = []
y1 = []
for aline in alines[1:] :
    if int(aline.split("\t")[3]) is not 0 and int(aline.split("\t")[4]) is not 0 :
        x1.append(log(int(aline.split("\t")[3])))
        y1.append(log(int(aline.split("\t")[4])))
(m1,b1) = polyfit(x1,y1,1)
yp1 = polyval([m1,b1],x1)

x2 = []
y2 = []
for aline in alines[1:] :
    if int(aline.split("\t")[3]) is not 0 and int(aline.split("\t")[5]) is not 0 :
        x2.append(log(int(aline.split("\t")[3])))
        y2.append(log(int(aline.split("\t")[5])))
(m2,b2) = polyfit(x2,y2,1)
yp2 = polyval([m2,b2],x2)

x3 = []
y3 = []
for aline in alines[1:] :
    if int(aline.split("\t")[3]) is not 0 and int(aline.split("\t")[6]) is not 0 :
        x3.append(log(int(aline.split("\t")[3])))
        y3.append(log(int(aline.split("\t")[6])))
(m3,b3) = polyfit(x3,y3,1)
yp3 = polyval([m3,b3],x3)

x4 = []
y4 = []
for aline in alines[1:] :
    if int(aline.split("\t")[3]) is not 0 and int(aline.split("\t")[7].strip("\n")) is not 0 :
        x4.append(log(int(aline.split("\t")[3])))
        y4.append(log(int(aline.split("\t")[7].strip("\n"))))
(m4,b4) = polyfit(x4,y4,1)
yp4 = polyval([m4,b4],x4)

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(2,2,1)
plt.xlim(0, 15)
plt.ylim(0, 15)
plt.plot(x1,yp1,color='red')
plt.scatter(x1,y1, color='blue')
#ax.set_yscale('log')
#ax.set_xscale('log')
plt.title("Sample1 vs. Control")
plt.ylabel("sample")

ax = fig.add_subplot(2,2,2)
plt.xlim(0, 15)
plt.ylim(0, 15)
plt.plot(x2,yp2,color='red')
plt.scatter(x2,y2, color='green')
plt.title("Sample2 vs. Control")

ax = fig.add_subplot(2,2,3)
plt.xlim(0, 15)
plt.ylim(0, 15)
plt.plot(x3,yp3,color='red')
plt.scatter(x3,y3, color='brown')
plt.title("Sample3 vs. Control")
plt.ylabel("sample")
plt.xlabel("control")

ax = fig.add_subplot(2,2,4)
plt.xlim(0, 15)
plt.ylim(0, 15)
plt.plot(x4,yp4,color='red')
plt.scatter(x4,y4, color='orange')
plt.title("Sample4 vs. Control")
plt.xlabel("control")

plt.savefig("samplevsControl.png")
plt.show()
