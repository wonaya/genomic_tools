import os,sys
import rpy2.robjects as robjects
from datetime import datetime
from rpy2.robjects.packages import importr
import multiprocessing
import readline, glob
import subprocess

def serial():
    ## with methylKit
    outfile1 = open("B73_all3/CpG_B73_all3_R1_bt202_merged_FR.methylKit",'w')
    for a in open("B73_all3/CpG_B73_all3_R1_bt202_merged.methylKit", 'r') :
        if a.split("\t")[3] == "1" :
            outfile1.write("chr"+a.split("\t")[1]+"."+a.split("\t")[0]+"\t")
            outfile1.write("chr"+a.split("\t")[1]+"\t")
            outfile1.write(a.split("\t")[2]+"\t")
            outfile1.write("\tF\t")
            outfile1.write("\t".join(a.split("\t")[4:]))
        elif a.split("\t")[3] == "0" :
            outfile1.write("chr"+a.split("\t")[1]+"."+a.split("\t")[0]+"\t")
            outfile1.write("chr"+a.split("\t")[1]+"\t")
            outfile1.write(a.split("\t")[2]+"\t")
            outfile1.write("\tR\t")
            outfile1.write("\t".join(a.split("\t")[4:]))
    outfile1.close()        
    """
    outfile2 = open("Mo17_all3/CpG_Mo17_all3_bt202_merged_FR.methylKit",'w')
    for a in open("Mo17_all3/CpG_Mo17_all3_bt202_merged.methylKit", 'r') :
        if a.split("\t")[3] == "1" :
            outfile2.write("chr"+a.split("\t")[1]+"."+a.split("\t")[0]+"\t")
            outfile2.write("chr"+a.split("\t")[1]+"\t")
            outfile2.write(a.split("\t")[2]+"\t")
            outfile2.write("\tF\t")
            outfile2.write("\t".join(a.split("\t")[4:]))
        elif a.split("\t")[3] == "0" :
            outfile2.write("chr"+a.split("\t")[1]+"."+a.split("\t")[0]+"\t")
            outfile2.write("chr"+a.split("\t")[1]+"\t")
            outfile2.write(a.split("\t")[2]+"\t")
            outfile2.write("\tR\t")
            outfile2.write("\t".join(a.split("\t")[4:]))
    outfile2.close()        
    """
    print "starting methylkit", datetime.now()
    importr('methylKit')
    robjects.r('file.list=list("B73_all3/CpG_B73_all3_R1_bt202_merged_FR.methylKit","Mo17_all3/CpG_Mo17_all3_bt202_merged_FR.methylKit")')
    robjects.r('myobj=read(file.list, sample.id=list("B73","Mo17"),assembly="B73",treatment=c(0,1),context="CpG")')
    print "done read file", datetime.now()
    robjects.r('meth=unite(myobj, destrand=FALSE)')
    print "unite done", datetime.now()
    robjects.r('''source("/work/02114/wonaya/scripts/diffMeth.R")''')
    robjects.r('''myDiff=calculateDiffMeth(meth)''')
    print "diffmeth calc done", datetime.now()
    robjects.r('''source("/work/02114/wonaya/scripts/eDMR_test.R")''')
    robjects.r('''myDMR=eDMR(myDiff,step=100,dist="none",DMC.qvalue=1,DMC.methdiff=25,num.DMCs=1,num.CpGs=3,DMR.methdiff=20,plot=FALSE,main="example",mode=1)''')
    print "edmr done", datetime.now()
    robjects.r('''myDMR.sig=filter.dmr(myDMR, mean.meth.diff=20, num.CpGs=3, num.DMCs=1)''')
    robjects.r('''myDMR.sig.df <- as.data.frame(myDMR.sig)''')
    print "sigdf dmr done", datetime.now()
    robjects.r('''write.table(myDMR.sig.df, file="speed_test.out", sep="\t")''')
    print "all done", datetime.now()

## with bigmemory
def parallel(chr):
    importr('GenomicRanges')
    importr('bigmemory')
    importr('data.table')
    robjects.r('''source('/work/02114/wonaya/scripts/base.r')''')
    robjects.r.assign('file1',"B73_all3/CpG_B73_all3_R1_bt202_chr"+str(chr)+".methylKit")
    robjects.r.assign('file2',"Mo17_all3/CpG_Mo17_all3_bt202_chr"+str(chr)+".methylKit")
    robjects.r.assign('outfile', "speed_test"+str(chr)+".out")
    robjects.r('file.list=list(file1,file2)')
    robjects.r('myobj=read(file.list, sample.id=list("B73","Mo17"), assembly="B73",treatment=c(0,1),context="CpG", resolution="base")')
    robjects.r('meth=unite(myobj, destrand=FALSE)')
    robjects.r('source("/work/02114/wonaya/scripts/diffMeth.R")')
    robjects.r('''myDiff=calculateDiffMeth(meth)''')
    robjects.r('''source("/work/02114/wonaya/scripts/eDMR_test.R")''')
    robjects.r('''myDMR=eDMR(myDiff,step=100,dist="none",DMC.qvalue=1,DMC.methdiff=25,num.DMCs=1,num.CpGs=3,DMR.methdiff=20,plot=FALSE,main="example",mode=1)''')
    robjects.r('''myDMR.sig=filter.dmr(myDMR, mean.meth.diff=20, num.CpGs=3, num.DMCs=1)''')
    robjects.r('''myDMR.sig.df <- as.data.frame(myDMR.sig)''')
    robjects.r('''write.table(myDMR.sig.df, file=outfile, sep="\t")''')

if sys.argv[1] == "serial" :
    serial()
elif sys.argv[1] == "parallel" :
    chr_list = {'1':[1],'2':[2],'3':[3],'4':[4],'5':[5],'6':[6],'7':[7],'8':[8],'9':[9],'10':[10]} 
    jobs = []
    print datetime.now()
    for chr in chr_list.keys() :
        s = multiprocessing.Process(target=parallel, args=(chr,  ))
        jobs.append(s)
        s.start()
    [x.join() for x in jobs]
    print datetime.now()