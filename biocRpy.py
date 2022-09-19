### Bioconductor, R, Python
### Started 16 April, 2013

import os,sys
from optparse import OptionParser
import rpy2.robjects as robjects
from rpy2.robjects import r
from rpy2.robjects.packages import importr
import csv

def main():
    ## specify R function name
    parser = OptionParser()
    parser.add_option("-p", "--prog", dest="prog",help="R function name, type --avail if not sure", type="string")
    parser.add_option("-a", "--avail", dest="avail") ## load available programs
    parser.add_option("-d", "--data", dest="data",help="Data in tab delimited text file", type="string")
    parser.add_option("--header", dest="header", help="Has header", type="string", default="yes")
    (options, args) = parser.parse_args()
    ## search for lower case of program
    #base = importr('base') ## use for library installation
    #base.source("http://www.bioconductor.org/biocLite.R")
    r.library(str(options.prog))
    #run_program(options.prog, read_in(options.header, options.data))
    run_program(options.prog, options.data)

def read_in(header, infile): ## read in txt file and convert it in appropriate way
    from rpy2.robjects.vectors import DataFrame
    f1 = open(infile, 'r') 
    f1lines = f1.readlines()
    if header == "yes" :
        data = []
        names = f1lines[0].strip().split("\t")
        for f1line in f1lines[1:] :
            data.append(f1line.strip().split("\t"))
    columns=zip(*data)
    d = {robjects.StrVector(columns[0]), robjects.IntVector(columns[1])}
    names=['col{i}'.format(i=i) for i in range(2)]
    dataf = r['data.frame'](**dict(zip(names,d)))
    return dataf

def run_program(prog, data):
    from rpy2.robjects.vectors import DataFrame
    dataf = DataFrame.from_csvfile(data, sep='\t',header=True)
    robjects.r('library('+str(prog)+')')
    try:
        robjects.r["exactTest"]
        is_13_plus = True
    except LookupError:
        is_13_plus = False
    robjects.r('group <- factor(c(1,2))')
    robjects.r('y <- DGEList(counts='+str(dataf)+',group=group)')

if __name__ == "__main__":
    main()
