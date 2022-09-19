import os,sys
import math
import numbers
import readline, glob
def complete(text, state):
    return (glob.glob(text+'*')+[None])[state]
readline.set_completer_delims(' \t\n;')
readline.parse_and_bind("tab: complete")
readline.set_completer(complete)

# no of bioreps
bioreps = raw_input("no of bioreps : ")        
controls = raw_input("no of controls : ")

# name of bioreps
biorep_list = []
for x in range(0, int(bioreps)) :
    biorep_list.append(raw_input("name of BAM biorep no. "+str(x)+" : "))
control_list = []
for x in range(0, int(controls)) :
    control_list.append(raw_input("name of BAM control no. "+str(x+1)+" : "))

# name affix for bioreps and controls 
biorep_name = []
for biorep in biorep_list :
    biorep_name.append(biorep.split("_")[0])
control_name = []
for control in control_list :
    control_name.append(control.split("_")[0])

# genome size    
genome_size = raw_input("genome size in bp (eg. 2.3e9 for maize) :")

# macs
print "#### running macs2"
for biorep in biorep_list :
    name = biorep_name[biorep_list.index(biorep)]
    os.system("/work/02114/wonaya/software/MACS-master/bin/macs2 callpeak -t "+biorep+" -c "+control_list[0]+" -f BAM -g "+genome_size+" --nomodel --nolambda --broad -n "+str(name))
os.system("/work/02114/wonaya/software/MACS-master/bin/macs2 callpeak -t "+control_list[0]+" -c "+control_list[0]+" -f BAM -g "+genome_size+" --nomodel --nolambda --broad -n Input")

print "#### calculating overlaps"
from scipy.stats.stats import pearsonr
biorep_pair = []
for biorep1 in biorep_list :
    for biorep2 in biorep_list :
        if biorep1 != biorep2 :
            if biorep2+"-"+biorep1 not in biorep_pair :
                name1 = biorep_name[biorep_list.index(biorep1)]
                name2 = biorep_name[biorep_list.index(biorep2)]
                os.system("intersectBed -a "+str(name1)+"_peaks.broadPeak -b "+str(name2)+"_peaks.broadPeak -wa -wb > "+str(name1)+"-"+str(name2)+".overlap")
                pval1 = []
                pval2 = []
                for a in open(str(name1)+"-"+str(name2)+".overlap", 'r') :
                    pval1.append(float(a.split("\t")[7]))
                    pval2.append(float(a.split("\t")[16]))
                print str(name1)+"-"+str(name2), "PCC:", pearsonr(pval1, pval2)[0]
                biorep_pair.append(biorep1+"-"+biorep2)

for biorep1 in biorep_list :
    name1 = biorep_name[biorep_list.index(biorep1)]
    os.system("intersectBed -a "+str(name1)+"_peaks.broadPeak -b Input_peaks.broadPeak -wa -wb > "+str(name1)+"-Input.overlap")
    pval1 = []
    pval2 = []
    for a in open(str(name1)+"-Input.overlap", 'r') :
        pval1.append(float(a.split("\t")[7]))
        pval2.append(float(a.split("\t")[16]))
    print str(name1)+"-Input", "PCC:", pearsonr(pval1, pval2)[0]