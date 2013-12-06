## usage python bam_to_window_bedGraph.py 
## input : bam files of both IP and control
## output : GFF file of given window size, output file from multiBamCov, bedGraphs of IP normalized by control. Gives log2(fold changes) as its score. 
## requires : multiBamCov and samtools
## currently only runs for maize and arabidopsis
##
## written by Jawon Song
## questions to jawon 'at' tacc 'dot' utexas.edu
