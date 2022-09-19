import os,sys

genome_list = []
chr_list = []
outfile = open("n_coord.bed", 'w')
for a in open(sys.argv[1], 'r') :
    if a[0] == ">" :
        if len(genome_list) == 0 :
            chr_list.append(a[1:3].strip(" "))
        else :
            #print chr_list[-1], genome_list ; sys.exit()
            index = 1
            n_list = []
            for b in genome_list : 
                if b == 'N' : 
                    n_list.append(index) ; print chr_list, b, index ; sys.exit()
                index += 1
            if len(n_list) > 0 : 
                for n in n_list :
                    outfile.write(chr_list[-1])
                    outfile.write("\t")
                    outfile.write(str(n))
                    outfile.write("\t")
                    outfile.write(str(n))
                    outfile.write("\n")
            n_list = []
            genome_list = []
            chr_list.append(a[1:3].strip(" "))
    else :
        genome_list.extend(a.strip("\n"))
        

