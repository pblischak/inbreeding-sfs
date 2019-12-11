#!/usr/bin/env python3

import numpy as np
from sys import argv

def get_AC(info):
    info_dict = {ii.split("=")[0]:ii.split("=")[1] for ii in info.split(";")}
    return int(info_dict["AC"])

if __name__ == "__main__":
    if len(argv) < 3:
        print("Not enough arguments...\n")
        print("Usage:")
        print("  python vcf2jsfs.py <vcf-file1> <vcf-file2> <outfile>")
        exit(-1)

    vcf1,vcf2 = argv[1],argv[2]
    outfile = argv[3]
    sfs = np.zeros((51,51))
    with open(vcf1) as f1, open(vcf2) as f2:
        vcf1_dict = {int(line.split()[1]):line.split()[7] for line in f1 if line[0] != "#"}
        vcf2_dict = {int(line.split()[1]):line.split()[7] for line in f2 if line[0] != "#"}
        vcf1_keys = list(vcf1_dict)
        vcf2_keys = list(vcf2_dict)
        mutations = np.union1d(list(vcf1_dict),list(vcf2_dict))
        print("Number of mutations: {}".format(len(mutations)))
        for s in mutations:
            try:
                AC1 = get_AC(vcf1_dict[s])
            except:
                AC1 = 0
                                
            try:
                AC2 = get_AC(vcf2_dict[s])
            except:
                AC2 = 0
            
            sfs[AC1,AC2] += 1.0

    with open(outfile,'w') as f_out:
        print("51 51 unfolded", file=f_out)
        for i in range(51):
            for j in range(51):
                print(sfs[i,j], " ", sep='', end='', file=f_out)
        print("", file=f_out)
