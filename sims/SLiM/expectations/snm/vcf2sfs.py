#!/usr/bin/env python3

import numpy as np
from sys import argv

if __name__ == "__main__":
    if len(argv) < 3:
        print("Not enough arguments...\n")
        print("Usage:")
        print("  python vcf2sfs.py <vcf-file> <nind>")
        exit(-1)

    vcf = argv[1]
    outfile = ".".join(vcf.split(".")[0:-1]) + ".fs"
    nind = int(argv[2])
    sfs = np.zeros((2*nind+1,))
    with open(vcf) as f:
        for line in f:
            if line[0] == "#":
                continue
            else:
                info = line.split()[7]
                info_dict = {ii.split("=")[0]:ii.split("=")[1] for ii in info.split(";")}
                AC = int(info_dict["AC"])
                sfs[AC] += 1.0
    
    
    with open(outfile,'w') as f_out:
        print("{} unfolded".format(2*nind+1), file=f_out)
        for i in range(2*nind):
            print(sfs[i], " ", sep='', end='', file=f_out)
        print(sfs[-1], file=f_out)
