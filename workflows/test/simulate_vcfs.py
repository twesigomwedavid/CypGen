#!/usr/bin/env python

import os 
import sys

infile = sys.argv[1]
outfile = sys.argv[2]
sample = sys.argv[3]

f = open(infile, "r")
g = open(outfile, "w")

g.write('##fileformat=VCFv4.1\n##contig=<ID=22,length=51304566>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="GenoType call. ./. is called if there is no coverage at the variant site.">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + sample + "\n")

for line in f:
    line = line.strip().split(";")
    for var in line:
        var = "22\t" +  var.replace("~","\t.\t").replace(">", "\t")
        if var[-3:] == "0/1":
            var = var.replace("0/1", ".\t.\tGT\t0/1")
            g.write(var + "\n")
        elif var[-3:] == "1/1":
            var = var.replace("1/1", ".\t.\tGT\t1/1")
            g.write(var + "\n")

f.close()
g.close()
