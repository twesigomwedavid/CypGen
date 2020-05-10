#!/bin/bash

import os 
import sys

infile = sys.argv[1]

f = open(infile, "r")

for line in f:
    line = line.strip().split("_")
    res = '_'.join(sorted(line))
    print(str(res))

f.close()
