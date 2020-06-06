#!/usr/bin/env python3

import os
import sys
import subprocess
from snv_def_caller import *

database = sys.argv[1]
infile = sys.argv[2]
infile_full = sys.argv[3]
infile_full_gt = sys.argv[4]
infile_spec = sys.argv[5]


snf_def_alleles = cand_snv_allele_calling(database, infile, infile_full, infile_full_gt, infile_spec)

print ("\n----------------------------\n")
print ("CypGen Test hg38")
print ("\n----------------------------\n")
print ("\nCore variants:")
print ("\n" + get_core_variants(infile))
print ("\nCandidate alleles:")
print ("\n" + snf_def_alleles[0])
print ("\nResult:")
print ("\n" + snf_def_alleles[1])
