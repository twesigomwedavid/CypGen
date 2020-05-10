#!/usr/bin/env python3

import os
import sys


print("--------------------------------------------\n")

print("CYP2D6 Star Allele Calling with CypGen\n")

print("--------------------------------------------\n")

print("Candidate SNV-defined alleles in sample:\n")


database = sys.argv[1]
infile = sys.argv[2]
infile_full = sys.argv[3]
#infile_full_gt = sys.argv[4]
#infile_spec = sys.argv[5]


#f = open(infile_spec, "r")


for line in open(infile_full, "r"):
    all_variants = line.strip().split(";")
#print(all_variants)


if os.stat(infile).st_size == 0:
    print("\n*2/*2")
    print("\nSupporting variants")
    print("\n" + "".join(all_variants))
    sys.exit()

core_variants = []

for line in open(infile, "r"):
    line = line.strip()
    core_variants.append(line)

core_variants = ";".join(sorted(core_variants))

# all_var_gt = []
# for line in open(infile_full_gt, "r"):
#     line = line.strip()
#     all_var_gt.append(line)

#print(core_variants)
    
dbs = []

for line in open(database, "r"):
    line = line.strip().split("\t")
    dbs.append(line)

soln_list1 = []
soln_list2 = []

for record in dbs:
    record_core_var = record[1].split(";")
    record_core_var = ";".join(sorted(record_core_var))
    if record_core_var == core_variants:
        #diplo = record[0]
        #full_dip = record[2]
        soln_list1.append(record[0])
        soln_list2.append(record[2])
    else:
        pass

print(soln_list1)

if len(soln_list1) == 1:
    diplo = "".join(soln_list1)
    res1 = [i for i in range(len(diplo)) if diplo.startswith("_", i)]
    res2 = [i for i in range(len(diplo)) if diplo.startswith(".", i)]
    hap1 = "*" + str (diplo[:res2[0]])
    hap2 = "*" + str (diplo[res1[0]+1:res2[1]])
    print("\n" + hap1 + "/" + hap2 + "\n")
    print ("\nSupporting variants:")
    print ("\n" + core_variants + "\n")

elif len(soln_list1) == 2:
    diplo1 = soln_list1[0]
    diplo2 = soln_list1[1]
    diplo1_supp_var = soln_list2[0].split(";")
#    print (diplo1_supp_var)
    diplo2_supp_var = soln_list2[1].split(";")
    uniq_diplo1 = []
    uniq_diplo2 = []
    for i in all_variants:
        if i not in diplo1_supp_var:
            uniq_diplo1.append(i)
        
        if i not in diplo2_supp_var:
            uniq_diplo2.append(i)

    print("\nUnique variants in soln 1: {}".format(len(uniq_diplo1)))
    print("\nUnique variants in soln 2: {}".format(len(uniq_diplo2)))
            
    if len(uniq_diplo1) < len(uniq_diplo2):
        res1 = [i for i in range(len(diplo1)) if diplo1.startswith("_", i)]
        res2 = [i for i in range(len(diplo1)) if diplo1.startswith(".", i)]
        hap1 = "*" + str (diplo1[:res2[0]])
        hap2 = "*" + str (diplo1[res1[0]+1:res2[1]])
        print("\n" + hap1 + "/" + hap2 + "\n")
        print ("Supporting variants:")
        print ("\n" + core_variants + "\n")

    elif len(uniq_diplo1) > len(uniq_diplo2):
        res1 = [i for i in range(len(diplo2)) if diplo2.startswith("_", i)]
        res2 = [i for i in range(len(diplo2)) if diplo2.startswith(".", i)]
        hap1 = "*" + str (diplo2[:res2[0]])
        hap2 = "*" + str (diplo2[res1[0]+1:res2[1]])
        print("\n" + hap1 + "/" + hap2 + "\n")
        print ("Supporting variants:")
        print ("\n" + core_variants + "\n")

    # else:
    #     tiebreak1 = []
    #     tiebreak2 = []
    #     score = []
    #     for line in f:
    #         line = line.strip().split()
    #         #print(line)
    #         if line[1] == core_variants:
    #             tiebreak1.append(line[0])
    #             tiebreak2.append(line[2])
    #     for full_dip in tiebreak2:
    #         diplo_supp_gt = full_dip.split(";")
    #         uniq_gt = []
    #         for i in all_var_gt:
    #             if i not in diplo_supp_gt:
    #                 uniq_gt.append(i)
    #         score_dip = len(uniq_gt)
    #         score.append(score_dip)

    #     minpos = score.index(min(score))
    #     best_diplo = tiebreak1[minpos]
    #     res1 = [i for i in range(len(best_diplo)) if best_diplo.startswith("_", i)]
    #     res2 = [i for i in range(len(best_diplo)) if best_diplo.startswith(".", i)]
    #     hap1 = "*" + str (best_diplo[:res2[0]])
    #     hap2 = "*" + str (best_diplo[res1[0]+1:res2[1]])
    #     print("\n" + hap1 + "/" + hap2 + "\n")
    #     print ("Supporting core variants:")
    #     print ("\n" + core_variants + "\n")


elif len(soln_list1) == 3:
    diplo1 = soln_list1[0]
    diplo2 = soln_list1[1]
    diplo3 = soln_list1[2]
    diplo1_supp_var = soln_list2[0].split(";") 
    diplo2_supp_var = soln_list2[1].split(";")
    diplo3_supp_var = soln_list2[2].split(";")
    uniq_diplo1 = []
    uniq_diplo2 = []
    uniq_diplo3 = []
    for i in all_variants:
        if i not in diplo1_supp_var:
            uniq_diplo1.append(i)

        if i not in diplo2_supp_var:
            uniq_diplo2.append(i)

        if i not in diplo3_supp_var:
            uniq_diplo3.append(i)


    print("\nUnique variants in soln 1: {}".format(len(uniq_diplo1)))
    print("\nUnique variants in soln 2: {}".format(len(uniq_diplo2)))
    print("\nUnique variants in soln 3: {}".format(len(uniq_diplo3)))


    if len(uniq_diplo1) < len(uniq_diplo2) and len(uniq_diplo1) < len(uniq_diplo3):
        res1 = [i for i in range(len(diplo1)) if diplo1.startswith("_", i)]
        res2 = [i for i in range(len(diplo1)) if diplo1.startswith(".", i)]
        hap1 = "*" + str (diplo1[:res2[0]])
        hap2 = "*" + str (diplo1[res1[0]+1:res2[1]])
        print("\n" + hap1 + "/" + hap2 + "\n")
        print ("Supporting variants:")
        print ("\n" + core_variants + "\n")

    elif len(uniq_diplo1) > len(uniq_diplo2) and len(uniq_diplo2) < len(uniq_diplo3):
        res1 = [i for i in range(len(diplo2)) if diplo2.startswith("_", i)]
        res2 = [i for i in range(len(diplo2)) if diplo2.startswith(".", i)]
        hap1 = "*" + str (diplo2[:res2[0]])
        hap2 = "*" + str (diplo2[res1[0]+1:res2[1]])
        print("\n" + hap1 + "/" + hap2 + "\n")
        print ("Supporting variants:")
        print ("\n" + core_variants + "\n")

    elif len(uniq_diplo1) > len(uniq_diplo2) and len(uniq_diplo2) > len(uniq_diplo3):
        res1 = [i for i in range(len(diplo2)) if diplo3.startswith("_", i)]
        res2 = [i for i in range(len(diplo2)) if diplo3.startswith(".", i)]
        hap1 = "*" + str (diplo3[:res2[0]])
        hap2 = "*" + str (diplo3[res1[0]+1:res2[1]])
        print("\n" + hap1 + "/" + hap2 + "\n")
        print ("Supporting variants:")
        print ("\n" + core_variants + "\n")

    else:
        res1 = [i for i in range(len(diplo2)) if diplo3.startswith("_", i)]
        res2 = [i for i in range(len(diplo2)) if diplo3.startswith(".", i)]
        hap1 = "*" + str (diplo3[:res2[0]])
        hap2 = "*" + str (diplo3[res1[0]+1:res2[1]])
        print("\n" + hap1 + "/" + hap2 + "\n")
        print ("Supporting variants:")
        print ("\n" + core_variants + "\n")
        
        
