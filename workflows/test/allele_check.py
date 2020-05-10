#!/usr/bin/env python37

import os
import sys

in_file = sys.argv[1]
#out1 = sys.argv[2]
#out2 = sys.argv[3]
#g = open(out_file, 'w')

counts = []

with open(in_file, 'r') as f: 
    for line in f:
        line = line.strip().split("\t")
#        print(line)
        line = ' '.join(line).split()
        #print(line)
        if line[2][-1:]=='#' or line[2][-4:]=='like':
            counts.append('amb')
        elif line[3][-1:]=='#' or line[3][-4:]=='like':
            counts.append('amb')
#        elif line[0]==line[1]:
#            continue        
        elif line[0]==line[2] and line[1]==line[3]:
            #g.write(str(2) + '\n')
            counts.append(2)
        elif line[0]==line[3] and line[0]!=line[2] and line[1]!=line[3] and line[1]==line[2]:
            #g.write(str(2)+ '\n')
            counts.append(2)
        elif line[0]==line[2] and line[0]!=line[3] and line[1]!=line[3] and line[1]!=line[2]:
            #g.write(str(1)+ '\n')
            counts.append(1)
        elif line[0]!=line[2] and line[0]==line[3] and line[1]!=line[3] and line[1]!=line[2]:
            #g.write(str(1)+ '\n')
            counts.append(1)
        elif line[0]==line[2] and line[0]==line[3] and line[1]!=line[3] and line[1]!=line[2]:
            #g.write(str(1)+ '\n')
            counts.append(1)
        elif line[0]!=line[3] and line[0]!=line[2] and line[1]==line[3] and line[1]!=line[2]:
            #g.write(str(1)+ '\n')
            counts.append(1)
        elif line[0]!=line[3] and line[0]!=line[2] and line[1]!=line[3] and line[1]==line[2]:
            #g.write(str(1)+ '\n')
            counts.append(1)
        elif line[0]!=line[3] and line[0]!=line[2] and line[1]==line[3] and line[1]==line[2]:
            #g.write(str(1)+ '\n')
            counts.append(1)
        elif line[0]!=line[3] and line[0]!=line[2] and line[1]!=line[3] and line[1]!=line[2]:
            #g.write(str(0)+ '\n')
            counts.append(0)
#        elif line[0]==line[1]==line[2]==line[3]
#            counts.append(2)
        elif line[0]==line[2] and line[0]!=line[3] and line[1]==line[2] and line[1]!=line[3]:
            counts.append(1)
        elif line[0]!=line[2] and line[0]==line[3] and line[1]!=line[2] and line[1]==line[3]:
            counts.append(1)
#        elif line[0]==line[1] and line[1]!=line[2] and line[2]!=line[3]:
#            counts.append(0)
        else:
            counts.append('check')


list1 = [ elem for elem in counts if elem != 'amb']
#list1 = [ elem for elem in counts if elem != 'check']


diplotypes = len(counts)
haplotypes = 2*diplotypes
ac_haps = sum(list1)
amb_dipl = counts.count('amb')
ac_dipl = counts.count(2)
inc_haps = counts.count(1) + counts.count(0)*2
debug = counts.count('check')
#debug_loc = counts.index('check')
debug_loc = [i for i, x in enumerate(counts) if x == 'check']
inc_dipl1 = [i for i, x in enumerate(counts) if x == 0]
inc_dipl2 = [i for i, x in enumerate(counts) if x == 1]
inc_dipl = (inc_dipl1 + inc_dipl2)

print("Total diplotypes: {}".format(diplotypes))
print("Concordant diplotypes: {}".format(ac_dipl))
print("Discordant diplotypes: {}".format(diplotypes-ac_dipl))
print("Ambiguous diplotypes: {}".format(amb_dipl))
print("Total haplotypes: {}".format(haplotypes))
print("Correct haplotypes: {}".format(ac_haps)) 
print("Incorrect haplotypes: {}".format(inc_haps))

print("Debug: {}".format (debug))
print("Debug_diplotypes: {}".format(debug_loc))

print("Discordant diplotypes: {}".format(inc_dipl))

f.close()
#g.close()

#print (counts[:100])

