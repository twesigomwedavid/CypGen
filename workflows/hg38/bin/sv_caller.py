#!/usr/bin/env python3

import os
import sys
import math
#import statistics as st


#sv_del = sys.argv[1]
#sv_dup = sys.argv[2]
#cov_file = sys.argv[3]


# def custom_round(m):
#     if "." not in str(m):
        

def get_total_CN(cov_file):

    all_reg =[]
    for line in open(cov_file, "r"):
        line = line.strip().split()
        all_reg.append(line)

    av_2d6_cov = float(all_reg[2][3])/(float(all_reg[2][2]) - float(all_reg[2][1]))
    av_2d8_cov = float(all_reg[3][3])/(float(all_reg[3][2]) - float(all_reg[3][1]))
    av_in1_3pr = float(all_reg[1][3])/(float(all_reg[1][2]) - float(all_reg[1][1]))
    av_ex9_3pr = float(all_reg[0][3])/(float(all_reg[0][2]) - float(all_reg[0][1]))
   # av_2d6_cov = float(all_reg[0[3]])/(float(all_reg[0[2]]) - float(all_reg[0[1]]))
   # av_2d8_cov = float(all_reg[1[3]])/(float(all_reg[1[2]]) - float(all_reg[1[1]]))
   # 2d6_covs.append(int(line[-1]))
   # av_2d6_cov = st.mean(2d6_covs)    
   # av_2d8_cov = st.mean(2d8_covs)

    comp_av = av_2d6_cov/av_2d8_cov
    temp_cn = 2 * comp_av
    total_cn = round(temp_cn)
    
    in1_3pr = round(2 * av_in1_3pr/av_2d8_cov) 
    ex9_3pr = (2 * av_ex9_3pr/av_2d8_cov)


    return [str(int(total_cn)), round(av_2d6_cov), str(int(in1_3pr)), round(av_2d8_cov), str(ex9_3pr)];

    # if "." not in test_av:
    #     cn_frac = test_av
    #     return  

samp_gt = ""
samp_gt_hap1 = ""

def del_test(sv_del):

    if os.stat(sv_del).st_size == 0:
        return "None"

    else:
        for line in open(sv_del, "r"):
            if "COVERAGE" in line:
                line = line.strip().split()
                #print (line)
                ABHom = line[-1]
                ABHet = line[-2]
                GT = line[2]
                DP = int(line[3])
        
                if float(ABHom) == 1.0:
                    return "*5/*5"
                elif float(ABHom) == -1.0:
                    return "*5"
            else:
                pass

het_hom_list = []
het_hom_list_new = []

def dup_test_init(sv_dup, av_cov):
    for line in open(sv_dup, "r"):
        if "COVERAGE" in line:
            continue
        elif "AGGREGATED" in line:
            continue

        else:
            fields = line.strip().split()
            het_hom_list.append(fields)

    test_list1 = []

    for i in het_hom_list:
        test_list1.append(int(i[2]))

    av_read_cov = sum(test_list1)/len(test_list1)
    norm_cov = (av_cov + av_read_cov)/2

    for i in het_hom_list:
        supp_reads = round(float(i[-2])*int(i[2]))
        i.append(round(supp_reads/norm_cov, 3))
        i.append(supp_reads)
        het_hom_list_new.append(i)


    return (het_hom_list_new)



hap_def_list = []
allele_cn_list = []

def dup_test_cn_3_4(sv_dup, hap_dbs, cand_allele1, cand_allele2, test_allele1, test_allele2, c_num, av_cov, in_list):

    #cn_amp = get_total_CN(bedcov_file)
    g = open(hap_dbs, "r")
    for line in g:
        line = line.strip().split()
        hap_def_list.append(line)
        
            
    test_list1 = []
    test_list2 = []
    het_list = []

    # return in_list

    for i in in_list:
        if i[1] == "0/1":
            het_list.append(i)

    for i in het_list:
        test_list1.append(i[0])
        test_list2.append(i[-2])

    # var = het_list[-1][0]
    # ABHom = het_list[-1][-1]
    # ABHet = het_list[-1][-2]
    # #GT = line[1]
    # #DP = line[2]
    # #return [var, test_allele1, test_allele2];

    max_het = max(test_list2)
    max_het_pos = test_list2.index(max_het)
    var = test_list1[max_het_pos]
    
    for elem in hap_def_list:
        if elem[1] == cand_allele1:
    #     #if line.startswith(test_allele1):
            list_3t = elem
            list_3t_2 = list_3t[2].split(';')
            l3 = len(list_3t_2)
            
        #if line.startswith(test_allele2):
        if elem[1] == cand_allele2:
            list_4t = elem
            list_4t_2 = list_4t[2].split(';')
            l4 = len(list_4t_2)

    hdb_list = list_3t_2 + list_4t_2

    #return hdb_list

    # for i in hdb_list:
    #     if i == test_list1[max_het_pos] and hdb_list.index(i) < l3:
    #     # if i in test_list1 and hdb_list.index(i) < l3:
    #     # var = i
    #     #     ind1 = test_list1.index(var1)
    #     #     abhet1 = het_list[ind1][-2]
    #     #     abhet1 = max_het
    #         break

    #     elif i == test_list1[max_het_pos] and hdb_list.index(i) >= l3:
    #         # var = i
    #         # ind2 = test_list1.index(var2)
    #         # abhet2 = het_list[ind2][-2]
    #         # abhet2 = max_het
    #         break

    # return [var, max_het];

    index_var = hdb_list.index(var)

    if index_var < l3:
        allele_cn_list.append(test_allele1)
        allele_cn_list.append(int(round(max_het*int(c_num))))

    elif index_var >= l3:
        allele_cn_list.append(test_allele2)
        allele_cn_list.append(int(round(max_het*int(c_num))))

    #return allele_cn_list
    
    if allele_cn_list[0] == test_allele1:
        rt_2 = int(c_num) - allele_cn_list[1]
        allele_cn_list.append(test_allele2)
        allele_cn_list.append(rt_2)

    elif allele_cn_list[0] == test_allele2:
        rt_2 = int(c_num) - allele_cn_list[1]
        allele_cn_list.append(test_allele1)
        allele_cn_list.append(rt_2)

    # return [allele_cn_list, av_cov, max_het, var];

    if allele_cn_list[1] == 1:
        res_dip = allele_cn_list[0] + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3])

    elif allele_cn_list[3] == 1:
        res_dip = allele_cn_list[2] + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1])  

    elif allele_cn_list[1] == 0:
        res_dip = allele_cn_list[0] + "/" + allele_cn_list[2] + "x" + str(allele_cn_list[3] - 1)

    elif allele_cn_list[3] == 0:
        res_dip = allele_cn_list[2] + "/" + allele_cn_list[0] + "x" + str(allele_cn_list[1] - 1)

    else:
        res_dip = 'check'

    return res_dip


def hybrid_test_68(sv_dup, c_num, av_cov, cn_in1_3pr1, in_list):
    
    test_list1 = []
    test_list2 = []
    test_list3 = []

    # return in_list

    for i in in_list:
        test_list1.append(i[0])
        test_list2.append(abs(float(i[-2])))
        test_list3.append(i[-1])

    index1 = test_list1.index('42526694~G>A')
    index2 = test_list1.index('42524947~C>T')
    
    val_68 = test_list3[index1]
    val_4 = test_list3[index2]
    
    #return [val_68, val_4];
    #cn_4 = int(round(val_4 * int(c_num)))
    #cn_68 = int(round(val_68 * int(c_num)))
    rt = round(val_68/val_4)

    #return rt

    # if cn_4 == cn_68:
    #     return 'norm_dup'
    # elif cn_4 < cn_68:
    #     return 'hyb_68'
    
    if rt <= 1.5:
        return 'norm_dup'

    elif rt > 1.5:
        # if int(cn_in1_3pr) == 2: 
        return 'hyb_68'

        # elif int(cn_in1_3pr) == int(c_num):
        #     return 'hyb_68_dup'

        # elif int(cn_in1_3pr) == 1:
        #     return 'hyb_68_del'
       
    else:
        return 'norm_dup'


def hyb_test_5_68_4(sv_del):
    test_del = []
    for line in open(sv_del, "r"):
        if "COVERAGE" in line:
            test_del.append(line.strip())

    if len(test_del) == 0:
        return 'norm_art'

    elif len(test_del) > 0:
        return 'del_hyb'

def hybrid_test_36(sv_dup, cn, av_cov, cn_ex9_3pr):

    if int(round(float(cn_ex9_3pr))) == int(cn):
        return 'norm_dup'
    # else:
    #     return cn_ex9_3pr

    elif ((int(cn) - 1) - 0.3) < float(cn_ex9_3pr) < ((int(cn) - 1) + 0.5):
        return 'hyb_36_10'
    elif (int(cn) - 2) <= float(cn_ex9_3pr) < (int(cn) - 2 + 0.7):
        return 'hyb_36_36'


def hybrid_test_36_mod(sv_dup, cn, av_cov, cn_ex9_3pr):

    if int(round(float(cn_ex9_3pr))) == int(cn):
        return 'norm_mt'
                                                                                                  
    elif ((int(cn) - 1) - 0.3) < float(cn_ex9_3pr) < ((int(cn) - 1) + 0.5):
        return 'hyb_36_10'
    elif (int(cn) - 2) <= float(cn_ex9_3pr) < (int(cn) - 2 + 0.7):
        return 'hyb_36_36'
