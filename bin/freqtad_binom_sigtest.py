import numpy as np
import itertools as it
import scipy
import statsmodels.api as sm
import os
from sys import argv

def read_chrom(chrom_file):
    """
    Read a list of chromosomes.
    Args:
       chrom_file (str): name of the file (ex: 'dm3.chrom.sizes')
    Returns:
       chrom_list (list): list of chromosome names
    """
    chrom_list = []
    with open(chrom_file) as f:
        for line in f:
            chrom_list.append(line.strip().split("\t")[0])
    return chrom_list

def read_tad(directory, file_name, chrom_list):
    """
    Read a TAD coordinates .
    Args:
       directory (str): directory of the file location (ex: '/Users/kimm/')
       file_name (str): name of the file (ex: 'tad.bed')
    Returns:
       in_tad (list): [chrom, start, end, tad_id, tad_len, tot_tad_in_chrom, tot_tad_len_in_chrom]
       tad_dict (dictionary): dictionary of tad coordinates by tad_id
    """
    in_tad = []
    tad_dict = {}
    with open(directory + file_name) as f:
        for line in f:
            tmplist = line.strip().split("\t")
            tmplist[1] = int(tmplist[1])
            tmplist[2] = int(tmplist[2])
            tmplist[3] = 'T'+tmplist[3].split('T')[1].zfill(3)
            tmplist.append(tmplist[2]-tmplist[1])
            in_tad.append(tmplist)
            tad_dict[tmplist[3]] = tmplist[0]+':'+str(tmplist[1])+'-'+str(tmplist[2])
    tad_summ = [[x] for x in chrom_list]
    for k in range(len(tad_summ)):
        subset = [x for x in in_tad if x[0]==tad_summ[k][0]]
        tad_summ[k].append(len(subset))
        tad_summ[k].append(sum([x[4] for x in subset]))
    for i in range(len(in_tad)):
        for m in range(len(tad_summ)):
            if in_tad[i][0]==tad_summ[m][0]:
                in_tad[i].extend(tad_summ[m][1:])
    return in_tad, tad_dict

def read_interactions(directory, annot_file, tad_coord):
    """
    Read a master file with TAD annotation.
    Args:
       directory (str): directory of the file location (ex: '/Users/kimm/')
       annot_file (str): name of the file 
       tad_coord (list): tad coordinates list
    Returns:
       tad_comb_list (list): unique tad combinations with counts
    """
    with open(directory + annot_file) as f:
        next(f)
        tad_comb = []
        for line in f:
            tad_ann = line.strip().split("\t")[13]
            tad_list = ['T'+x.split('T')[1].zfill(3) for x in tad_ann.split(";") if x != '-']
            tad_dict = {x:tad_list.count(x) for x in tad_list}
            sub_dict = dict((k,v) for k,v in tad_dict.items() if v > 1)
            if len(sub_dict) > 1:
                tad_str = ','.join([key for key, value in sub_dict.items()])
                tad_comb.append(tad_str)
    tad_comb_dict = {x:tad_comb.count(x) for x in tad_comb}
    tad_comb_list = []
    for key, value in tad_comb_dict.items():
        tad_indx = int(key.split(',')[0][1:])-1
        chrom = tad_coord[tad_indx][0]
        num_tad = len(key.split(','))
        tad_comb_list.append([key, value, chrom, num_tad])
    return tad_comb_list

def compute_sig(intr_list, tad_coord, intr_num, chrom_list, pval_thresh):
    """
    Compute significance for <intr_num> TADs.
    Args:
       intr_list (list): list of interactions
       tad_coord (list): tad coordinates list
       intr_num (int): number of TADs in interactions
       chrom_list (list): list of chromosome names
       pval_thresh (float): significance cut-off threshold (ex: 0.05)
    Returns:
       None
    """
    subset_tads = [x for x in intr_list if x[3]==intr_num]
    gem_summ = {key: 0 for key in chrom_list}
    for key in gem_summ.keys():
        gem_summ[key] += sum([x[6] for x in subset_tads if x[2]==key])
    comb_summ = {key: 0 for key in chrom_list}
    for key in comb_summ.keys():
        comb_summ[key] += len([x for x in subset_tads if x[2]==key])
    for k in range(len(subset_tads)):
        subset_tads[k].append(gem_summ[subset_tads[k][2]])
        subset_tads[k].append(comb_summ[subset_tads[k][2]])
    for k in range(len(subset_tads)):
        tmp = subset_tads[k]
        pval_uniobs = scipy.stats.binom_test(tmp[6], n=tmp[7], p=1/tmp[8], alternative='greater')
        subset_tads[k].append(pval_uniobs)
        if pval_uniobs < pval_thresh:
            subset_tads[k].append('PASS')
        else:
            subset_tads[k].append('FAIL')

if __name__ == '__main__':

    ### Set directory and input file name ###
    
    directory = argv[1] # ex: '/Users/kimm/Documents/MultiChIA/'
    chrom_file = argv[2] # ex: 'dm3.chrom.sizes'
    tad_file = argv[3] # ex: 'GSM3347523_FDR_0.1_ratiothresh_2_pseudoGEM_100000_distTest_PASS_subgem_cov_gem_wgt.1000bp_binned_TAD.bed'
    annot_file = argv[4] # ex: 'GSM3347523_FDR_0.1_ratiothresh_2_pseudoGEM_100000_distTest_master_PASS_annot.txt'
    prefix = argv[5] # ex: 'GSM3347523'
    pval_thresh = float(argv[6]) # ex: 0.05
    out_file = prefix + '_interTAD_BinomTest_sig.tsv'
    
    #### Log file ####
    out = open(directory + prefix + "_BinomTest_logFile.txt", "a")
    
    out.write("Software version: v0.1 (2019-05-21, Kim)" + "\n")
    out.write("Directory: " + directory + "\n")
    out.write("Chrom file name: " + chrom_file + "\n")
    out.write("TAD file name: " + tad_file + "\n")
    out.write("Library name: " + prefix + "\n")
    out.write("p-value threshold: " + str(pval_thresh) + "\n")
    out.write("Started processing frequency-based binomial test for inter-TAD contacts. \n")
    out.write("================================= \n")
    
    ### Read input GEM file ###
    chrom_list = read_chrom(chrom_file)
    tad_coord, tad_dictionary = read_tad(directory, tad_file, chrom_list)
    tad_intrx = read_interactions(directory, annot_file, tad_coord)
    
    out.write("Finished reading files. \n")
    out.write("Chromosomes: " + ','.join(map(str, chrom_list)) + '\n')
    out.write("Total " + str(len(tad_coord)) + " TADs." + "\n")
    out.write("Total " + str(len(tad_intrx)) + " combinations of TAD interactions.\n")
    out.write("================================= \n")
    
    tad_intrx2 = tad_intrx
    
    for i in range(len(tad_intrx2)):
        tads = tad_intrx2[i][0].split(",")
        all_pairs = list(it.combinations(tads,2))
        gem_cnts = []
        for k in range(len(all_pairs)):
            gem_cnts.extend([x[1] for x in tad_intrx if all_pairs[k][0] in x[0] and all_pairs[k][1] in x[0] and x[3]==len(tads)])
        tad_intrx2[i].append(sum(gem_cnts))
        gem_count2 = [x[1] for x in tad_intrx if x[3]>len(tads) and set(tads) < set(x[0].split(","))]
        tad_intrx2[i].append(sum(gem_count2))
    
    for x in tad_intrx2:
        x.append(int(x[4]+x[5]*x[3]*(x[3]-1)/2))
    
    tad_num = [x[3] for x in tad_intrx2]
    for i in list(set(tad_num)):
        compute_sig(tad_intrx2, tad_coord, i, chrom_list, pval_thresh)
        subset = [x for x in tad_intrx2 if x[3]==i]
        subset_pass = [x for x in subset if x[10] == 'PASS']
        subset_fail = [x for x in subset if x[10] == 'FAIL']
        out.write("== Interactions among " + str(i) + " TADs == \n")
        out.write("Total: " + str(len(subset)) + ' combinations. \n')
        out.write("Pass: " + str(len(subset_pass)) + ' combinations. \n')
        if len(subset_pass) > 0:
            out.write("Pass avg. complex cnt: " + str(round(np.mean([x[1] for x in subset_pass]), 2)) + ' \n')
            out.write("Pass avg. combs in higher class: " + str(round(np.mean([x[5] for x in subset_pass]), 2)) + ' \n')
            out.write("Pass avg. norm. complex cnt : " + str(round(np.mean([x[6] for x in subset_pass]), 2)) + ' \n')
        
        out.write("Fail: " + str(len(subset_fail)) + ' combinations. \n')
        if len(subset_fail) > 0:
            out.write("Fail avg. complex cnt: " + str(round(np.mean([x[1] for x in subset_fail]), 2)) + ' \n')
            out.write("Fail avg. combs in higher class: " + str(round(np.mean([x[5] for x in subset_fail]), 2)) + ' \n')
            out.write("Fail avg. norm. complex cnt: " + str(round(np.mean([x[6] for x in subset_fail]), 2)) + ' \n')
    
    for j in range(len(tad_intrx2)):
        tadlist = tad_intrx2[j][0].split(',')
        tad_intrx2[j].append(','.join([tad_dictionary[t] for t in tadlist]))
    
    ### Write output file ###
    header = ['TAD combination', 'ComplexCnt', 'Chrom', '# of TADs', 'Pairs in same class', 'Combs in higher class', 'Norm. cnt', 'Tot. norm. cnts by tadn in chrom', 'Tot. # of combs by tadn in chrom', 'p-val (uniform obs.)', 'Decision', 'TAD coord.']
    with open(directory + out_file, 'a') as file1:
        file1.write('\t'.join(map(str,header)) + '\n')
        for i in range(len(tad_intrx2)):
            file1.write('\t'.join(map(str,tad_intrx2[i])) + '\n')
    file1.close()
    out.write("================================= \n")
    out.write("DONE. \n")
    out.close()

