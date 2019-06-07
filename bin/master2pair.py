import numpy as np
import random
import itertools as it
from sys import argv
import os

def read_mf(directory, file_name):
    """
    Read a master file after significance test.
    Args:
       directory (str): directory of the file location (ex: '/Users/kimm/')
       file_name (str): name of the file (ex: 'SHG0008H.Fragnum_PlinePgem')
    Returns:
       in_gems (list): tab-separated lists of lists
    """
    with open(directory + file_name) as f:
        next(f) # skip the first line, since it's a header
        in_gems = [line.strip().split("\t") for line in f]
    return in_gems

def extract_frags(gem_list):
    """
    Extract necessary information from input gem.
    Args:
       gem_list (list): list with 10 items as encoded in PlinePgem file
    Returns:
       frags (list of list): [chrom,start,end] for each fragment 
    """
    raw_frags = gem_list[4].split(";")
    frags = []
    for j in range(len(raw_frags)):
        bedentry = raw_frags[j].split("(")[0].replace("-", ":").split(":")
        bedentry[1] = int(bedentry[1])
        bedentry[2] = int(bedentry[2])
        frags.append(bedentry)
    return frags

def all_pairs(frags, plec_weight, plisrs_iter, plisrs_dedup):
    """
    Enumerate all pairs.
    Args:
       frags (list of list): [chrom,start,end] for each fragment in a GEM
    Returns:
       pair_all (list of pair): all pairs of frags in format [pair1, pair2]
    """
    allPairs = list(it.combinations(frags, 2))
    pair_all = []
    for j in range(len(allPairs)):
        pair1 = allPairs[j][0]
        pair2 = allPairs[j][1]
        pair_all.append([pair1, pair2, -1])
    return pair_all

def pairs2other(pair_list, file_type):
    """
    Conbert pairs to short.
    Args:
       pair_list (list of pairs): [[chrom,start,end], [chrom,start,end], weight] for each interaction
       file_type (string): 'short' or 'bedpe'
    Returns:
       short (list of list): all pairs of frags in short format (0, chrom1, mid1, 0, 0, chrom2, mid2, 1)
    """
    pair1 = pair_list[0]
    pair2 = pair_list[1]
    weight = pair_list[2]
    other_list = []
    if file_type == 'bedpe':
        other_list.extend(pair1)
        other_list.extend(pair2)
    elif file_type == 'short':
        other_list = [0, pair1[0], int((pair1[1] + pair1[2])/2), 0, 0, pair2[0], int((pair2[1] + pair2[2])/2), 1]
        if weight != -1:
            other_list.append(weight)
    return other_list

def run_plec(frags, plec_weight, plisrs_iter, plisrs_dedup):
    """
    Run PLEC.
    Args:
       frags (list of list): [chrom,start,end] for each fragment in a GEM
       plec_weight (string): factor to add as a function of n, the number of fragments in a GEM
                                ex: '4/(2*n-1)'
    Returns:
       short_plec (list of list): all pairs of frags in short format (0, chrom1, mid1, 0, 0, chrom2, mid2, 1, weight)
    """
    allPairs = list(it.combinations(frags, 2))
    pair_plec = []
    n = len(frags)
    for j in range(len(allPairs)):
        pair1 = allPairs[j][0]
        pair2 = allPairs[j][1]
        pair_plec.append([pair1, pair2, eval(plec_weight)])
    return pair_plec

def run_plisrs(frags, plec_weight, plisrs_iter, plisrs_dedup):
    """
    Run PLISRS.
    Args:
       frags (list of list): [chrom,start,end] for each fragment in a GEM
       iteration (int): number of iterations 
       dedup (boolean): indicator for de-duplication (options: True or False)
    Returns:
       short_pisrs (list of list): all pairs of frags in short format (0, chrom1, mid1, 0, 0, chrom2, mid2, 1, weight)
    """
    pair_dict = {}
    pair_plisrs = []
    n = len(frags)
    for i in range(plisrs_iter):
        frags_double = np.repeat([i for i in range(n)], 2, axis = 0)
        np.random.shuffle(frags_double)
        for y in range(n):
            sorted_pair = ','.join(str(e) for e in sorted(frags_double[2*y:2*y+2]))
            pair_dict[sorted_pair] = pair_dict.setdefault(sorted_pair, 0) + 1
        if plisrs_dedup == True:
            pair_dict = dict.fromkeys(pair_dict, 1)
    for key in pair_dict:
        pair1 = frags[int(key.split(",")[0])]
        pair2 = frags[int(key.split(",")[1])]
        pair_plisrs.append([pair1, pair2, pair_dict[key]])
    return pair_plisrs

def write_output_file(pair_list, out_name):
    """ 
    Write out fragments in short format (regardless PLEC, PLISRS, all pairs).
    Args: 
       pair_list (list): list of pairs of fragments in short or bedpe format
       out_name (string): output file name
    Returns:
       None
    """
    with open(out_name, 'a') as file1:
        file1.write('\t'.join(map(str, pair_list)) + '\n')
    file1.close()

if __name__ == '__main__':

    ### Set directory and input file name ###
    directory = argv[1] # ex: '/Users/kimm/Documents/MultiChIA/'
    file_name = argv[2] # ex: 'P2MC7N8HCE3K_chr2L_FDR_0.1_pseudoGEM_100_enrichTest_master.txt'
    catgr = argv[3] # options: 'PASS', 'FAIL', 'ALL'
    method = argv[4] # options: 'allpairs', 'plec', 'plisrs'
    bedpe_flag = argv[5] # options: 'True', 'False'
    prefix = file_name.split("_master.txt")[0]+'_'+catgr+'_'+method
    
    # method dictionary #
    call_method = {
        'allpairs': all_pairs,
        'plec': run_plec,
        'plisrs': run_plisrs
    }
    
    plec_weight = '4/(2*n-1)'
    plisrs_iter = 0
    plisrs_dedup = False
    
    #### Remove existing file ####
    if os.path.exists(directory+prefix+'.short'):
        os.remove(directory+prefix+'.short')
    if os.path.exists(directory+prefix+'.bedpe'):
        os.remove(directory+prefix+'.bedpe')
    
    #### Log file ####
    out = open(directory + prefix + "_logFile.txt", "a")
    
    out.write("Software version: v0.1 (2019-04-23, Kim)" + "\n")
    out.write("Directory: " + directory + "\n")
    out.write("File name: " + file_name + "\n")
    out.write("================================= \n")
    
    ### Read input master file ###
    tot_gems = read_mf(directory, file_name)
    out.write("Finished reading the input master file. \n")
    out.write(str(len(tot_gems)) + " total complexes. \n")
    out.write("================================= \n")
    
    ### Subset GEMs ###
    if catgr != "ALL":
        in_gems = [x for x in tot_gems if x[9]==catgr or x[12]==catgr]
    else:
        in_gems = tot_gems
    out.write("Finished subsetting master file. \n")
    out.write(str(len(in_gems)) + " complexes in the " + catgr + " category. \n")
    out.write("================================= \n")
    
    fragnum = []
    for k in range(len(in_gems)):
        frags = extract_frags(in_gems[k])
        fragnum.append(len(frags))
        pairs = call_method[method](frags, plec_weight, plisrs_iter, plisrs_dedup)
        for j in range(len(pairs)):
            pair_short = pairs2other(pairs[j], 'short')
            write_output_file(pair_short, directory+prefix+'.short')
            if bedpe_flag == 'True':
                pair_bedpe = pairs2other(pairs[j], 'bedpe')
                write_output_file(pair_bedpe, directory+prefix+'.bedpe')

    exp_all = int(sum([x*(x-1)/2 for x in fragnum]))
    
    out.write("Finished converting " +str(k+1)+ " GEMs to pairs using the method " + method + ". \n")
    out.write("Expected number of lines for allpairs/plec method: " + str(exp_all)+". \n")
    if method in ('allpairs','plec') and sum(1 for line in open(directory+prefix+'.short'))==exp_all:
        out.write("Verified that the line number matches the expectation.")
    out.write("DONE. \n")
    out.close()
