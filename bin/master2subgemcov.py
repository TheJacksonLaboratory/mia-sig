import numpy as np
import random
import itertools as it
import os
from sys import argv

def read_mf(directory, file_name, category):
    """
    Read a master file after significance test.
    Args:
       directory (str): directory of the file location (ex: '/Users/kimm/')
       file_name (str): name of the file (ex: 'SHG0008H.Fragnum_PlinePgem')
    Returns:
       in_gems (list): tab-separated lists of lists
    """
    in_gems = []
    with open(directory + file_name) as f:
        next(f) # skip the header
        for line in f:
            dat = line.strip().split("\t")
            if category == 'ALL':
                in_gems.append(dat)
            else:
                if dat[9]==category or dat[12]==category:
                    in_gems.append(dat)
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

def split_gems(frags, thr):
    """
    Split GEMs into sub-GEMs depending on cutoff threshold thr.
    Args:
       frags (list of list): [chrom,start,end] for each fragment 
       thr (int): cutoff threshold in bps
    Returns:
       subgems (list of list)
    """
    subgems = []
    tmpgem = [frags[0]]
    for i in range(len(frags)-1):
        fragdist = frags[i+1][1]-frags[i][2]
        if fragdist < thr:
            tmpgem.append(frags[i+1])
        else:
            subgems.append(tmpgem)
            tmpgem = [frags[i+1]]
    subgems.append(tmpgem)
    fin_subgems = [x for x in subgems if len(x)>1]
    return fin_subgems

def all_pairs(frags, thr):
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
        gem_cov = allPairs[j][0][0:2]
        gem_cov.append(allPairs[j][1][2])
        if thr > 0:
            if (gem_cov[2]-gem_cov[1]) < thr:
                pair_all.append(gem_cov)
        else:
            pair_all.append(gem_cov)
    return pair_all

def linear_tri(frags, thr):
    """
    Triangle (pyramid) of linear increase then linear decrease weights.
    Args:
       frags (list of list): [chrom,start,end] for each fragment in a GEM
    Returns:
       lin_all (list of bed entries): bed entries in format [chrom,start,end]
    """
    lin_all = []
    n = len(frags)
    for i in range(int(n/2)):
        gem_cov = frags[i][0:2]
        gem_cov.append(frags[n-i-1][2])
        if thr > 0:
            if (gem_cov[2]-gem_cov[1]) < thr:
                lin_all.append(gem_cov)
        else:
            lin_all.append(gem_cov)
    return lin_all

def gem_cov(frags, weight_flag):
    """
    GEM coverage.
    Args:
       frags (list of list): [chrom,start,end] for each fragment in a GEM
       weight_flag (binary): True if weighted by fragment number; False otherwise
    Returns:
       gem_all (list of bed entries): bed entries in format [chrom,start,end]
    """
    gem_all = frags[0][0:2]
    gem_all.append(frags[-1][2])
    if weight_flag == True:
        gem_fin = [gem_all] * len(frags)
    else:
        gem_fin = [gem_all]
    return gem_fin

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

def tadsize_chart(genome_name):
    """
    Determine the distance threshold to build coverage tracks.
    Args:
       genome_name (string): name of the reference genome; 
                           ex: mammals, drosophila, c_elegans, s_pombe, c_crescentus
    Returns:
       dist_thresh (int): integer specifying distance threshold in basepairs
    """
    low_bound = {
        "mammals": 100000,
        "drosophila": 10000,
        "c_elegans": 1000000,
        "s_pombe": 50000,
        "c_crescentus": 30000
    }
    upp_bound = {
        "mammals": 2000000,
        "drosophila": 100000,
        "c_elegans": 2000000,
        "s_pombe": 100000,
        "c_crescentus": 400000
    }
    typ_res = {
        "mammals": 1000000,
        "drosophila": 250000,
        "c_elegans": 3000000,
        "s_pombe": 300000,
        "c_crescentus": 250000
    }
    return low_bound[genome_name], upp_bound[genome_name], typ_res[genome_name]

if __name__ == '__main__':

    ### Set directory and input file name ###
    directory = argv[1]
    file_name = argv[2]
    category = argv[3] # options: 'PASS', 'FAIL', 'ALL'
    prefix = file_name.split("_master.txt")[0] + '_' + category + '_'
    species = argv[4] # options: 'mammals', 'drosophila', etc.
    thr = tadsize_chart(species)[1]

    #### Log file ####
    out = open(directory + prefix + "gemcov_logFile.txt", "a")
    
    out.write("Directory: " + directory + "\n")
    out.write("File name: " + file_name + "\n")
    out.write("Category: " + category + "\n")
    out.write("Species: " + species + "\n")
    out.write("Threshold: " + str(thr)+ "\n")
    out.write("================================= \n")
    
    ### Read input master file ###
    tot_gems = read_mf(directory, file_name, category)
    out.write("Finished reading the input master file. \n")
    out.write(str(len(tot_gems)) + " total complexes. \n")
    out.write("================================= \n")

    for i in range(len(tot_gems)):
        frags = extract_frags(tot_gems[i])
        subgems = split_gems(frags, thr) # threshold
        for j in range(len(subgems)):
            gcov_wgt = gem_cov(subgems[j], True)
            for n in range(len(gcov_wgt)):
                write_output_file(gcov_wgt[n], directory+prefix+'subgem_cov_gem_wgt.bed')

    out.write("DONE. \n")
    out.close()
