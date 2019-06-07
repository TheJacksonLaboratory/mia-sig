import os
from sys import argv

def read_gems(directory, file_name):
    """
    Read a GEM file of the form "PlinePgem".
    Args:
       directory (str): directory of the file location (ex: '/Users/kimm/')
       file_name (str): name of the file (ex: 'SHG0008H.Fragnum_PlinePgem')
    Returns:
       in_gems (list): tab-separated lists of lists
    """
    with open(directory + file_name) as f:
        next(f)
        in_gems = [line.strip().split("\t") for line in f]
    return in_gems

def extract_frags_gemid(gem_list):
    """
    Extract necessary information from input gem.
    Args:
       gem_list (list): list with 10 items as encoded in PlinePgem file
    Returns:
       frags (list of list): [chrom,start,end] for each fragment 
       gem_id (string): GEM ID 
    """
    raw_frags = gem_list[4].split(";")
    gem_id = gem_list[0]
    frags = []
    for j in range(len(raw_frags)):
        bedentry = raw_frags[j].replace("-", ":").split(":")
        frags.append(bedentry)
    return frags, gem_id

def write_frag_file(frag_result, gem_id, out_name):
    """ 
    Write out fragments in fragment format to be visualized in ChIA-viewer.
    Args: 
       frag_result (list of list): [chrom,start,end] for each fragment 
       gem_id (string): GEM ID
       out_name (string): output file name
    Returns:
       None
    """
    n = len(frag_result)
    with open(out_name, 'a') as file1:
        for i in range(len(frag_result)):
            file1.write('\t'.join(map(str, frag_result[i])) + '\t')
            file1.write(str(n) + '\t')
            file1.write(gem_id + '\n')
    file1.close()

if __name__ == '__main__':

    ### Set directory and input file name ###
    directory = argv[1]
    master_file = argv[2]
    prefix = argv[3]

    ### Read input GEM file ###
    in_gems = read_gems(directory, master_file)

    for k in range(len(in_gems)):
        frags, gem_id = extract_frags_gemid(in_gems[k])
        write_frag_file(frags, gem_id, directory + prefix + '.region')

