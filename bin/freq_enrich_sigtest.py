import numpy as np
import random
import statsmodels.api as sm
import os
from statsmodels.sandbox.stats.multicomp import multipletests
import pybedtools
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
        in_gems = [line.strip().split("\t") for line in f]
    return in_gems

def read_chroms(directory, genome_name):
    """
    Read a tab-delimited text file with list of chromosomes and their sizes.
    Args:
       directory (str): directory of the file location (ex: '/Users/kimm/')
       genome_name (str): name of reference genome (ex: 'dm3', 'mm10', 'hg38', 'hg19', etc.)
                        Note: must have a text file <genome_name>.chrom.sizes in the directory.
    Returns:
       chrom_dict (dictionary): tab-separated lists of lists
    """
    chrom_dict = {}
    with open(directory + genome_name + '.chrom.sizes') as f:
        for line in f:
            tmp_list = line.strip().split("\t")
            chrom_dict[tmp_list[0]] = int(tmp_list[1])
    return chrom_dict

def read_bedgraph(directory, library_name, bin_size, chrom):
    """
    Read a bedgraph file averaged over fixed windows.
    Args:
       directory (str): directory of the file location (ex: '/Users/kimm/')
       library_name (str): name of the file (ex: 'P2MC7N8HCE3K_binned_10bp_')
       bin_size (int): fixed bin size in which bedgraph was generated (ex: 10)
    Returns:
       cov_list (list): list recording only 4th column of input file (i.e., coverage values)
    """
    with open(directory + library_name + '_binned_' + str(bin_size) + 'bp_' + chrom + '.bedgraph') as f:
        cov_list = [int(float(line.strip().split("\t")[3])) for line in f]
    return cov_list

def get_mean_cov(start, end, coverage, bin_size):
    """
    Extract mean (average) coverage of a given region.
    Args:
       start (int): start location of region (ex: 250)
       end (int): end location of region (ex: 400)
       coverage (str): coverage list name
       bin_size (int): bin size used for bedgraph generation
    Returns:
       mean_cov (int): mean coverage value
    """
    interval_values = coverage[int(start/bin_size):int(end/bin_size)+1]
    return sum(interval_values)/len(interval_values)

def random_loc(chrom_size, gem_span, sample_size):
    """
    Randomly locate GEM with same span in a chromosome sample_size times.
    Args:
       chrom_size (int): size of a given chromosome (ex: 23011544 if 'chr2L')
       gem_span (int): GEM from start of first fragment to end of last fragment (ex: 5000)
       sample_size (int): number of random locations to sample
    Returns:
       startpos (array): random integers array of length 'sample_size'
    """
    np.random.seed(12345)
    startpos = np.random.randint(low = 0, high = chrom_size - gem_span - 1, size = sample_size)
    return(startpos)

def pseudo_gem(start_loc, frag_lengths, f2f_dist):
    """
    Create a pseudo-GEM on a random start location.
    Args:
       start_loc (int): random start location on the genome (ex: 390)
       frag_lengths (list): lengths of each fragments (ex: [100, 200, 150])
       f2f_dist (list): fragment to fragment distances (ex: [1000, 4000])
    Returns:
       pseudo_gem (list): start and end of fragments in a list
                     (ex: [390,490,1490,1690,5690,5840])
    """
    f2f_dist.append(0)
    pseudo_gem = []
    left_coord = start_loc
    for i in range(len(frag_lengths)):
        right_coord = left_coord + frag_lengths[i]
        pseudo_gem.extend([left_coord, right_coord])
        left_coord = right_coord + f2f_dist[i]
    return(pseudo_gem)

def extract_info(gem_list):
    """
    Extract necessary information from input gem.
    Args:
       gem_list (list): list with 10 items as encoded in PlinePgem file
    Returns:
       gem_id (str): GEM ID (ex: SHG0008H-1000-10000799-AAACACCCACTAGTAC-K00384-FA-1-0)
       chrom (str): chromosome name (ex: chr2L)
       frags (list of list): [start,end] for each fragment 
                        (ex: [3251935, 3252561], [3413857, 3414485])
       span (int): start of leftmost fragment to end of rightmost fragment (ex: 162550)
       fraglen (list): length of each fragment (ex: [626, 628])
       f2fdist (list): distance between neighboring fragments (ex: [161296])
       fragnum (int): number of fragments
       fragstr (str): fragment coordinates stringi
       coord (str): gem coordinates
    """
    raw_frags = gem_list[4].split(";")
    frags = []
    for j in range(len(raw_frags)):
        bedentry = raw_frags[j].split("(")[0].replace("-", ":").split(":")
        frags.append([int(x) for x in bedentry[1:]])
    gem_id = gem_list[0]
    chrom = bedentry[0]
    span = frags[-1][1] - frags[0][0]
    fraglen = [x[1]-x[0] for x in frags]
    f2fdist = [frags[x][0]-frags[x-1][1] for x in range(1,len(frags))]
    fragnum = j + 1
    fragstr = ';'.join([chrom+':'+str(x[0])+'-'+str(x[1]) for x in frags])
    coord = chrom+':'+str(frags[0][0])+'-'+str(frags[-1][1])
    return gem_id, chrom, frags, span, fraglen, f2fdist, fragnum, fragstr, coord

def get_raw_pval(obs_frags, obs_span, obs_fraglen, obs_f2f, sample_size, cov, bin_size, chrom_size):
    """
    Compute raw p-value for a GEM.
    Args:
       obs_frags (list of list): [start,end] for each fragment 
       obs_span (int): start of leftmost fragment to end of rightmost fragment
       obs_fraglen (list): length of each fragment
       obs_f2f (list): distance between neighboring fragments
       sample_size (int): number of pseudo-GEMs to sample
       cov: bedgraph coverage track
       bin_size (int): size of the bin used to generate cov
       chrom_size (int): size of the chromosome
    Returns:
       raw_pval (float): raw p-value computed
       obs_enrich (int): observed enrichment
    """
    # observed enrichment
    frag_cov = []
    for i in range(len(obs_frags)):
        frag_cov.append(get_mean_cov(obs_frags[i][0], obs_frags[i][1], cov, bin_size))
    obs_enrich = sum(frag_cov)/len(frag_cov)
    
    # expected enrichment (sampling background)
    startpos = random_loc(chrom_size, obs_span, sample_size)
    exp_enrich = []
    for k in range(sample_size):
        pseudo = pseudo_gem(startpos[k], obs_fraglen, obs_f2f)
        exp = []
        for j in range(int(len(pseudo)/2)):
            exp.append(get_mean_cov(pseudo[2*j], pseudo[2*j+1], cov, bin_size))
        exp_enrich.append(sum(exp)/len(exp))
    raw_pval = sum(i > obs_enrich for i in exp_enrich)/sample_size
    return(raw_pval, obs_enrich)

def get_adj_pval(raw_pval_list, fdr_thresh, method):
    """ 
    Adjust raw pvalues by Benjamini Hochberg multiple testing adjustment. 
    Args:
       raw_pval_list (list): list of raw p-values (ex: [0.1, 0.04, 0.1])
       fdr_thresh (float): false discovery rate (ex: 0.05)
       method (string): adjustment method (ex: 'fdr_bh')
    Returns:
       adj_pval_list (array): array of booleans and adjusted p-values (ex: [0.1, 0.1, 0.1])
    """
    adj_pval_list = multipletests(raw_pval_list, alpha = fdr_thresh, method = method)
    return(adj_pval_list)

def write_master_result(out_gem_list, out_name):
    """ 
    Write out significance test results.
    Args: 
       out_gem_list (list): list with 11 items including gem_id, frag_str, p-values, etc.
       out_name (string): output file name
    Returns:
       None
    """
    with open(out_name, 'a') as file1:
        header = ['GEM_ID', 'GEM_coord', 'GEM_span', 'Frag_number', 'List_of_frag_coord', 'Category', 'Obs', 'rawpval1', 'adjpval1', 'decis1', 'rawpval2', 'adjpval2', 'decis2']
        file1.write('\t'.join(map(str, header)) + '\n')
        for i in range(len(out_gem_list)):
            file1.write('\t'.join(map(str, out_gem_list[i])) + '\n')

if __name__ == '__main__':
    ### Set parameters ###
    library_name = argv[1] ## Library name of our data ##
    genome_name = argv[2] ## Name of the reference genome ##
    fdr_thresh = float(argv[3])  # should be argument; Benjamini-Hochberg FDR; p-value cutoff ##
    chr_name = argv[4] # should be argument
    samp_size = int(argv[5]) ## Number of pseudo-GEMs ##
    bin_size = int(argv[6])
    prefix = library_name + "_" + chr_name + "_FDR_" + str(fdr_thresh) + "_pseudoGEM_" + str(samp_size)
    
    ### Set directory and input file name ###
    directory = argv[7]
    file_name = argv[8]
    out_directory = directory + library_name + "_EnrichTest_FDR_" + str(fdr_thresh) + '/'
    if not os.path.exists(out_directory):
        os.mkdir(out_directory)
    
    #### Log file ####
    out = open(out_directory + prefix + "_enrichTest_logFile.txt", "a")
    
    out.write("Software version: v0.1 (2019-04-23, Kim)" + "\n")
    out.write("Input directory: " + directory + "\n")
    out.write("Input file name: " + file_name + "\n")
    out.write("Output directory: " + out_directory + "\n")
    out.write("Library name: " + library_name + "\n")
    out.write("Reference genome: " + genome_name + "\n")
    out.write("FDR threshold: " + str(fdr_thresh) + "\n")
    out.write("Number of pseudo-GEMs: " + str(samp_size) + "\n")
    out.write("Started processing frequency-based enrichment test. \n")
    out.write("================================= \n")
    out.write("===== Chromosome is: " + chr_name + ". ===== \n")
    out.write("================================= \n")
    
    ### Read input GEM file ###
    in_gems = read_gems(directory, file_name)
    out.write("Finished reading the input GEM file. \n")
    out.write(str(len(in_gems)) + " total GEMs. \n")
    out.write("================================= \n")
    
    ### Read chrom sizes file ###
    chrom_size = read_chroms(directory, genome_name)
    
    ### Subset input GEM files by chromosome ###
    subset_gems = [x for x in in_gems if chr_name == x[1].split(":")[0]]
    out.write("Finished subsetting GEMs. \n")
    out.write(str(len(subset_gems)) + " GEMs in " + chr_name + " (of " + str(len(in_gems)) + " total GEMs). \n")
    out.write("================================= \n")
    
    del in_gems
    
    ### Read bedgraph coverage file ###
    bg_name = library_name + '_binned_' + str(bin_size) + 'bp_' + chr_name + '.bedgraph'
    if not os.path.isfile(directory+bg_name):
        out.write("Error: no bedgraph coverage file." + "\n")
    coverage = read_bedgraph(directory, library_name, bin_size, chr_name)
    out.write("Finished reading the coverage file " + library_name + '_binned_' + str(bin_size) + 'bp_' + chr_name + '.bedgraph'+ ". \n")
    out.write("================================= \n")

    ### Initialize output GEM list ###
    out_gems = [] 
    
    ### Calculate raw p-values for subset_gems ###
    pvals = []
    #for k in range(100):
    for k in range(len(subset_gems)):
        gem_id, chrom, frags, span, fraglen, f2fdist, fragnum, frag_str, coord = list(extract_info(subset_gems[k]))
        pval, enr = get_raw_pval(frags, span, fraglen, f2fdist, samp_size, coverage, bin_size, chrom_size[chr_name])
        pvals.append(pval)
        out_gems.append([gem_id, coord, span, fragnum, frag_str, 'Orig', round(enr,1), round(pval,3)])
    out.write("Finished calculating raw p-values for " + str(k + 1) + " GEMs. \n")
    raw_pass = sum(i <= fdr_thresh for i in pvals)
    out.write(str(raw_pass) + " GEMs have raw p-val <= " + str(fdr_thresh) + "(i.e., " + "{0:.2f}".format(raw_pass/(k+1)*100) +" %). \n")
    out.write("================================= \n")
    
    ### Adjust p-values ###
    adj_pval = get_adj_pval(pvals, fdr_thresh, 'fdr_bh')
    out.write("Finished adjusting p-values for " + str(k + 1) + " GEMs. \n")
    adj_pass = 0
    for i in range(len(adj_pval[1])):
        out_gems[i].append(round(adj_pval[1][i], 3))
        if adj_pval[0][i]:
            out_gems[i].extend(['PASS','.','.','.'])
            adj_pass += 1
        else:
            out_gems[i].extend(['FAIL','.','.','.'])
    out.write(str(adj_pass) + " GEMs have adjusted p-val <= " + str(fdr_thresh) + "(i.e.," + "{0:.2f}".format(adj_pass/(k+1)*100) +" %). \n")
    out.write("================================= \n")
    
    ### Write results ###
    write_master_result(out_gems, out_directory+prefix+'_enrichTest_master.txt')
    out.write("Finished writing files. \n")
    out.write("DONE. \n")
    out.close()
