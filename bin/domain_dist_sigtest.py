import numpy as np
import random
import statsmodels.api as sm
import itertools as it
import os
from statsmodels.sandbox.stats.multicomp import multipletests
from itertools import compress
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
        for line in f:
            entries = line.strip().split("\t")
            frags = entries[4]
            gem_id = entries[0]
            in_gems.append([frags+'|'+gem_id])
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
       fragstr (str): fragment coordinates string
       coord (str): gem coordinates
    """
    raw_frags = gem_list[0].split("|")[0].split(";")
    frags = []
    for j in range(len(raw_frags)):
        bedentry = raw_frags[j].split("(")[0].replace("-", ":").split(":")
        frags.append([int(x) for x in bedentry[1:]])
    if '|' in gem_list[0]:
        gem_id = gem_list[0].split("|")[1]
    else:
        gem_id = ''
    chrom = bedentry[0]
    span = frags[-1][1] - frags[0][0]
    fraglen = [x[1]-x[0] for x in frags]
    f2fdist = [frags[x][0]-frags[x-1][1] for x in range(1,len(frags))]
    fragnum = j + 1
    fragstr = ';'.join([chrom+':'+str(x[0])+'-'+str(x[1]) for x in frags])
    coord = chrom+':'+str(frags[0][0])+'-'+str(frags[-1][1])
    return gem_id, chrom, frags, span, fraglen, f2fdist, fragnum, fragstr, coord

def norm_shannon_ent(prob_list):
    """
    Compute the normalized Shannon entropy of a list of probabilities.
    Args:
       prob_list (list): list of probabilities (ex: [0.1, 0.4, 0.5])
    Returns:
       entropy (float): tab-separated lists of lists (ex: 0.85867)
    """
    if sum(prob_list) != 1:   ## input is count or frequency instead of probability
        prob_list = [i/sum(prob_list) for i in prob_list]
    entropy = sum([x*np.log2(1/x) for x in prob_list])/np.log2(len(prob_list))
    return float(format(entropy, '.5f'))

def create_pseudo_gems(f2f_bucket, fragnum_bucket, samp_size, file_name):
    """
    Create <samp_size> pseudoGEMs for each unique fragment number.
    Args:
       f2f_bucket (list): list of fragment-to-fragment distances, from which samples will be drawn
       fragnum_bucket (list): number of fragments in each GEM for all GEMs
       samp_size (int): number of pseudoGEMs to create
       file_name (str): file to write 1000 pseudoGEMs per frag_num
    Returns:
       pseudo_gem_tot_dist (dictionary): dictionary of a list of <samp_size> total distances for fragnum
       pseudo_gem_mean_ent (dictionary): dictionary of mean entropy for fragnum
    """
    unique_frag_num = sorted(list(set(fragnum_bucket)))
    f2f_bucket_size = len(f2f_bucket)
    bucket_array = np.array(f2f_bucket)
    pseudo_gem_tot_dist = {}
    pseudo_gem_mean_ent = {}
    np.random.seed(12345)
    with open(file_name, 'a') as file1:
        for i in unique_frag_num:
            rand_ind = list(np.random.choice(range(f2f_bucket_size), (i-1)*samp_size))
            count = 0
            sum_ent = 0
            tmp_tot_dist = []
            for j in range(0,len(rand_ind), (i-1)):
                tmp_dist = bucket_array[rand_ind[j:(j+i-1)]].tolist()
                # get normalized Shannon entropy of first 1000 pseudoGEMs if a given GEM has > 2 fragments
                if (count < 1000): 
                    file1.write(str(i) + '\t')
                    file1.write(';'.join(map(str, tmp_dist)) + '\t')
                    if (i > 2):
                        tmp_ent = norm_shannon_ent(tmp_dist)
                        file1.write(str(tmp_ent) + '\n')
                        sum_ent = sum_ent + tmp_ent
                    else:
                        file1.write("NA" + '\n')
                    count += 1
                tmp_tot_dist.append(sum(tmp_dist))
            pseudo_gem_tot_dist[str(i)] = np.sort(tmp_tot_dist)
            pseudo_gem_mean_ent[str(i)] = float(format(sum_ent/1000, '.5f'))
        pseudo_gem_tot_dist['1'] = [] # empty for singleton
    file1.close()
    return pseudo_gem_tot_dist, pseudo_gem_mean_ent


def get_raw_pval(f2f_dist, pseudo_gem_dist, samp_size):
    """
    Get raw p-value by performing a distance test.
    Args:
       f2f_dist (list): fragment-to-fragment distances for a GEM (ex: [446,21231,249])
       pseudo_gem_dist (list): pseudo gem total distances of a given <fragnum>
       samp_size (int): number of pseudoGEMs to create (ex: 10000)
    Returns:
       raw_pval (float): raw p-value from a total distance test; '.' if singleton GEM
    """
    if not f2f_dist:
        pval = '.'
    else:
        tot_dist = sum(f2f_dist)
        count = max(1, np.searchsorted(pseudo_gem_dist, tot_dist))
        pval = count/samp_size
    return pval

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

def entropy_filter(in_gem, pseudo_gem_ent, cutoff_ent_filter):
    """ 
    Filter a GEM into subGEMs if potentially doublet or triplet (metric: norm shannon Entropy). 
    Args:
       in_gem (string): input GEM with fragments and GEM ID (ex: 'chr2L:11-33;chr2L:44-55|tmp_gem_id')
       pseudo_gem_ent (float): mean of normalized Shannon entropies in 1000 pseudoGEM (ex: 0.12)
       cutoff_ent_filter (int): max_dist/sec_max_dist ratio cutoff threshold (ex: 2)
    Returns:
       num_cuts (int): number of cuts made in this filter (ex: 0,1, or 2)
       sub_gems (list of strings): [sub_gem1, sub_gem2, etc.]
    """
    gem_id, chrom, frags, span, fraglen, f2fdist, fragnum, frag_str, coord = list(extract_info(in_gem))
    num_cuts = 0
    sub_gems = []
    obs_ent = norm_shannon_ent(f2fdist)
    if obs_ent < pseudo_gem_ent[str(fragnum)]: # cut once
        max_dist = max(f2fdist)
        max_dist_ind = f2fdist.index(max_dist)
        sec_max_dist = sorted(set(f2fdist))[-2]
        sec_max_dist_ind = f2fdist.index(sec_max_dist)
        sub_frags = [frags[:(max_dist_ind+1)], frags[(max_dist_ind+1):]]
        if max_dist/sec_max_dist < cutoff_ent_filter: # cut twice
            sorted_ind = sorted([max_dist_ind, sec_max_dist_ind])
            sub_frags = [frags[:(sorted_ind[0]+1)], frags[(sorted_ind[0]+1):(sorted_ind[1]+1)], frags[(sorted_ind[1]+1):]]
        num_cuts = len(sub_frags)-1
        for i in range(len(sub_frags)):
            sub_gems.append(';'.join([chrom+':'+str(x[0])+'-'+str(x[1]) for x in sub_frags[i]])+'|'+gem_id+'-sub-'+str(len(sub_frags))+'-'+str(i+1))
    else:
        sub_gems = [frag_str+'|'+gem_id]
    return sub_gems, num_cuts

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
        for i in range(len(out_gem_list)):
            file1.write('\t'.join(map(str, out_gem_list[i])) + '\n')
    file1.close()

if __name__ == '__main__':
    """ 
    Final results with following columns:
    1) GEM ID 2) GEM coordinate 3) GEM span 4) fragnum 5) frag coord 6) orig/postEF 7) totdist 
    8) pval1 9) adj.pval1 10) status1 11) pval2 12) adj.pval2 13) status2
    """
    ### Set parameters ###
    library_name = argv[1] ## Library name of our data ##
    genome_name = argv[2] ## Name of the reference genome ##
    fdr_thresh = float(argv[3])  ## Benjamini-Hochberg FDR; p-value cutoff ##
    cutoff_ent_filter = int(argv[4]) ## Ratio of largest to second largest F2F cutoff in our entropy filter ##
    samp_size = int(argv[5]) ## Number of pseudo-GEMs ##
    prefix = library_name + "_FDR_" + str(fdr_thresh) + "_ratiothresh_"+ str(cutoff_ent_filter) + "_pseudoGEM_"+str(samp_size)

    ### Set directory and input file name ###
    directory = argv[6]
    file_name = argv[7]
    chrom_dir = argv[8]
    out_directory = directory + library_name + "_DistTest_FDR_" + str(fdr_thresh) + '/'
    if not os.path.exists(out_directory):
        os.mkdir(out_directory)

    #### Log file ####
    out = open(out_directory + prefix + "_distTest_logFile.txt", "a")
    
    out.write("Software version: v0.1 (2019-04-26, Kim)" + "\n")
    out.write("Directory: " + directory + "\n")
    out.write("File name: " + file_name + "\n")
    out.write("Library name: " + library_name + "\n")
    out.write("Reference genome: " + genome_name + "\n")
    out.write("FDR threshold: " + str(fdr_thresh) + "\n")
    out.write("Ratio threshold in entropy filter: " + str(cutoff_ent_filter) + "\n")
    out.write("Number of pseudo-GEMs: " + str(samp_size) + "\n")
    out.write("Started processing domain-based distance test. \n")
    out.write("================================= \n")
        
    ### Read input GEM file ###
    in_gems = []
    in_gems = read_gems(directory, file_name)
    out.write("Finished reading the input GEM file. \n")
    out.write(str(len(in_gems)) + " total GEMs. \n")
    out.write("================================= \n")
    
    ### Read chromosome names and sizes ###
    chrom_dict = read_chroms(chrom_dir, genome_name)

    header = ['GEM_ID', 'GEM_coord', 'GEM_span', 'Frag_number', 'List_of_frag_coord', 'category', 'dist', 'rawpval1', 'adjpval1', 'decis1', 'rawpval2', 'adjpval2', 'decis2']
    with open(out_directory+prefix+'_distTest_master.txt', 'a') as file1:
        file1.write('\t'.join(map(str, header)) + '\n')
    file1.close()

    tot_pass = 0
    for chr_name in chrom_dict.keys():
        ### Subset input GEM files by chromosome ###
        out.write("================================= \n")
        out.write("===== Chromosome is: " + chr_name + ". ===== \n")
        out.write("================================= \n")
        subset_gems = [x for x in in_gems if chr_name == x[0].split(":")[0]]
        out.write("Finished subsetting GEMs. \n")
        out.write(str(len(subset_gems)) + " GEMs in " + chr_name + " (of " + str(len(in_gems)) + " total GEMs). \n")
        if len(subset_gems) < 100:
            out.write("Skipped: less than 100 GEMs. \n")
            out.write("================================= \n")
            continue
            
        ### Initialize output GEM list ###
        out_gems = [] 
        
        ### Create distance bucket ###
        f2f_bucket = []
        fragnum_bucket = []
        for i in range(len(subset_gems)):
            gem_id, chrom, frags, span, fraglen, f2fdist, fragnum, frag_str, coord = extract_info(subset_gems[i])
            f2f_bucket.extend(f2fdist)
            fragnum_bucket.append(fragnum)
        out.write("Maximum number of fragments in a GEM: "+str(max(fragnum_bucket))+"\n")
        
        ### Create pseudoGEMS ###i
        pseudo_file = out_directory + prefix + '_distTest_pseudo_' + chr_name + '.txt'
        pseudo_gem_tot_dist, pseudo_gem_mean_ent = create_pseudo_gems(f2f_bucket, fragnum_bucket, samp_size, pseudo_file)
        out.write("Finished creating "+str(samp_size)+" pseudoGEMs and total distance for each fragment number class. \n")
        del f2f_bucket
        
        ### Calculate raw p-values for subset_gems via distance test ###
        out.write("================================= \n")
        out.write(" 1) Distance test. \n")
        out.write("================================= \n")
        for k in range(len(subset_gems)):
            gem_id, chrom, frags, span, fraglen, f2fdist, fragnum, frag_str, coord = list(extract_info(subset_gems[k]))
            raw_pval = get_raw_pval(f2fdist, pseudo_gem_tot_dist[str(fragnum)], samp_size)
            out_gems.append([gem_id, coord, span, fragnum, frag_str, 'Orig', sum(f2fdist), raw_pval])
        out.write("Finished computing raw p-values for "+str(k+1)+" GEMs. \n")
        
        ### Adjust p-values by fragment number classes in each chromosome ###
        deferred_gems = []
        for fn in sorted(list(set(fragnum_bucket))):
            indx = [i for i,val in enumerate(fragnum_bucket) if val==fn]
            pvals1 = [out_gems[k][7] for k in indx]
            adj_pval1 = get_adj_pval(pvals1, fdr_thresh, 'fdr_bh')
            for k in range(len(indx)):
                out_gems[indx[k]].append(round(adj_pval1[1][k], 4))
                # Subcategories from first distance test #
                if adj_pval1[0][k]==True:
                    out_gems[indx[k]].extend(['PASS', '.', '.', '.'])
                elif adj_pval1[0][k]==False:
                    if fn==2: # GEMs with 2 fragments fail
                        out_gems[indx[k]].extend(['FAIL', '.', '.', '.'])
                    else: # GEMs with > 2 fragments go to entropy filter
                        out_gems[indx[k]].extend(['DEFER', '.', '.', '.'])
                        deferred_gems.append(subset_gems[indx[k]])
        out.write("Finished adjusting p-values. \n")
        out.write("First test PASS: "+str(len([x for x in out_gems if x[9] =='PASS']))+"\n")
        out.write("First test FAIL: "+str(len([x for x in out_gems if x[9] =='FAIL']))+"\n")
        out.write("First test DEFER: "+str(len([x for x in out_gems if x[9] =='DEFER']))+"\n")
        pass1 = len([x for x in out_gems if x[9] =='PASS'])
        
        del pvals1
        del adj_pval1
        
        ### Entropy Filter for those in 'DEFER' category ###
        out.write("================================= \n")
        out.write(" 2) Entropy filter. \n")
        out.write("================================= \n")
        ef_gems = []
        num_frag_list = []
        num_cut_list = []
        for k in range(len(deferred_gems)):
            sub_gems, num_cuts = entropy_filter(deferred_gems[k], pseudo_gem_mean_ent, cutoff_ent_filter)
            num_cut_list.append(num_cuts)
            for i in range(len(sub_gems)): 
                gem_id, chrom, frags, span, fraglen, f2fdist, fragnum, fragstr, coord = list(extract_info([sub_gems[i]]))
                num_frag_list.append(fragnum)
                raw_pval = get_raw_pval(f2fdist, pseudo_gem_tot_dist[str(fragnum)], samp_size)
                ef_gems.append([gem_id, coord, span, fragnum, fragstr, 'EF-'+str(num_cuts), sum(f2fdist),'.','.', 'DEFER', raw_pval])
        out.write("Finished entropy filter for removing doublets. \n")
        out.write("Number of GEMs not split: "+str(len([x for x in num_cut_list if x == 0]))+"\n")
        out.write("Number of GEMs split once: "+str(len([x for x in num_cut_list if x == 1]))+"\n")
        out.write("Number of GEMs split twice: "+str(len([x for x in num_cut_list if x == 2]))+"\n")
        out.write("Number of non-singleton GEMs after cutting: "+str(len([x for x in num_frag_list if x > 1]))+"\n")
        out.write("Number of singleton GEMs after cutting: "+str(len([x for x in num_frag_list if x == 1]))+"\n")
        
        ### Adjust p-values by fragment number classes in each chromosome ###
        out.write("================================= \n")
        out.write(" 3) Distance test.  \n")
        out.write("================================= \n")
        sub_fragnum_bucket = [x[3] for x in ef_gems]
        for fn in sorted(list(set(sub_fragnum_bucket))):
            indx = [i for i,val in enumerate(sub_fragnum_bucket) if val==fn]
            if fn==1:
                for k in range(len(indx)):
                    ef_gems[indx[k]].extend(['.','.'])
            else:
                pvals2 = [ef_gems[k][10] for k in indx]
                adj_pval2 = get_adj_pval(pvals2, fdr_thresh, 'fdr_bh')
                for k in range(len(indx)):
                    ef_gems[indx[k]].append(round(adj_pval2[1][k], 4))
                    if adj_pval2[0][k]==True:
                        ef_gems[indx[k]].append('PASS')
                    elif adj_pval2[0][k]==False:
                        ef_gems[indx[k]].append('FAIL')
        out.write("Finished adjusting p-values. \n")
        out.write("Second test PASS: "+str(len([x for x in ef_gems if x[12] =='PASS']))+"\n")
        out.write("Second test FAIL: "+str(len([x for x in ef_gems if x[12] =='FAIL']))+"\n")
        pass2 = len([x for x in ef_gems if x[12] =='PASS'])
        out.write("================================= \n")
        
        del pvals2
        del adj_pval2
        
        ### Write results ###
        write_master_result(out_gems, out_directory+prefix+'_distTest_master.txt')
        write_master_result(ef_gems, out_directory+prefix+'_distTest_master.txt')
        out.write("Finished writing files. \n")
        
        tot_pass += pass1+pass2
        
        out.write("Total " +  str(pass1+pass2) + " GEMs out of " + str(len(subset_gems))+ " GEMs passed. \n")
        out.write("i.e., " + str(round((pass1+pass2)*100/len(subset_gems), 3)) + "% passed in " +chr_name+ ". \n")
        out.write("================================= \n")
        del out_gems
        del ef_gems
        del deferred_gems
    out.write("================================= \n")
    out.write("======= Summary Statistics ======= \n")
    out.write("Total " +  str(tot_pass) + " GEMs out of " + str(len(in_gems))+ " GEMs passed. \n")
    out.write("i.e., " + str(round((tot_pass)*100/len(in_gems), 3)) + "% passed. \n")
    out.write("================================= \n")
    out.write("DONE. \n")
    
    out.close()
