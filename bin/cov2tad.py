import numpy as np
import os
from numpy import linspace
from numpy import diff
from scipy.signal import argrelextrema
import scipy
from scipy import signal
import pywt
from itertools import chain
from sys import argv

def read_bg(directory, file_name):
    """
    Read a bedgraph file representing potential TAD track.
    NOTE: bedgraph must be binned at a pre-defined constant window.
    Args:
       directory (str): directory of the file location (ex: '/Users/kimm/')
       file_name (str): name of the file (ex: 'SHG0008H.bedgraph')
    Returns:
       cov_dict (dict): dictionary of coverage
    """
    cov_dict = {}
    with open(directory + file_name) as f:
        in_cov = [line.strip().split("\t") for line in f]
    chrom_list = list(set([x[0] for x in in_cov]))
    for chrom in chrom_list:
        tmp_cov = [int(float(x[3])) for x in in_cov if x[0]==chrom]
        cov_dict[chrom] = tmp_cov
    return cov_dict

def read_chroms(chrom_file):
    """
    Read a tab-delimited text file with list of chromosomes and their sizes.
    Args:
       chrom_file (str): location of chrom sizes file (ex: '/Users/kimm/dm3.chrom.sizes')
    Returns:
       chrom_dict (dictionary): tab-separated lists of lists
    """
    chrom_dict = {}
    with open(chrom_file) as f:
        for line in f:
            tmp_list = line.strip().split("\t")
            chrom_dict[tmp_list[0]] = int(tmp_list[1])
    return chrom_dict

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

def takederiv(y, dx):
    """
    Take the derivative of a signal y and return dy/dx.
    Args:
       y (list): list of numbers representing time-series signal
       dx (int): stepsize to compute the slope
    Returns:
       dy (list): list of length (len(y)-1)
    """
    dy = list(diff(y)/dx)
    return dy

def waveletdec(s, wavelet, lvl, coord):
    """
    Decompose .
    Args:
       s (list): list of numbers representing time-series signal
       wavelet (string): name of the wavelet; ex: 'sym2', 'bior1.3', etc.
       lvl (int): level of decomposition; length of coeffs will be lvl+1
       coord (str): genomic coordinates of the signal s
    Returns:
       coeffs (list of list): first list being approximation and the rest details
       Also plots approximation
    """
    coeffs = pywt.wavedec(s, wavelet, level = lvl)
    return coeffs

def leftlimit(minpoint, ds, tau):
    """
    Find left limit.
    Args:
       minpoint (int): x coordinate of the local minimum
       ds (list): list of numbers representing derivative of the time-series signal s
       tau (int): threshold to determine a slope 'too sharp'
    Returns:
       minpoint (int): x coordinate of the left-modified local minimum
    """
    slope = ds[minpoint-1]
    while slope > tau and slope < 0 and minpoint > 0:
        minpoint -= 1
        slope = ds[minpoint-1]
    return minpoint

def rightlimit(minpoint, ds, tau):
    """
    Find right limit.
    Args:
       minpoint (int): x coordinate of the local minimum
       ds (list): list of numbers representing derivative of the time-series signal s
       tau (int): threshold to determine a slope 'too sharp'
    Returns:
       minpoint (int): x coordinate of the right-modified local minimum
    """
    slope = ds[minpoint]
    while slope < tau and slope > 0 and minpoint > 0 and minpoint < len(ds)-1:
        minpoint += 1
        slope = ds[minpoint]
    return minpoint

def adjustmin(minpoint, ds, tau1, tau2):
    """
    Find right limit.
    Args:
       minpoint (int): x coordinate of the local minimum
       ds (list): list of numbers representing derivative of the time-series signal s
       tau1 (int): threshold to determine a slope 'too sharp' to the left
       tau2 (int): threshold to determine a slope 'too sharp' to the right
    Returns:
       interval (list): left- and right- adjusted interval
    """
    if minpoint == len(ds):
        return [minpoint, minpoint]
    else:
        left = leftlimit(minpoint, ds, tau1)
        right = rightlimit(minpoint, ds, tau2)
        interval = [left, right]
    return interval

def get_locext(s):
    """
    Obtain local minima and maxima.
    Args:
       s (list): list of numbers representing time-series signal
    Returns:
       locmax (list): indices (i.e., x-coordinates) of local maxima
       locmin (list): indices (i.e., x-coordinates) of local minima
    """
    locmax = list(argrelextrema(s, np.greater, order = 1)[0])
    locmin = list(argrelextrema(s, np.less_equal, order = 1)[0])
    locmin.append(0)
    locmin.append(len(s)-1)
    locmin = sorted(list(set(locmin)))
    return locmax, locmin

def extreme_min(s, locmax, locmin, thr_rel, thr_abs):
    """
    Obtain local minima and maxima.
    Args:
       s (list): list of numbers representing time-series signal
       locmax (list): indices (i.e., x-coordinates) of local maxima
       locmin (list): indices (i.e., x-coordinates) of local minima
       thr_rel (float): relative threshold controlling 'how extreme' local minima need to be
       thr_abs (int): absolute threshold controlling 'how extreme' local minima need to be
    Returns:
       locmin_ext (list): indices (i.e., x-coordinates) of extreme local minima
    """
    locmin2 = []
    thresh = int(np.median(s))
    for item in locmax:
        indx = np.searchsorted(locmin, item)
        max_s = int(s[item])
        indx_l = locmin[indx-1]
        min_l = int(s[indx_l])
        indx_r = locmin[indx]
        min_r = int(s[indx_r])
        if min_l/max_s < thr_rel or min_l < thr_abs:
            locmin2.append(indx_l)
        if min_r/max_s < thr_rel or min_r < thr_abs:
            locmin2.append(indx_r)
    locmin_ext = sorted(list(set(locmin2)))
    return locmin_ext

def segment(locmax, locmin):
    """
    Segment the signal by selecting local maximum bounded by neighboring local minima.
    Args:
       locmax (list): indices (i.e., x-coordinates) of local maxima
       locmin (list): indices (i.e., x-coordinates) of local minima
    Returns:
       region (list): indices (i.e., x-coordinates) of extreme local minima
    """
    region = []
    loc_pre = 0
    for x in locmax:
        loc = np.searchsorted(locmin, x)
        if loc < len(locmin) and loc != loc_pre:
            region.append([locmin[loc-1],locmin[loc]])
            loc_pre = loc
    return region

def get_regions(s):
    """
    Segment the signal by selecting local maximum bounded by neighboring local minima.
    Args:
       s (list): list of numbers representing time-series signal
       locmin (list): indices (i.e., x-coordinates) of local minima
    Returns:
       region (list): indices (i.e., x-coordinates) of extreme local minima
    """
    locmax, locmin = get_locext(s)
    ds = takederiv(s, 1)
    tau1 = int(np.median([x for x in ds if x < 0]))
    tau2 = int(np.median([x for x in ds if x > 0]))
    upperbnd = int(np.median(s))
    lowerbnd = int(np.median([s[i] for i in locmin]))
    locmin_ext = extreme_min(s, locmax, locmin, 0.5, lowerbnd)
    adjlocmin = sorted(list(set(list(chain(*[adjustmin(x, ds, tau1, tau2) for x in locmin_ext])))))
    refined_locmax = [x for x in locmax if s[x] > upperbnd]
    regions = segment(refined_locmax, adjlocmin)
    return regions

def write_tads(final_tads, file_name):
    """
    Segment the signal by selecting local maximum bounded by neighboring local minima.
    Args:
       s (list): list of numbers representing time-series signal
       locmin (list): indices (i.e., x-coordinates) of local minima
    Returns:
       region (list): indices (i.e., x-coordinates) of extreme local minima
    """
    jbox = [x*2 for x in final_tads]
    with open(file_name, 'a') as file1:
        file1.write('\t'.join(map(str, ['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'color', 'comment'])) + '\n')
        for i in range(len(jbox)):
            file1.write('\t'.join(map(str, jbox[i])) + '\t' + '0,0,0' + '\t' + 'TADs' + '\n')
    file1.close()

def write_beds(final_tads, file_name):
    """ 
    Write out fragments in short format (regardless PLEC, PLISRS, all pairs).
    Args: 
       pair_list (list): list of pairs of fragments in short or bedpe format
       out_name (string): output file name
    Returns:
       None
    """
    with open(file_name, 'a') as file1:
        for i in range(len(final_tads)):
            file1.write('\t'.join(map(str, final_tads[i])) + '\n')
    file1.close()

if __name__ == '__main__':
    directory = argv[1]
    file_name = argv[2]
    cov_dict = read_bg(directory, file_name)
    bin_size = int(argv[3])
    species = argv[4]
    chrom_file = argv[5]
    tadsize_ref = tadsize_chart(species)
    chrom_dict = read_chroms(chrom_file)
    prefix = file_name.split('.bedgraph')[0]
    
    #### Log file ####
    out = open(directory + prefix + "_tadcalling_logFile.txt", "a")
    
    out.write("Program version: v0.1 (2019-04-30, Kim)" + "\n")
    out.write("Directory: " + directory + "\n")
    out.write("File name: " + file_name + "\n")
    out.write("Species genome: " + species + "\n")
    out.write("Started calling TADs. \n")
    out.write("================================= \n")
    
    final_tads = []
    for chr_name in chrom_dict.keys():
        out.write("================================= \n")
        out.write("===== Chromosome is: " + chr_name + ". ===== \n")
        out.write("================================= \n")
        chr_size = len(cov_dict[chr_name])
        win_size = min(int(tadsize_ref[2]/bin_size*10), chr_size-1)
        out.write("Window size: " + str(win_size) + "\n")
        start = 0
        end = start + win_size
        coord = chr_name+":"+str(start)+"-"+str(end)
        tad_idx = []
        while end < chr_size+1:
            sig = cov_dict[chr_name][start:end]
            if len(np.nonzero(sig)[0]) > 0:
                coef = waveletdec(sig, 'bior1.1', 3, coord)
                sig_appr = coef[0]
                regions = get_regions(sig_appr)
                ratio = len(sig)/len(sig_appr)
                tad_idx.extend([[x[0]*ratio+start, x[1]*ratio+start] for x in regions])
            if end == chr_size:
                break
            start = int(regions[-1][1]*ratio)+start-1
            end = min(start + win_size, chr_size)
        tad_bed = [[chr_name, int(round(x[0]*bin_size, -3)), int(round(x[1]*bin_size, -3))] for x in tad_idx]
        final_tads.extend(tad_bed)
        out.write("Number of TADs: " + str(len(tad_bed)) + "\n")
        if len(tad_bed)>0:
            tad_sizes = [x[2]-x[1] for x in tad_bed]
            out.write("Minimum TAD size: " + str(min(tad_sizes)) + "\n")
            out.write("Maximum TAD size: " + str(max(tad_sizes)) + "\n")
            out.write("Mean TAD sizes: " + str(int(np.mean(tad_sizes))) + "\n")
            out.write("Median TAD sizes: " + str(int(np.median(tad_sizes))) + "\n")
    write_tads(final_tads, directory + prefix + '_TAD_JBox.txt')
    for i in range(len(final_tads)):
        final_tads[i].extend(['T'+str(i+1)])
    write_beds(final_tads, directory + prefix + '_TAD.bed')
    out.write("================================= \n")
    out.write("Total number of TADs: " + str(len(final_tads)) + "\n")
    final_tad_sizes = [x[2]-x[1] for x in final_tads]
    out.write("Mean TAD sizes: " + str(int(np.mean(final_tad_sizes))) + "\n")
    out.write("Median TAD sizes: " + str(int(np.median(final_tad_sizes))) + "\n")
    out.write("DONE.")
    out.close()
