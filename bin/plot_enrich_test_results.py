import os
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

def read_file(file_name):
    """
    Read a file.
    """
    with open(file_name) as f:
        in_data = [line.strip().split("\t") for line in f]
    return in_data

def plot_ecdf(x1, x2, dlab, clr1, clr2, tit, xlab, fig_name):
    fig, ax = plt.subplots(figsize=(6, 6))
    
    n1, bins1, patches1 = ax.hist(x1, bins = 1000, density=True, histtype='step', cumulative=True, label=dlab[0], color = clr1, linewidth = 4)
    n2, bins2, patches2 = ax.hist(x2, bins = 1000, density=True, histtype='step', cumulative=True, label=dlab[1], color = clr2, linewidth = 4)
    patches1[0].set_xy(patches1[0].get_xy()[:-1])
    patches2[0].set_xy(patches2[0].get_xy()[:-1])
    ax.legend(loc = 'lower right', fontsize = 18)
    plt.title(tit, fontsize = 18)
    plt.xlabel(xlab, fontsize = 18)
    plt.ylabel("Empirical Cum. Dist. Func. (ECDF)", fontsize = 18)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.title(tit, fontsize = 18)
    plt.savefig(fig_name+'.pdf', dpi=300, bbox_inches="tight")
    plt.close()  

def plot_2hist(x1, x2, bin_lims, lab1, lab2, clr1, clr2, tit, xlab, fig_name):
    bin_centers = 0.5*(bin_lims[:-1]+bin_lims[1:])
    bin_widths = bin_lims[1:]-bin_lims[:-1]
    hist1, _ = np.histogram(x1, bins=bin_lims)
    hist2, _ = np.histogram(x2, bins=bin_lims)
    
    ##normalizing
    hist1b = hist1/np.max(hist1)
    hist2b = hist2/np.max(hist2)

    fig, (ax2) = plt.subplots(nrows = 1, ncols = 1, figsize=(8, 6))

    ax2.bar(bin_centers, hist1b, width = bin_widths, align = 'center', label = lab1, color = clr1, alpha = 0.5)
    ax2.bar(bin_centers, hist2b, width = bin_widths, align = 'center', label = lab2, color = clr2, alpha = 0.5)
    ax2.legend(loc = 'upper right', fontsize = 18)    
    plt.title(tit, fontsize = 18)
    plt.xlabel(xlab, fontsize = 18)
    plt.ylabel("Relative Proportion", fontsize = 18)
    #plt.show()
    plt.savefig(fig_name+'.pdf', dpi=300)

if __name__ == '__main__':
    ### Set parameters ###
    lib_name = argv[1] # ex: 'GSM3347525NR'
    in_dir = argv[2] # '/Users/kimm/Documents/miasig_analysis/GSM3347525NR_EnrichTest_FDR_0.1/'
    #out_dir = '/Users/kimm/Documents/miasig_analysis/raw_figures_RNAPII_Oct9/'
    nf = argv[3] # ex: 'GSM3347525NR_chr2L_FDR_0.1_pseudoGEM_1000_enrichTest_null.txt'
    mf = argv[4] # ex: 'GSM3347525NR_chr2L_FDR_0.1_pseudoGEM_1000_enrichTest_master.txt'

    # change to in_dir 
    os.chdir(in_dir)

    null = read_file(nf)
    obs_master = read_file(mf)
    #null = read_file('GSM3347525NR_chr2L_FDR_0.1_pseudoGEM_1000_enrichTest_null.txt')
    #obs_master = read_file('GSM3347525NR_chr2L_FDR_0.1_pseudoGEM_1000_enrichTest_master.txt')

    obs = [float(x[6]) for x in obs_master if x[0]!='GEM_ID']

    tot_null = []
    for i in range(len(null)):
        tot_null.extend([float(x) for x in null[i][0].split(',')])

    plot_ecdf(tot_null, obs, ['Empirical null', 'Observed'], 'gray', 'brown', 'Empirical null vs. observed', 'Enrichment score', in_dir+lib_name+'_enrichment_score_null_obs_ecdf')

    plot_2hist(tot_null, obs, np.linspace(0,100,101), 'Null', 'Observed', 'black', 'firebrick', 'Empirical null vs. observed', 'Enrichment score', in_dir+lib_name+'_enrichment_score_null_obs_hist_max100')


