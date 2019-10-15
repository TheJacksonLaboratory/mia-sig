import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.kde import gaussian_kde
from numpy import linspace
from sys import argv

def modLog(num, denom):
    if num==0 or denom==0:
        return 0
    else:
        return float(format(np.log2(num/denom), '.4f'))

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
    entropy = sum([x*modLog(1,x) for x in prob_list])/np.log2(len(prob_list))
    return float(format(entropy, '.5f'))

def read_file(file_name):
    """
    Read a file.
    """
    with open(file_name) as f:
        in_data = [line.strip().split("\t") for line in f]
    return in_data

def plot_hist(data, binr, tit, xlab, fig_name):
    """
    Plot histogram.
    """
    plt.hist(data, bins=binr)
    plt.title(tit)
    plt.xlabel(xlab)
    plt.ylabel("Frequency")
    plt.savefig(fig_name+'.pdf', dpi=300)
    plt.close()

def plot_dist_byfrag(fragn, frags, category, fig_name):
    """
    Plot neighboring fragment distance distribution categorized by fragment number.
    """
    maxfn = max(fragn)
    minfn = min(fragn)
    if maxfn < 5:
        breakPoints = [range(2,max(fragn)+1)]
        legendLab = ['2-max']
    elif maxfn > 4 and maxfn < 10:
        if minfn < 5:
            breakPoints = [range(2,5), range(5,maxfn+1)]
            legendLab = ['2-4', '5-max']
        else:
            breakPoints = [range(5,maxfn+1)]
            legendLab = ['5-max']
    else:
        if minfn < 5:
            breakPoints = [range(2,5), range(5, 10), range(10,maxfn+1)]
            legendLab = ['2-4', '5-9', '10-max']
        elif minfn > 4 and minfn < 10:
            breakPoints = [range(5, 10), range(10,maxfn+1)]
            legendLab = ['5-9', '10-max']
        else:
            breakPoints = [range(10,maxfn+1)]
            legendLab = ['10-max']
    fdistbyfnum = [[] for i in range(max(fragn)+1)] # distByFragNum[i] contains all distances for i fragments / GEM 
    fdist = []
    for i in range(len(frags)):
        coords = [x.split(':')[1].split('-') for x in frags[i]]
        init_dist = [int(coords[j+1][0])-int(coords[j][1]) for j in range(0, len(coords)-1)]
        dist = [x for x in init_dist if x > 3000]
        if len(dist)>0:
            fdistbyfnum[len(dist)+1].extend(dist)
            fdist.extend(dist)
        
    neigh_distFrag = [[] for i in range(len(breakPoints))]
    for k in range(len(breakPoints)):
        for x in breakPoints[k]:
            neigh_distFrag[k].extend(fdistbyfnum[x])
            
    dist_space = linspace( 3, 8, 100 )
    for y in range(len(breakPoints)):
        plt.plot(dist_space, gaussian_kde(np.log10(neigh_distFrag[y]))(dist_space), linewidth = 4)
    plt.legend(legendLab, title="Fragment #")
    plt.title("F2F distance in "+str(len(fragn)) + " complexes (" + category + ")")
    plt.xlabel("Log10(Fragment-to-fragment distance)")
    plt.ylabel("Relative Density")
    plt.savefig(fig_name+'f2f_by_fnum.pdf', dpi=300)
    plt.close()
    
    dist_space = linspace( 3, 8, 100)
    plt.plot(dist_space, gaussian_kde(np.log10(fdist))(dist_space))
    plt.title("F2F distance in "+str(len(fragn)) + " complexes (" + category + ")")
    plt.xlabel("Log10(Fragment-to-fragment distance)")
    plt.ylabel("Relative Density")
    #plt.show()
    plt.savefig(fig_name+'f2f_all.pdf', dpi=300)
    plt.close()
    
    del fdistbyfnum
    del neigh_distFrag
    return fdist

def plot_ent_byfrag(entbyfrag, category, fig_name):
    """
    Plot normalized Shannon entropy categorized by fragment number.
    """
    maxfn = len(entbyfrag)
    breakPoints = [range(3,6), range(6,10), range(10,maxfn)]
    legendLab = ['3-5', '6-9', '10-max']
    
    ent_cat = [[] for i in range(len(breakPoints))]
    for k in range(len(breakPoints)):
        for x in breakPoints[k]:
            ent_cat[k].extend(entbyfrag[x])
    
    dist_space = linspace( 0, 1, 100 )
    for y in range(len(breakPoints)):
        plt.plot(dist_space, gaussian_kde(ent_cat[y])(dist_space))
    plt.legend(legendLab, title="Fragment #")
    plt.title("Norm. Shannon entropy in "+str(sum([len(x) for x in entbyfrag])) + " complexes (" + category + ")")
    plt.xlabel("Norm. Shannon entropy")
    plt.ylabel("Relative Density")
    plt.savefig(fig_name+'entropy.pdf', dpi=300)
    plt.close()

def plot_entfilt(frag_trace, keys, colcode, fig_name):
    fig = plt.figure(figsize = (30,6))
    medianprops = dict(linewidth = 3, color='gray')
    boxprops = dict(linewidth = 1.5)
    toplot = [np.asarray([]) for i in range(len(keys))]
    i=0
    for key in keys:
        datax = toplot
        datax[i] = np.asarray(frag_trace[key])
        plt.boxplot(datax, widths = 0.6, medianprops = medianprops, boxprops = boxprops)
        i+=1
    plt.xticks([i for i in range(1, 19)], [x+'\n ('+str(len(frag_trace[x])) + ')' for x in keys], fontsize=18)
    plt.yticks(fontsize = 18)
    plt.ylabel('Normalized Shannon entropy', fontsize = 18)
    plt.title('Entropy filter effects', fontsize = 18)
    for ticklabel, tickcolor in zip(plt.gca().get_xticklabels(), colcode):
        ticklabel.set_color(tickcolor)
    plt.savefig(fig_name+'.pdf', dpi=300, bbox_inches="tight")
    plt.close()

def plot_ecdf(x1, x2, dlab, clr1, clr2, tit, xlab, fig_name):
    fig, ax = plt.subplots(figsize=(9, 6))
    
    n1, bins1, patches1 = ax.hist(x1, bins = 1000, density=True, histtype='step', cumulative=True, label=dlab[0], color = clr1, linewidth = 4)
    n2, bins2, patches2 = ax.hist(x2, bins = 1000, density=True, histtype='step', cumulative=True, label=dlab[1], color = clr2, linewidth = 4)
    patches1[0].set_xy(patches1[0].get_xy()[:-1])
    patches2[0].set_xy(patches2[0].get_xy()[:-1])
    ax.legend(loc = 'lower right', fontsize = 18)
    plt.title(tit, fontsize = 18)
    plt.xlabel(xlab, fontsize = 18)
    plt.ylabel("Empirical CDF (ECDF)", fontsize = 18)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.title(tit, fontsize = 18)
    plt.savefig(fig_name+'.pdf', dpi=300, bbox_inches="tight")
    plt.close()    

if __name__ == '__main__':
    ### Set parameters ###
    lib_name = argv[1] # ex: 'GSM3347523'
    in_dir = argv[2]
    nf = argv[3] # null file
    mf = argv[4] # master file

    # change to in_dir 
    os.chdir(in_dir)

    ### Read files ###
    master_data = read_file(mf)
    pseudo_data = read_file(nf)

    pass_data = [x for x in master_data if x[9]=='PASS' or x[12]=='PASS']
    fail_data = [x for x in master_data if x[9]=='FAIL' or x[12]=='FAIL']

    #### Part 1.1) Plot fragment number and fragment-to-fragment distance distributions ####

    # Pass category
    pass_fragn = [int(x[3]) for x in pass_data]
    plot_hist(pass_fragn, range(2,40,1), "Number of fragments in "+str(len(pass_fragn)) + " complexes (pass)", "Number of fragments", in_dir+lib_name+"_pass_fragnum")
    pass_frags = [x[4].split(';') for x in pass_data]
    pass_fdist = plot_dist_byfrag(pass_fragn, pass_frags, "pass", in_dir+lib_name+"_pass_")
    del pass_fragn
    del pass_frags

    # Fail category
    fail_fragn = [int(x[3]) for x in fail_data]
    plot_hist(fail_fragn, range(2,40,1), "Number of fragments in "+str(len(fail_fragn)) + " complexes (fail)", "Number of fragments", in_dir+lib_name+"_fail_fragnum")
    fail_frags = [x[4].split(';') for x in fail_data]
    fail_fdist = plot_dist_byfrag(fail_fragn, fail_frags, "fail", in_dir+lib_name+"_fail_")
    del fail_fragn
    del fail_frags

    # Original category
    orig_fragn = [int(x[3]) for x in master_data if x[5]=='Orig']
    plot_hist(orig_fragn, range(2,40,1), "Number of fragments in "+str(len(orig_fragn)) + " complexes (orig)", "Number of fragments", in_dir+lib_name+"_orig_fragnum")
    orig_frags = [x[4].split(';') for x in master_data if x[5]=='Orig']
    orig_fdist = plot_dist_byfrag(orig_fragn, orig_frags, "orig", in_dir+lib_name+"_orig_")
    del orig_fragn
    del orig_frags

    plot_ecdf(pass_fdist, orig_fdist, ['Significant', 'Original'], 'firebrick', 'black', 'Original vs. significant', 'Fragment-to-fragment distances', in_dir+lib_name+"_original_vs_pass_distance_ecdf")

    # Pseudo category
    pseudo_fragn = [int(x[0]) for x in pseudo_data]
    pseudo_frags = []
    for x in pseudo_data:
        init = ['chr2L:0-100']
        end = 100
        dists = x[1].split(';')
        for j in range(len(dists)):
            start = end + int(dists[j])
            init.append('chr2L:'+str(start)+'-'+str(start+100))
            end = start + 100
        pseudo_frags.append(init)
    pseudo_fdist = plot_dist_byfrag(pseudo_fragn, pseudo_frags, "pseudo null", in_dir+lib_name+"_pseudonull_")
    del pseudo_fragn
    del pseudo_frags

    #### Part 1.2) Shannon entropy filter on deferred complexes after first distance test   ####

    # After first distance test: effects of entropy filter
    deferred_data = [x for x in master_data if x[9]=='DEFER']

    deferred_bybc = {}
    for x in deferred_data:
        bc = x[0].split('-sub')[0]
        if bc not in deferred_bybc.keys():
            deferred_bybc[bc]= [x]
        else:
            deferred_bybc[bc].append(x)

    frag_trace = {}
    for key, value in deferred_bybc.items():
        orig_coord = [x.split(':')[1] for x in value[0][4].split(';')]
        dist = [int(orig_coord[i+1].split('-')[0])-int(orig_coord[i].split('-')[1]) for i in range(len(orig_coord)-1)]
        ent = norm_shannon_ent(dist)
        fncomb = ','.join(str(e) for e in sorted([int(x[3]) for x in value[1:]]))
        if fncomb not in frag_trace.keys():
            frag_trace[fncomb] = [ent]
        else: 
            frag_trace[fncomb].append(ent)

    keys1 = ['3', '1,2', '4', '1,3', '2,2', '1,1,2', '5', '1,4', '2,3', '1,1,3', '1,2,2', '6', '1,5', '2,4', '3,3', '1,1,4', '1,2,3', '2,2,2']
    color1 = ['red', 'red', 'orange','orange','orange','orange','green','green','green','green','green','blue','blue','blue','blue','blue','blue','blue']

    keys2 = ['7', '1,6', '2,5', '3,4', '1,1,5', '1,2,4', '1,3,3', '2,2,3', '8', '1,7', '2,6', '3,5', '4,4', '1,1,6', '1,2,5', '1,3,4', '2,2,4', '2,3,3']
    color2 = ['navy','navy','navy','navy','navy','navy','navy','navy','purple','purple','purple','purple','purple','purple','purple','purple','purple','purple']

    plot_entfilt(frag_trace, keys1, color1, in_dir+lib_name+'_entropy_filt_boxplot_set1')
    plot_entfilt(frag_trace, keys2, color2, in_dir+lib_name+'_entropy_filt_boxplot_set2')



