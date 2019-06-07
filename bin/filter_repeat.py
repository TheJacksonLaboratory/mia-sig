import numpy as np
from sys import argv
import os
import pybedtools


if __name__ == '__main__':

    ### Set directory and input file name ###
    repeat_directory = argv[1]
    input_file = argv[2]
    output_file = argv[3]
    chr_name = argv[4]    

    #### Repeat region ####
    rptregion = pybedtools.BedTool(repeat_directory+chr_name+'.bed')

    with open(output_file, 'a') as f1:
        with open(input_file) as f2:
            #f1.write(f2.readline() + '\n')
            next(f2) # skip the first line, since it's a header
            for line in f2:
                origline = line.strip().split('\t')
                origfrags = origline[4].split(';')
                if origfrags[0].split(':')[0]==chr_name:
                    #gemid = line.strip().split('\t')[0].split('|')[1]
                    filtfrags = []
                    counter = 0
                    for k in range(len(origfrags)):
                        frag_list = origfrags[k].replace('-', ':').replace('(', ':').split(':')
                        frag_list[2] = str(int(frag_list[2])-500)
                        frag_region = pybedtools.BedTool(' '.join(frag_list[0:3]), from_string=True)
                        intersected = rptregion.intersect(frag_region)
                        if intersected.count()==0: # not falling in the repeat region
                            filtfrags.append(origfrags[k].split('(')[0]+"("+str(counter)+")")
                            counter += 1
                    #finalfrags = ';'.join(filtfrags)
                    #finalline = [origline[0], filtfrags[0].split('-')[0]+filtfrags[-1].split('-')[1].split('(')[0], counter]
                    #finalline.append(finalfrags)
                    if(counter>1):
                        finalfrags = ';'.join(filtfrags)
                        start = filtfrags[0].split(':')[1].split('-')[0]
                        end = filtfrags[-1].split('-')[1].split('(')[0]
                        crd = chr_name + ':' + start + '-' + end
                        #crd = filtfrags[0].split('-')[0] + '-'
                        #crd = crd + filtfrags[-1].split('-')[1].split('(')[0]
                        finalline = [origline[0], crd, int(end)-int(start), counter]
                        finalline.append(finalfrags)
                        f1.write('\t'.join(map(str, finalline)) + '\n')
        f2.close()
    f1.close()

