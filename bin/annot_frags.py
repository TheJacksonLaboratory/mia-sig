import pybedtools
from sys import argv

if __name__ == '__main__':
    
    ### Set directory and input file name ###
    directory = argv[1]
    in_file = argv[2]
    tad_file = argv[3]
    out_file = in_file.split('.txt')[0]+'_annot.txt'
    tads = pybedtools.BedTool(directory+tad_file)
    
    with open(directory+out_file, 'a') as f1:
        with open(directory+in_file) as f2:
            header = f2.readline()
            header = header.rstrip('\n') + '\t' + 'TAD_annot'
            f1.write(header + '\n')
            for line in f2:
                origline = line.strip().split('\t')
                origfrags = origline[4].split(';')
                tad_annot = []
                for k in range(len(origfrags)):
                    frag_list = origfrags[k].replace('-', ':').replace('(', ':').split(':')
                    frag_region = pybedtools.BedTool(' '.join(frag_list[0:3]), from_string=True)
                    intersected = tads.intersect(frag_region, wa=True)
                    if intersected.count() > 0:
                        tad_annot.append(intersected[0][3])
                    else:
                        tad_annot.append('-')
                origline.append(';'.join(tad_annot))
                f1.write('\t'.join(map(str, origline)) + '\n')
        f2.close()
    f1.close()
