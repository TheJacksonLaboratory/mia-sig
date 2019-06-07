#!/bin/bash 

## The help message:
function usage
{
    echo -e "usage:
    bash gems2cov.sh --conf conf_file --dir directory --i in_file --lib libid --r ref --bs binsize"
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -s | --conf )         shift
                                conf_file=$1
                                ;;
        -s | --dir )           shift
                                directory=$1
                                ;;
        -s | --i )           shift
                                in_file=$1
                                ;;
        -s | --lib )           shift
                                libid=$1
                                ;;
        -s | --r )           shift
                                ref=$1
                                ;;
        -s | --bs )           shift
                                binsize=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

module load gcc/4.9.2
module load bedtools/2.27.0

cd ${directory}
source ${conf_file}

# convert gems to fragments
if [ ! -f "${libid}.frags.bed" ]
then
    bash ${miasig_dir}gems2frags.sh --i ${in_file} --o ${libid}.frags.bed
fi

# sort
if [ ! -f "${libid}.frags.sorted.bed" ]
then
    sort -k1,1V -k2,2n ${libid}.frags.bed > ${libid}.frags.sorted.bed
fi

# generate coverage
if [ ! -f "${libid}.bedgraph" ]
then
    bedtools genomecov -bg -i ${libid}.frags.sorted.bed -g ${ref}.chrom.sizes > ${libid}.bedgraph
fi

# create binned genome if it doesn't exist
if [ ! -f "${ref}.${binsize}bp.windows.bed" ]
then
    bedtools makewindows -g ${ref}.chrom.sizes -w ${binsize} > ${ref}.${binsize}bp.windows.bed
fi

# bin coverage by binsize
if [ ! -f "${libid}_binned_${binsize}bp.bedgraph" ]
then
    bedtools map -a ${ref}.${binsize}bp.windows.bed -b ${libid}.bedgraph -c 4 -o mean -null 0 > ${libid}_binned_${binsize}bp.bedgraph
fi

# separate by chromosome
for chrom in $(cut -f1 ${ref}.chrom.sizes); do
    awk '{ if ($1=="'${chrom}'") {print}}' ${libid}_binned_${binsize}bp.bedgraph > ${libid}_binned_${binsize}bp_${chrom}.bedgraph; 
done

