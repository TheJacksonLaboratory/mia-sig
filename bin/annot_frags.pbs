#!/bin/bash 
# Set parameters for the cluster
#PBS -m ae
#PBS -l nodes=1:ppn=4
#PBS -l walltime=3:00:00
#PBS -l mem=16GB
#PBS -l vmem=16GB
#PBS -j oe

## The help message:

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -s | --conf )         shift
                                conf_file=$1
                                ;;
        -s | --dir )         shift
                                directory=$1
                                ;;
        -s | --in_f )           shift
                                in_file=$1
                                ;;
        -s | --tad )         shift
                                tad_file=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done


# Load Python 

#module load python/3.6.0
module load python/3.6.6 
module load gcc/4.9.2
module load bedtools/2.27.0

# Change directory

cd ${PBS_O_WORKDIR}
source ${conf_file}

#### Arguments for python code: ####

python ${miasig_dir}annot_frags.py ${directory} ${in_file} ${tad_file}


