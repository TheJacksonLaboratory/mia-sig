#!/bin/bash 

## The help message:
function usage
{
    echo -e "usage:
    bash gems2frags.sh --i in_file --o out_file"
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -s | --i )           shift
                                in_file=$1
                                ;;
        -s | --o )           shift
                                out_file=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done


awk '{if (NR!=1) {print $5}}' ${in_file} | tr ';' $'\n' | tr ':' $'\t' | tr '-' $'\t' | cut -f1 -d'(' > ${out_file}

