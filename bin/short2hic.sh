#!/bin/bash 

## The help message:
function usage
{
    echo -e "usage:
    bash short2hic.sh --conf conf_file --dir directory --fn file_name --gf genome_file"
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -s | --conf )         shift
                                conf_file=$1
                                ;;
        -s | --dir )         shift
                                directory=$1
                                ;;
        -s | --fn )           shift
                                file_name=$1
                                ;;
        -s | --gf)            shift
                                genome_file=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

cd ${directory}
source ${conf_file}

#juicer1="/opt/compsci/juicer/1.7.5/CPU/common/"
#juicer2="/opt/compsci/juicer/1.7.5/"

module load java

sort -k2,2d -k6,6d ${file_name} > "${file_name//.short/}.sorted.txt"

sorted_file_name="${file_name//.short/}.sorted.txt"

uniq ${sorted_file_name} > "${sorted_file_name//.sorted.txt/}.uniq.sorted.txt"

uniq_file_name="${sorted_file_name//.sorted.txt/}.uniq.sorted.txt"

# Drosophila: use 2.5M, 1M, 500K, 250K, 100K, 50K, 25K, 10K, 5K, 1K
if [[ ${genome_file} =~ "dm3" ]]
then
    java -Xmx8g -jar ${juicer_loc}juicer.jar pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 "${uniq_file_name}" "${file_name//.short/}.hic" ${genome_file};
fi

# Mouse: use 2.5M, 1M, 500K, 250K, 100K, 50K, 25K
if [[ ${genome_file} =~ "mm9" ]] || [[ ${genome_file} =~ "mm10" ]]
then
    java -Xmx8g -jar ${juicer_loc}juicer.jar pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000 "${uniq_file_name}" "${file_name//.short/}.hic" ${genome_file};
fi

# Human: use 2.5M, 1M, 500K, 250K, 100K, 50K, 25K
if [[ ${genome_file} =~ "hg19" ]] || [[ ${genome_file} =~ "hg38" ]]
then
    java -Xmx8g -jar ${juicer_loc}juicer.jar pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000 "${uniq_file_name}" "${file_name//.short/}.hic" ${genome_file};
fi

# remove irrelevant files

rm ${sorted_file_name}
rm ${uniq_file_name}
