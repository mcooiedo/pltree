#!/bin/bash
# 2 4 8 16 24 32 48
thread_num=(0 1 2 4 8 16 24 28)
for j in 1
do
    rm -f /mnt/pmem0/template.data
    numactl --cpunodebind=0 --membind=0 ./cmake-build-debug/pltree \
        --keys_file=$1 \
        --keys_file_type=$2 \
        --keys_type=$3 \
        --init_num_keys=1000000 \
        --workload_keys=$4 \
        --total_num_keys=$5 \
        --operation=$6 \
        --insert_frac=${9} \
        --lookup_distribution=uniform \
        --theta=0.99 \
        --using_epoch=$7 \
        --thread_num=${thread_num[$j]} \
        --index=$8 \
        --random_shuffle \
        --sort_bulkload=$8 \
        --positive_search=${10}
done
