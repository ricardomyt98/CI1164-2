#!/bin/bash
# Command to use the count registers inside the processor.
# sudo modprobe msr
modprobe msr

# Cleaning up files.
rm -rf ./likwidPerformance/*

array=(32 50 64)
# array=(32 50 64 100 128 200 256 300 400 512 1000 1024 2000 2048 3000 4000 4096 5000 10000)

# [L3].
for nx_ny in ${array[*]}
do
    printf "%s " $nx_ny >> ./likwidPerformance/L3.dat
    likwid-perfctr -f -g L3 -C 0 ./pdeSolver -nx $nx_ny -ny $nx_ny -i 10 -o arquivo_saida | grep "load bandwidth" | grep -o -P "[0-9]+\.[0-9]+" >> ./likwidPerformance/L3.dat
done

# [L2CACHE].
for nx_ny in ${array[*]}
do
    printf "%s " $nx_ny >> ./likwidPerformance/L2_CACHE.dat
    likwid-perfctr -f -g L2CACHE -C 0 ./pdeSolver -nx $nx_ny -ny $nx_ny -i 10 -o arquivo_saida | grep "L2 miss ratio" | grep -o -P "[0-9]+\.[0-9]+" >> ./likwidPerformance/L2_CACHE.dat
done

# [FLOPS_DP].
for nx_ny in ${array[*]}
do
    printf "%s " $nx_ny >> ./likwidPerformance/DP.dat
    likwid-perfctr -f -g FLOPS_DP -C 0 ./pdeSolver -nx $nx_ny -ny $nx_ny -i 10 -o arquivo_saida | grep "DP \[MFLOP\/s\]" | head -1 | grep -o -P "[0-9]+\.[0-9]+" >> ./likwidPerformance/DP.dat
    printf "%s " $nx_ny >> ./likwidPerformance/AVX_DP.dat
    likwid-perfctr -f -g FLOPS_DP -C 0 ./pdeSolver -nx $nx_ny -ny $nx_ny -i 10 -o arquivo_saida | grep "AVX DP \[MFLOP\/s\]" | grep -o -P "[0-9]+" >> ./likwidPerformance/AVX_DP.dat
done
