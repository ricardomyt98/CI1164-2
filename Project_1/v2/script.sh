#!/bin/bash
# Command to use the count registers inside the processor.
modprobe msr
export PATH=/home/soft/likwid/bin:/home/soft/likwid/sbin:$PATH
export LD_LIBRARY_PATH=/home/soft/likwid/lib:$LD_LIBRARY_PATH

# Set the processor frequency.
echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

# Cleaning up files.
rm -rf ./likwidPerformance/*

array=(32 50 64)
# array=(32 50 64 100 128 200 256 300 400 512 1000 1024 2000 2048 3000 4000 4096 5000 10000)

# Get the processor topology.
# likwid-topology -g -c

# [L3].
for nx_ny in ${array[*]}
do
    printf "%s " $nx_ny >> ./likwidPerformance/L3.dat
    likwid-perfctr -m -f -g L3 -C 0 ./pdeSolver -nx $nx_ny -ny $nx_ny -i 10 -o arquivo_saida | grep "L3 bandwidth" | grep -o -P "[0-9]+\.[0-9]+" >> ./likwidPerformance/L3.dat
done

# [L2CACHE].
for nx_ny in ${array[*]}
do
    printf "%s " $nx_ny >> ./likwidPerformance/L2_CACHE.dat
    likwid-perfctr -m -f -g L2CACHE -C 0 ./pdeSolver -nx $nx_ny -ny $nx_ny -i 10 -o arquivo_saida | grep "L2 miss ratio" | grep -o -P "[0-9]+\.[0-9]+" >> ./likwidPerformance/L2_CACHE.dat
done

# [FLOPS_DP].
for nx_ny in ${array[*]}
do
    printf "%s " $nx_ny >> ./likwidPerformance/DP.dat
    likwid-perfctr -m -f -g FLOPS_DP -C 0 ./pdeSolver -nx $nx_ny -ny $nx_ny -i 10 -o arquivo_saida | grep "DP MFLOP\/s" | head -1 | grep -o -P "[0-9]+\.[0-9]+" >> ./likwidPerformance/DP.dat
    printf "%s " $nx_ny >> ./likwidPerformance/AVX_DP.dat
    likwid-perfctr -m -f -g FLOPS_DP -C 0 ./pdeSolver -nx $nx_ny -ny $nx_ny -i 10 -o arquivo_saida | grep "AVX DP MFLOP\/s" | grep -o -P "[0-9]+" >> ./likwidPerformance/AVX_DP.dat
done

# Set the processor frequency.
echo "powersave" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor