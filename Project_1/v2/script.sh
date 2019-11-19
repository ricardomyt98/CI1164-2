#!/bin/bash
# Command to use the count registers inside the processor.
modprobe msr
export PATH=/home/soft/likwid/bin:/home/soft/likwid/sbin:$PATH
export LD_LIBRARY_PATH=/home/soft/likwid/lib:$LD_LIBRARY_PATH

# Set the processor frequency.
echo "performance" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

# Cleaning up files.
rm -rf ./likwidPerformance/*

# array=(32 50 64 100 128 200 256 300 400 512 1000 1024 2000 2048 3000 4000 4096 5000)
array=(32 50 64 100 128 200 256 300 400 512 1000 1024 2000 2048 3000)

# Get the processor topology.
# likwid-topology -g -c

# [L3].
for nx_ny in ${array[*]}
do
    printf "%s " $nx_ny >> ./likwidPerformance/L3_Gauss_Seidel_Likwid_Performance.dat
    printf "%s " $nx_ny >> ./likwidPerformance/L3_L2_Norm_Likwid_Performance.dat
    likwid-perfctr -m -f -g L3 -C 0 ./pdeSolver -nx $nx_ny -ny $nx_ny -i 10 -o arquivo_saida | grep "L3 bandwidth" | head -1 | grep -o -P "[0-9]+\.[0-9]+" >> ./likwidPerformance/L3_Gauss_Seidel_Likwid_Performance.dat
    likwid-perfctr -m -f -g L3 -C 0 ./pdeSolver -nx $nx_ny -ny $nx_ny -i 10 -o arquivo_saida | grep "L3 bandwidth" | tail -1 | grep -o -P "[0-9]+\.[0-9]+" >> ./likwidPerformance/L3_L2_Norm_Likwid_Performance.dat
done

# [L2CACHE].
for nx_ny in ${array[*]}
do
    printf "%s " $nx_ny >> ./likwidPerformance/L2_CACHE_Gauss_Seidel_Likwid_Performance.dat
    printf "%s " $nx_ny >> ./likwidPerformance/L2_CACHE_L2_Norm_Likwid_Performance.dat
    likwid-perfctr -m -f -g L2CACHE -C 0 ./pdeSolver -nx $nx_ny -ny $nx_ny -i 10 -o arquivo_saida | grep "L2 miss ratio" | head -1 | grep -o -P "[0-9]+\.[0-9]+" >> ./likwidPerformance/L2_CACHE_Gauss_Seidel_Likwid_Performance.dat
    likwid-perfctr -m -f -g L2CACHE -C 0 ./pdeSolver -nx $nx_ny -ny $nx_ny -i 10 -o arquivo_saida | grep "L2 miss ratio" | tail -1 | grep -o -P "[0-9]+\.[0-9]+" >> ./likwidPerformance/L2_CACHE_L2_Norm_Likwid_Performance.dat
done

# [FLOPS_DP].
for nx_ny in ${array[*]}
do
    printf "%s " $nx_ny >> ./likwidPerformance/DP_MFLOPs_Gauss_Seidel_Likwid_Performance.dat
    printf "%s " $nx_ny >> ./likwidPerformance/DP_MFLOPs_L2_Norm_Likwid_Performance.dat
    likwid-perfctr -m -f -g FLOPS_DP -C 0 ./pdeSolver -nx $nx_ny -ny $nx_ny -i 10 -o arquivo_saida  | grep -P "^[^\w]+DP MFLOP" | head -1 | grep -o -P "[0-9]+\.[0-9]+" >> ./likwidPerformance/DP_MFLOPs_Gauss_Seidel_Likwid_Performance.dat
    likwid-perfctr -m -f -g FLOPS_DP -C 0 ./pdeSolver -nx $nx_ny -ny $nx_ny -i 10 -o arquivo_saida  | grep -P "^[^\w]+DP MFLOP" | tail -1 | grep -o -P "[0-9]+\.[0-9]+" >> ./likwidPerformance/DP_MFLOPs_L2_Norm_Likwid_Performance.dat

    printf "%s " $nx_ny >> ./likwidPerformance/AVX_DP_MFLOPs_Gauss_Seidel_Likwid_Performance.dat
    printf "%s " $nx_ny >> ./likwidPerformance/AVX_DP_MFLOPs_L2_Norm_Likwid_Performance.dat
    likwid-perfctr -m -f -g FLOPS_DP -C 0 ./pdeSolver -nx $nx_ny -ny $nx_ny -i 10 -o arquivo_saida | grep "AVX DP MFLOP\/s" | head -1 | grep -o -P "[0-9]+" >> ./likwidPerformance/AVX_DP_MFLOPs_Gauss_Seidel_Likwid_Performance.dat
    likwid-perfctr -m -f -g FLOPS_DP -C 0 ./pdeSolver -nx $nx_ny -ny $nx_ny -i 10 -o arquivo_saida | grep "AVX DP MFLOP\/s" | tail -1 | grep -o -P "[0-9]+" >> ./likwidPerformance/AVX_DP_MFLOPs_L2_Norm_Likwid_Performance.dat
done

# Set the processor frequency.
echo "powersave" > /sys/devices/system/cpu/cpufreq/policy3/scaling_governor

gnuplot gnuplotScript
