#!/bin/bash

versionArray=(v1 v2)

# Deleting files.
for version in ${versionArray[@]};
do
    rm -rf ./likwidPerformance/$version/*
done

# Creating files.
for version in ${versionArray[@]};
do
    cd $version
    make clean
    reset
    make
    ./script.sh
    cd likwidPerformance
    cp AVX_DP_MFLOPs.eps DP_MFLOPs.eps L2_CACHE.eps L3.eps ./../../likwidPerformance/$version/
    cd ../..
done
