#!/bin/bash


psr_list=('J0751+1807' 'J0900-3144' 'J1012+5307' 'J1022+1001' 'J1600-3053' 'J1640+2224' 'J1713+0747' 'J1744-1134' 'J1909-3744' 'J1918-0642');
#psr_list=('J0900-3144')

for index in "${!psr_list[@]}"; do

    psr="${psr_list[index]}";

    echo ${psr};

    echo "Getting chain.h5...";
    rsync -sarv wolfrock:/nvme1/initu/EPTA_25PSRs_runs/QP_runs/${psr}/chain.h5 ${psr}/chain.h5;

done
