#!/bin/bash

. /packages/pulsarsoft/current_version/login/psr.bashrc;

readarray -t a <25PSRs.txt

for psr in "${a[@]}"
do
    echo $psr;

    # make psr_fitP.tim by adding ../../../psr/
    cp ${psr}/${psr}_all.tim ${psr}/${psr}_fitP.tim;
    sed -i 's/tims/..\/..\/..\/'${psr}'\/tims/' ${psr}/${psr}_fitP.tim;

    echo ${psr} > ${psr}.txt
done

