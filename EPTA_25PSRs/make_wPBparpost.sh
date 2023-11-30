#!/bin/bash

#psrs=('J0751+1807' 'J0900-3144' 'J1012+5307' 'J1022+1001' 'J1600-3053' 'J1640+2224' 'J1713+0747' 'J1744-1134' 'J1909-3744' 'J1918-0642')
psrs=('J1713+0747')

for i in ${!psrs[@]}
do
    psr=${psrs[$i]};
    echo ${psr};
    binaryfile=${psr}_BINARY.txt;
    postfile0=${psr}.par.post;
    postfilewPB=${psr}_wPB.par.post;
    wPBwfit=${psr}_wPB_post_wfit.par;
    noPB0=${psr}_noPB_post_wfit.par;
    noPB0_nofit=${psr}_noPB_nomorefit.par;


    cd /home/mbcx4in2/EPTA_25PSRs_runs/planet_runs/dynesty-run_enterprise/withP/${psr};

    pwd;

    rm ${postfilewPB} ${wPBwfit} ${noPB0} ${noPB0_nofit};
    sed '/BINARY/d' ${postfile0} > ${postfilewPB};

    cat ${psr}_BINARY.txt  >> ${postfilewPB};

    tempo2 -f ${postfilewPB} ${psr}_all.tim -newpar;
    mv new.par ${wPBwfit};

    cp ${wPBwfit} ${noPB0};

    if grep -q '^PB_2' ${noPB0}; then
        echo 'Already had binary\n';
        sed -i '/^PB_2/d' ${noPB0};
        sed -i '/^T0_2/d' ${noPB0};
        sed -i '/^A1_2/d' ${noPB0};
        sed -i '/^OM_2/d' ${noPB0};
        sed -i '/^ECC_2/d' ${noPB0};

    else
        echo 'Did NOT already have binary\n';
        sed -i '/^PB/d' ${noPB0};
        sed -i '/^T0/d' ${noPB0};
        sed -i '/^A1/d' ${noPB0};
        sed -i '/^OM/d' ${noPB0};
        sed -i '/^ECC/d' ${noPB0};
        sed -i '/^BINARY/d' ${noPB0};
        
    fi

    cp ${noPB0} dummy.txt;

    sed -i 's/ 1 / /g' dummy.txt;
    ../../disableparam.py dummy.txt 'JUMP' >> ${noPB0_nofit};
    rm dummy.txt;


    tempo2 -output general2 -f ${noPB0_nofit} ${psr}_all.tim -s '{sat} {post} {err} {freq} zzz\n' | grep 'zzz' >> residuals_${psr}_noPB_nofit.txt;

done
