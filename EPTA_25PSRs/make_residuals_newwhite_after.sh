#!/bin/bash
. /packages/pulsarsoft/current_version/login/psr.bashrc

#psrs=('J0751+1807' 'J0900-3144' 'J1012+5307' 'J1022+1001' 'J1600-3053' 'J1640+2224' 'J1713+0747' 'J1744-1134' 'J1909-3744' 'J1918-0642')
#psrs=('J1744-1134')
psrs=('J1600-3053')

for i in ${!psrs[@]}
do
    psr=${psrs[$i]};
    echo ${psr};
    binaryfile=${psr}_p2_BINARY.txt;
    postfile0=${psr}.par.post;
    notnew=notnew.par;
    new=new.par;   
    newnoP=new_noP.par;
    newnoPnofit=new_noP_nofit.par;   
    noglf2=nogl_nof2.par;
    white=white.par;
    newwhite=newwhite.par;
    newwhiteP=newwhite_planet.par;


    cd /nvme1/initu/EPTA_25PSRs_runs/planet_runs/dynesty-run_enterprise/withP/${psr}_p2;

    pwd;

    rm residuals{0,_newwhite,_planet}.txt;
#    sed '/BINARY/d' ${postfile0} > ${notnew};

#    cat ../${psr}_BINARY.txt  >> ${notnew};

#    tempo2 -f ${notnew} ${psr}_all.tim -newpar;
#    mv new.par ${new};

#    cp ${new} ${newnoP};

#    if grep -q '^PB_2' ${newnoP}; then
#        echo 'Already had binary\n';
#        sed -i '/^PB_2/d' ${newnoP};
#        sed -i '/^T0_2/d' ${newnoP};
#        sed -i '/^A1_2/d' ${newnoP};
#        sed -i '/^OM_2/d' ${newnoP};
#        sed -i '/^ECC_2/d' ${newnoP};

#    else
#        echo 'Did NOT already have binary\n';
#        sed -i '/^PB/d' ${newnoP};
#        sed -i '/^T0/d' ${newnoP};
#        sed -i '/^A1/d' ${newnoP};
#        sed -i '/^OM/d' ${newnoP};
#        sed -i '/^ECC/d' ${newnoP};
#        sed -i '/^BINARY/d' ${newnoP};
        
#    fi

#    cp ${newnoP} dummy.txt;

#    sed -i 's/ 1 / /g' dummy.txt;
#    ../../disableparam.py dummy.txt 'JUMP' >> ${newnoPnofit};
#    rm dummy.txt;


#    tempo2 -output general2 -f ${newnoP} ${psr}_all.tim -s '{sat} {pre} {err} {freq} zzz\n' | grep 'zzz' >> residuals0.txt;

#    grep -vE '^F2|^GL' ${new} > ${noglf2}; # shouldn't really do anything here..

#    ../../make_pulsar_plots.py ${new} ${psr}_all.tim ${noglf2};

#    cp ${white} ${newwhite};


# ONLY AFTER EDITING NEWWHITE!

    tempo2 -output general2 -f ${newnoP} ${psr}_all.tim -s '{sat} {pre} {err} {freq} zzz\n' | grep 'zzz' >> residuals0.txt;
    
    if grep -q '^PB_2' ${newwhite}; then
        echo 'Already had binary\n';
        sed -i '/^PB_2/d' ${newwhite};
        sed -i '/^T0_2/d' ${newwhite};
        sed -i '/^A1_2/d' ${newwhite};
        sed -i '/^OM_2/d' ${newwhite};
        sed -i '/^ECC_2/d' ${newwhite};

    else
        echo 'Did NOT already have binary\n';
        sed -i '/^PB/d' ${newwhite};
        sed -i '/^T0/d' ${newwhite};
        sed -i '/^A1/d' ${newwhite};
        sed -i '/^OM/d' ${newwhite};
        sed -i '/^ECC/d' ${newwhite};
        sed -i '/^BINARY/d' ${newwhite};
        
    fi

    tempo2 -output general2 -f ${newwhite} tndm.tim -s '{sat} {pre} {err} {freq} zzz\n' | grep 'zzz' >> residuals_newwhite.txt;

    sed '/BINARY/d' ${newwhite} > ${newwhiteP};

    cat ../${psr}_p2_BINARY.txt  >> ${newwhiteP};

    tempo2 -output general2 -f ${newwhiteP} tndm.tim -s '{sat} {pre} {err} {freq} zzz\n' | grep 'zzz' >> residuals_planet.txt;

        

done
