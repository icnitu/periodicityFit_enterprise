#!/bin/bash

. /packages/pulsarsoft/next_version/login/psr.bashrc;


readarray -t psrs <../25PSRs.txt
readarray -t dmNCs < dmNC_new.txt
readarray -t rnNCs < rnNC_new.txt
readarray -t Aredmins < Aredmin.txt
readarray -t Aredmaxs < Aredmax.txt
readarray -t Admmins < Admmin.txt
readarray -t Admmaxs < Admmax.txt
readarray -t pmins < pmin.txt
readarray -t pmaxs < pmax.txt

for i in 8 13 22 2 9 19 # psr=1600, 1744, 1918, 0751, 1640, 1909
do
    psr=${psrs[i]}

    dmNC=${dmNCs[i]}
    rnNC=${rnNCs[i]}

    Aredmin=${Aredmins[i]}
    Aredmax=${Aredmaxs[i]}
    Admmin=${Admmins[i]}
    Admmax=${Admmaxs[i]}
    pmax=${pmaxs[i]}
    pmin=${pmins[i]}

    echo ${psr}": RN "${rnNC}" & DM "${dmNC}
    echo "Ared: "${Aredmin}" : "${Aredmax}
    echo "Adm: "${Admmin}" : "${Admmax}
    echo "p: "${pmin}" : "${pmax}

#    nice -n 13 ./run_enterprise.py ../${psr}/${psr}.par ../${psr}/${psr}_all.tim -f group --all-corner --plot-chain --dynesty --dynesty-plots -t 32 -P --planets 1 --period-min ${pmin} --period-max ${pmax} --mass-max 1e-1  --mass-log-prior --mass-min 1e-5 --dm --red-ncoeff ${rnNC} --dm-ncoeff ${dmNC} --tspan-mult 1 --ecc-log-prior --white-prior-log --red-prior-log --Ared-max ${Aredmax} --Ared-min ${Aredmin} --dm-prior-log --Adm-max ${Admmax} --Adm-min ${Admmin};


    if [ "$rnNC" -eq 0 ]; then

        nice -n 13 ./run_enterprise.py ../${psr}/${psr}.par ../${psr}/${psr}_all.tim -f group --all-corner --plot-chain --dynesty --dynesty-plots -t 32 -P --planets 1 --period-min ${pmin} --period-max ${pmax} --mass-max 1e-1  --mass-log-prior --mass-min 1e-5 --dm --dm-ncoeff ${dmNC} --tspan-mult 1 --ecc-log-prior --white-prior-log --dm-prior-log --Adm-max ${Admmax} --Adm-min ${Admmin} --no-red-noise;


    else

	if [ "$dmNC" -eq 0 ]; then

            nice -n 13 ./run_enterprise.py ../${psr}/${psr}.par ../${psr}/${psr}_all.tim -f group --all-corner --plot-chain --dynesty --dynesty-plots -t 32 -P --planets 1 --period-min ${pmin} --period-max ${pmax} --mass-max 1e-1  --mass-log-prior --mass-min 1e-5 --red-ncoeff ${rnNC} --tspan-mult 1 --ecc-log-prior --white-prior-log --red-prior-log --Ared-max ${Aredmax} --Ared-min ${Aredmin};

	else

	    nice -n 13 ./run_enterprise.py ../${psr}/${psr}.par ../${psr}/${psr}_all.tim -f group --all-corner --plot-chain --dynesty --dynesty-plots -t 32 -P --planets 1 --period-min ${pmin} --period-max ${pmax} --mass-max 1e-1  --mass-log-prior --mass-min 1e-5 --dm --red-ncoeff ${rnNC} --dm-ncoeff ${dmNC} --tspan-mult 1 --ecc-log-prior --white-prior-log --red-prior-log --Ared-max ${Aredmax} --Ared-min ${Aredmin} --dm-prior-log --Adm-max ${Admmax} --Adm-min ${Admmin};

        fi

    fi

done


