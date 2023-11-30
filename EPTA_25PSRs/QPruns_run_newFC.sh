#!/bin/bash

. /packages/pulsarsoft/current_version/login/psr.bashrc;


readarray -t psrs <../25PSRs.txt
readarray -t dmNCs < ../planet_runs/run_enterprise/dmNC_new.txt
readarray -t rnNCs < ../planet_runs/run_enterprise/rnNC_new.txt
readarray -t Aredmins < ../planet_runs/run_enterprise/Aredmin.txt
readarray -t Aredmaxs < ../planet_runs/run_enterprise/Aredmax.txt
readarray -t Admmins < ../planet_runs/run_enterprise/Admmin.txt
readarray -t Admmaxs < ../planet_runs/run_enterprise/Admmax.txt

#for psr in "${a[i]}"

for i in 8 13 22 2 9 19 # psr=1600, 1744, 1918, 0751, 1640, 1909
do
    psr=${psrs[i]}

    dmNC=${dmNCs[i]}
    rnNC=${rnNCs[i]}

    Aredmin=${Aredmins[i]}
    Aredmax=${Aredmaxs[i]}
    Admmin=${Admmins[i]}
    Admmax=${Admmaxs[i]}

    echo ${psr}": RN "${rnNC}" & DM "${dmNC}
    echo "Ared: "${Aredmin}"-"${Aredmax}
    echo "Adm: "${Admmin}"-"${Admmax}

    nice -n 12 ../planet_runs/run_enterprise/run_enterprise.py ${psr}/${psr}.par ${psr}/${psr}_all.tim -f group --white-prior-log --red-prior-log --Ared-max  ${Aredmax} --Ared-min ${Aredmin} --red-ncoeff 100 --dm-prior-log --Adm-max ${Admmax} --Adm-min ${Admmin} --dm-ncoeff ${dmNC} --qp --emcee -N 4000 --nwalkers 500 -t 32 --all-corner --plot-chain --plot-derived --outdir ${psr}/ --dm --tspan-mult 1 --cont;

:'
    if [ "$rnNC" -eq 0 ]; then


	nice -n 12 ../planet_runs/run_enterprise/run_enterprise.py ${psr}/${psr}.par ${psr}/${psr}_all.tim -f group --white-prior-log --dm-prior-log --Adm-max ${Admmax} --Adm-min ${Admmin} --dm-ncoeff ${dmNC} --qp --emcee -N 2000 --nwalkers 500 -t 32 --all-corner --plot-chain --plot-derived --outdir ${psr}/ --dm --tspan-mult 1 --no-red-noise;
       
    else

        if [ "$dmNC" -eq 0 ]; then

	nice -n 12 ../planet_runs/run_enterprise/run_enterprise.py ${psr}/${psr}.par ${psr}/${psr}_all.tim -f group --white-prior-log --red-prior-log --Ared-max  ${Aredmax} --Ared-min ${Aredmin} --red-ncoeff ${rnNC} --qp --emcee -N 2000 --nwalkers 500 -t 32 --all-corner --plot-chain --plot-derived --outdir ${psr}/ --tspan-mult 1;



	else

	nice -n 12 ../planet_runs/run_enterprise/run_enterprise.py ${psr}/${psr}.par ${psr}/${psr}_all.tim -f group --white-prior-log --red-prior-log --Ared-max  ${Aredmax} --Ared-min ${Aredmin} --red-ncoeff ${rnNC} --dm-prior-log --Adm-max ${Admmax} --Adm-min ${Admmin} --dm-ncoeff ${dmNC} --qp --emcee -N 2000 --nwalkers 500 -t 32 --all-corner --plot-chain --plot-derived --outdir ${psr}/ --dm --tspan-mult 1;


 	fi

    fi    
'
done

