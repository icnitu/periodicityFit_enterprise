. /packages/pulsarsoft/current_version/login/psr.bashrc;

readarray -t psrs <25PSRs.txt


for psr in "${psrs[@]}"
do

        echo "***** "${psr}" *****"

        rm residuals_${psr}_tndm.txt;

        tempo2 -f ${psr}/${psr}.par.post ${psr}/${psr}_all.tim -newpar;
        
        mv new.par ${psr}/${psr}_post_wfit.par;

        tempo2 -output storeTNDM -f ${psr}/${psr}_post_wfit.par ${psr}/${psr}_all.tim;
        mv tndm.tim ${psr}/${psr}_tndm.tim;
        mv tndm_offset.txt ${psr}/${psr}_tndm_offset.txt;

        tempo2 -output general2 -f ${psr}/${psr}_post_wfit.par ${psr}/${psr}_tndm.tim -s '{sat} {post} {err} {freq} zzz\n' | grep 'zzz' >> residuals_${psr}_tndm.txt;

done

