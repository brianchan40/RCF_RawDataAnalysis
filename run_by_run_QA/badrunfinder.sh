#!/bin/bash
# Run-by-run QA script developed by Fudan-BNL group for Isobar data analysis
# author (c) Y. Hu, P. Tribedy, S. Choudhury Feb 12, 2020
# last edit on Mar 01, 2020

file=$1
input=$2
region=$3
round=$4
nsigma=$5
badrunlist=$6



if [ -z "$1" ]
  then
    printf "\e[31m ERROR: (Argument 1) No ROOT file supplied, will exit \n"
     printf "\e[39m "
     printf "\e[34m Try: bash badrunfinder.sh ROOTFILE [e.g qahit.root] HISTNAME [e.g runidvsrefmult]  RUN_REGIONS [e.g run_regions.list] ROUND [e.g 1] NSIGMA [e.g 5, 10] BADRUNLIST\n" 
     printf "\e[39m "
    exit
fi

if [ -z "$2" ]
  then
    printf "\e[31m ERROR: (Argument 2) No HISTNAME supplied, will exit \n"
     printf "\e[39m "
     printf "\e[34m Try: bash badrunfinder.sh ROOTFILE [e.g qahit.root] HISTNAME [e.g runidvsrefmult]  RUN_REGIONS [e.g run_regions.list] ROUND [e.g 1] NSIGMA [e.g 5, 10] BADRUNLIST \n" 
     printf "\e[39m "
    exit
fi

if [ -z "$3" ]
  then
    printf "\e[31m ERROR: (Argument 3) No RUN-REGIONS supplied, will exit \n"
     printf "\e[39m "
     printf "\e[34m Try: bash badrunfinder.sh ROOTFILE [e.g qahit.root] HISTNAME [e.g runidvsrefmult]  RUN_REGIONS [e.g run_regions.list] ROUND [e.g 1]  NSIGMA [e.g 5, 10] BADRUNLIST\n" 
     printf "\e[39m "
    exit
fi


if [ -z "$4" ]
then
    printf "\e[33m WARNING: (Argument 4) No ROUND supplied, round=0 reject runs with 10-RMS limit region-by-region, round=1 reject runs with 5-Weighted-Error limit region-by-region. We will set round=0 as the default \n"
    printf "\e[39m "
    round=0 
fi

if [ -z "$5" ]
then
    printf "\e[33m WARNING: (Argument 5) No NSIGMA supplied, round=0 reject runs with 10-RMS limit region-by-region, round=1 reject runs with 5-RMS limit region-by-region. \n"
    printf "\e[39m "
    if [ ${round} -eq 1 ]
    then
	n_sigma=5
	printf "\e[39m Will set NSIGMA as 5 \n"
	printf "\e[39m "
    else
	n_sigma=10
	printf "\e[39m Will set NSIGMA as 10 \n"
	printf "\e[39m "
    fi
else
    n_sigma=${nsigma}
fi

if [ ${round} -eq 1 ]
then
    if [ -z "$6" ]
    then
	printf "\e[33m WARNING: (Argument 5) This is not the first round, but No BADRUNSLIST supplied. \n"
	printf "\e[39m "
	nobadruns=1
    else
	nobadruns=0
    fi  
else
    nobadruns=1
fi

if [ -f ${input}_raw_0.dat ]; then
    rm ${input}_raw_0.dat
fi

root -l "readoutdata.C("'"'"${file}"'"'","'"'"${input}"'"'","'"'"${region}"'"'")" -q | awk '{if ($1 == ($1+0)) print $0}'>"${input}_raw_0.dat"
cat ${region} | awk '{if ($1 == ($1+0)) print $0}' | sort -n -k1,1 > temp_run_regions
NR_runregions=$(($(cat "temp_run_regions" | wc -l )-1))

if [ ${nobadruns} -eq 1 ]
then
    mv ${input}_raw_0.dat ${input}_newraw.dat
else
    if [ -f badruns_runid ]; then
	rm badruns_runid
    fi
    cat "${badrunlist}" | awk '{if ($1 == ($1+0)) print $1}' > badruns_runid
    grep -vf badruns_runid ${input}_raw_0.dat  >"${input}_newraw.dat"
    rm badruns_runid ${input}_raw_0.dat
fi


if [ -f local_mean_sigma.txt ]; then
    rm local_mean_sigma.txt
fi

if [ -f badrunlist_${input}_round${round}.list ]; then
    rm badrunlist_${input}_round${round}.list
fi

 
for i in `seq 1 $NR_runregions`
do
    tempnumber=$i
    lowlimit=$(awk 'NR=='${i}' {print($1)}' ${region})
    highlimit=$(awk 'NR=='$((${i}+1))' {print($1)}' ${region})
    awk '{if (($1>='${lowlimit}')&&($1<'${highlimit}')) print $0}' ${input}_newraw.dat > region_${tempnumber}.dat
    runnb=$(cat region_${tempnumber}.dat | wc -l )
#    echo $runnb 
#    if we use weighted sigma as the badruns limit
#    awk '{if ($3!=0){ u += $2*(1/($3*$3)); w += (1/($3*$3)); m = u/w; s += (1/($3*$3)); sig = (1/s)^0.5} } END {print '$tempnumber',m, sig*'$runnb'^0.5}' "region_${tempnumber}.dat" >> local_mean_sigma.txt

#    if we use fitting RMS as the badruns limit 
    awk '{if ($3!=0){ u += $2*(1/($3*$3)); w += (1/($3*$3)); m = u/w; mean=$4; rms = $5;  s += (1/($3*$3)); sig = (1/s)^0.5 } } END {print '$tempnumber',m, rms, sig*'$runnb'^0.5, mean}' "region_${tempnumber}.dat" >> local_mean_sigma.txt 

    #echo ${n_sigma}
    mean=$(awk 'NR=='${i}' {print($2)}' local_mean_sigma.txt)
    sigma=$(awk 'NR=='${i}' {print($3)}' local_mean_sigma.txt)
    if [ ${round} -eq 0 ]
    then
	awk '{if (($2>'${mean}'+'${n_sigma}'*$5)||($2<'${mean}'-'${n_sigma}'*$5)) print( $1, $2, $3,'${mean}'-'${n_sigma}'*$5, '${mean}'+'${n_sigma}'*$5)  }' "region_${tempnumber}.dat" > badruns_${tempnumber}.dat
    else
	awk '{if (($2>'${mean}'+'${n_sigma}'*'${sigma}')||($2<'${mean}'-'${n_sigma}'*'${sigma}')) print( $1, $2, $3,'${mean}'-'${n_sigma}'*'${sigma}', '${mean}'+'${n_sigma}'*'${sigma}')  }' "region_${tempnumber}.dat"> badruns_${tempnumber}.dat
    fi
    if [ ${round} -eq 0 ]; then
    awk '{print($1,"#'${input}'_badruns", "'${n_sigma}'SIGMA" )}' "badruns_${tempnumber}.dat" >> badrunlist_${input}_round${round}.list
    else
    awk '{print($1,"#'${input}'_badruns", "'${n_sigma}'SIGMA" )}' "badruns_${tempnumber}.dat" >> badrunlist_${input}_round${round}.list
    fi
    #cache file clean up
    rm badruns_${tempnumber}.dat
    rm region_${tempnumber}.dat
done

if [ -f badrunlist_${input}_round${round}.list ]; then
    Nb_badruns=$(cat "badrunlist_${input}_round${round}.list" | wc -l)
    firstrun=$(cat ${input}_newraw.dat | sort -n -k1,1 | awk 'NR==1 {print ($1)}')
    lastrun=$(cat ${input}_newraw.dat | sort -n -r -k1,1 | awk 'NR==1 {print ($1)}')
    printf "\e[34m We found ${Nb_badruns} badruns for ${input} in round${round} from Runid = ${firstrun} to Runid = ${lastrun} with NSIGMA = ${n_sigma}\e[39m \n"
    cat badrunlist_${input}_round${round}.list
    else
    printf "\e[34m We found 0 badruns for ${input} in round${round} from Runid = ${firstrun} to Runid = ${lastrun} with NSIGMA = ${n_sigma}\e[39m \n"
fi


#cache file clean up
rm temp_run_regions
rm ${input}_newraw.dat
#exit
if [ -f local_mean_sigma.txt ]; then
    rm local_mean_sigma.txt
fi
