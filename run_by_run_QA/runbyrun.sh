#!/bin/bash
# Run-by-run QA script developed by Fudan-BNL group for Isobar data analysis
# author (c) Y. Hu, P. Tribedy, S. Choudhury Feb 12, 2020
# last edit on Mar 01, 2020
ROOTFILE=$1
HISTNAME=$2
NSIGMA_JUMP=$3
NSIGMA_BADRUNS_0=$4
NSIGMA_BADRUNS_1=$5



if [ -z "$1" ]
  then
    printf "\e[31m ERROR: (Argument 1) No ROOT file supplied, will exit \n"
     printf "\e[39m "
     printf "\e[34m Try: bash superjumpcheck.sh FILENAME [e.g qahist.root ] HISTNAME [e.g runidvsrefmult ] \n" 
     printf "\e[39m "
    exit
fi

if [ -z "$2" ]
  then
    printf "\e[31m ERROR: (Argument 2) No HISTNAME supplied, will exit \n"
     printf "\e[39m "
     printf "\e[34m Try: bash superjumpcheck.sh FILENAME [e.g qahist.root ] HISTNAME [e.g runidvsrefmult ] \n" 
     printf "\e[39m "
    exit
fi

if [ -z "$3" ]
  then
    NSIGMA_JUMP=5
fi
if [ -z "$4" ]
  then
    NSIGMA_BADRUNS_0=10
fi
if [ -z "$5" ]
  then
    NSIGMA_BADRUNS_1=5
fi


echo 'WILL USE N-SIGMA =' ${NSIGMA_JUMP} 'AND 1% DIFFERENCE TO DETERMINED JUMPS'
if [ -f jumpcheck_${HISTNAME}_log ]; then
    rm jumpcheck_${HISTNAME}_log
fi

echo 'DOING STEP-1: Precheck the jumps without any cut'
bash jumpcheck.sh ${ROOTFILE} ${HISTNAME} 0 ${NSIGMA_JUMP}>jumpcheck_${HISTNAME}_log
echo 'STEP-1 DONE.'
echo 'DOING STEP-2: 1st round BADRUNS check, '${NSIGMA_BADRUNS_0}'-RMS'
bash badrunfinder.sh ${ROOTFILE} ${HISTNAME} Cut_for_${HISTNAME}.dat 0 ${NSIGMA_BADRUNS_0}>>jumpcheck_${HISTNAME}_log
echo 'STEP-2 DONE.'
echo 'DOING STEP-3: Check the jumps with 1st round BADRUNS rejected'
if [ -f badrunlist_${HISTNAME}_round0.list ]; then
    bash jumpcheck.sh ${ROOTFILE} ${HISTNAME} badrunlist_${HISTNAME}_round0.list ${NSIGMA_JUMP} >>jumpcheck_${HISTNAME}_log
else
    bash jumpcheck.sh ${ROOTFILE} ${HISTNAME} 0 ${NSIGMA_JUMP}>>jumpcheck_${HISTNAME}_log
fi
echo 'STEP-3 DONE.'

echo 'DOING STEP-4: Recheck every region '
#Here we recheck every region to make sure we don't miss any jumps
mv Cut_for_${HISTNAME}.dat Run_region_${HISTNAME}.list
NR_jumps=$(($(cat "Run_region_${HISTNAME}.list" | wc -l)-1))

if [ $NR_jumps -lt 1 ]; then
    echo 'There is no jump, will exit'
    if [ -f badrunlist_${HISTNAME}_round0.list ]; then
    mv badrunlist_${HISTNAME}_round0.list badrunlist_${HISTNAME}_10SIGMA.list
    fi
    exit
fi

root -l "readoutroot.C("'"'"${ROOTFILE}"'"'","'"'"${HISTNAME}"'"'")" -q | awk '{if ($1 == ($1+0)) print $0}'> recheck_${HISTNAME}_raw.dat

if [ -f newcuts_in_recheck ]; then
    rm newcuts_in_recheck
fi
for i in `seq 1 $NR_jumps`
do
    tempnumber=$i
    lowlimit=$(awk 'NR=='${i}' {print($1)}' Run_region_${HISTNAME}.list)
    highlimit=$(awk 'NR=='$((${i}+1))' {print($1)}' Run_region_${HISTNAME}.list)
    awk '{if (($1>='${lowlimit}')&&($1<'${highlimit}')) print $0}' recheck_${HISTNAME}_raw.dat > recheck_region_${tempnumber}.dat

    bash jumpcheck.sh ${ROOTFILE} ${HISTNAME} badrunlist_${HISTNAME}_round0.list ${NSIGMA_JUMP} recheck_region_${tempnumber}.dat >>jumpcheck_${HISTNAME}_log
    
    if [ -f "Cut_for_${HISTNAME}.dat" ]; then
	NR_newcut=$(($(cat "Cut_for_${HISTNAME}.dat" | wc -l)-1))
	if [ $NR_newcut -gt 1 ]; then
	    for ii in `seq 1 $NR_newcut`
	    do
	     awk '{if(NR=='$((${ii}+1))') print ($1)}' Cut_for_${HISTNAME}.dat >> newcuts_in_recheck
	    done
	fi
	rm Cut_for_${HISTNAME}.dat
    fi
    rm recheck_region_${tempnumber}.dat
done

rm recheck_${HISTNAME}_raw.dat

    if [ -f "newcuts_in_recheck" ]; then
	NR_newrecheck=$(cat newcuts_in_recheck | wc -l )
	#echo 'Found '$NR_newrecheck' more jumps in this step'
	cat Run_region_${HISTNAME}.list newcuts_in_recheck | awk '{if($1>1)print($0)}' | sort -n -k1,1 > Final_region_${HISTNAME}.list
	rm Run_region_${HISTNAME}.list newcuts_in_recheck
    else
	mv Run_region_${HISTNAME}.list Final_region_${HISTNAME}.list
    fi
echo 'STEP-4 DONE.'
echo 'DOING STEP-5: Strict BADRUNS check, '${NSIGMA_BADRUNS_1}'-RMS '
###Get the badruns based on the final region
bash badrunfinder.sh ${ROOTFILE} ${HISTNAME} Final_region_${HISTNAME}.list 1 ${NSIGMA_BADRUNS_1} >>jumpcheck_${HISTNAME}_log
echo 'STEP-5 DONE.'

echo '######################'
echo 'Here is the final jumps, saved in Final_region_'${HISTNAME}'.list'
cat Final_region_${HISTNAME}.list
echo '######################'
echo 'Here is the final BADRUNS, saved in badrunlist_'${HISTNAME}'_10SIGMA.list and badrunlist_'${HISTNAME}'_5SIGMA.list '
if [ -f badrunlist_${HISTNAME}_round0.list ]; then
    mv badrunlist_${HISTNAME}_round0.list badrunlist_${HISTNAME}_10SIGMA.list
    cat badrunlist_${HISTNAME}_10SIGMA.list
fi
if [ -f badrunlist_${HISTNAME}_round1.list ]; then
    mv badrunlist_${HISTNAME}_round1.list badrunlist_${HISTNAME}_5SIGMA.list
    cat badrunlist_${HISTNAME}_5SIGMA.list
fi
echo '######################'
echo 'Runing details please find in jumpcheck_'${HISTNAME}'_log' 
