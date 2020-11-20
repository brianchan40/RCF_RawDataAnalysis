#!/bin/bash
# Run-by-run QA script developed by Fudan-BNL group for Isobar data analysis
# authors (c) Fudan-BNL blind analysis team (Y. Hu, P. Tribedy, S. Choudhury)
### PLEASE WRITE TO <yuhu@bnl.gov> , <ptribedy@bnl.gov> , <ptribedy@icloud.com>
# first created Yu Hu, Feb 12, 2020
# edited by Yu Hu, Mar 01, 2020
# edited by Yu Hu, Mar 30, 2020
# last edit by P. Tribedy, Mar 30, 2020, ready to be frozen, this is the version 3.0 
# 



ROOTFILE=$1
NAMELIST=$2
NSIGMA_JUMP=$3
NSIGMA_BADRUNS_0=$4
NSIGMA_BADRUNS_1=$5


    echo '-------------------------------------------------------------------------------'
    echo '-------------------------------------------------------------------------------'
    printf "\e[35m-------------------Run-by-Run QA script for STAR data analysis ----------------\n"
    printf "\e[39m-----------------------------------version 3.1---------------------------------\n"
    echo '-------------------Input1: ROOT file with profile histograms vs runid ---------'
    echo '-------------------Input2: a list of profile histograms -----------------------'
    echo '-------------------Output: lists of bad runs & stable regions -----------------'
    echo '-------------------------------------------------------------------------------'
    echo '--------------------------------Example to run --------------------------------'
    printf "\e[35m--------------bash runbyrun_advanced.sh 27gev_qahist.root qa.list--------------\n"
    printf "\e[39m-------------------------------------------------------------------------------\n"
    echo '-----------------------Developed for Isobar blind analysis --------------------'
    echo '-----------------------By: Fudan-BNL analysis team, Feb 2020 ------------------'
    printf "\e[35m-------------------Contact: <yuhu@bnl.gov>, <ptribedy@bnl.gov> -----------------\n"
    printf "\e[39m-------------------------------------------------------------------------------\n"
    echo '-------------------------------------------------------------------------------'
    echo '-------------------------------------------------------------------------------'

if [ -z "$1" ]
  then
    printf "\e[31m ERROR: (Argument 1) No ROOT file supplied, will exit \n"
     printf "\e[39m "
     printf "\e[34m Try: bash runbyrun_advanced.sh FILENAME [e.g qahist.root ] NAMELIST [e.g name.list ] \n" 
     printf "\e[35m e.g.:  bash runbyrun_advanced.sh 27gev_qahist.root qa.list \n"
     printf "\e[39m "
    exit
fi

if [ -z "$2" ]
  then
    printf "\e[31m ERROR: (Argument 2) No NAMELIST supplied, will exit \n"
     printf "\e[39m "
     printf "\e[34m Try: bash runbyrun_advanced.sh FILENAME [e.g qahist.root ] NAMELIST [e.g name.list ] \n" 
     printf "\e[35m e.g.:  bash runbyrun_advanced.sh 27gev_qahist.root qa.list \n"
     printf "\e[39m "
    exit
fi

if [ -z "$4" ]
then
    NSIGMA_BADRUNS_0=10
fi

if [ -z "$5" ]
then
    NSIGMA_BADRUNS_1=5
fi

# need to check if you have the readoutnrms script
if [ ! -f "readoutnrms.C" ]
then
    printf "\e[31m ERROR: the file readoutnrms.C does not exists , will exit \n"
    printf "\e[39m "
    exit
fi


# Check how many quantities we have in the list
NR_HISTNAME=$(cat ${NAMELIST} | wc -l)

#Need to make sure the first step have jumps, use a until loop to do it
NR_jumps=0
NR_qalist=1
until [ $NR_jumps -gt 0 ]
do
    HISTNAME=$(awk 'NR=='${NR_qalist}' {print($2)}' ${NAMELIST})
    echo ${HISTNAME}
    #Use the 27 GeV data as the standard. check how many times of weighted error we should use
    if [ -z "$3" ]
    then
	NSIGMA_JUMP=$(root -l "readoutnrms.C("'"'"${ROOTFILE}"'"'","'"'"${HISTNAME}"'"'","'"'"${ROOTFILE}"'"'")" -q | awk '{if ($1 == ($1+0)) print $1}' )
	echo 'FOR AUTO CHECK, WE WILL USE N-SIGMA =' ${NSIGMA_JUMP} 'AND 1% DIFFERENCE TO DETERMINE JUMPS'
	#      NSIGMA_JUMP=5
    fi

    #Here we can give the data sanity index for this data set
    SANITYINDEX=$(root -l "readoutnrms.C("'"'"${ROOTFILE}"'"'","'"'"${HISTNAME}"'"'","'"'"${ROOTFILE}"'"'")" -q | awk '{if ($2 == ($2+0)) print $2}' )

if [ $NR_qalist == 1 ]
then
    echo '--------------------------------------------------------------------------'
    echo '---------------------How good/bad is this dataset ?-----------------------'
    echo '--------------------------------------------------------------------------'
    echo 'SanityIndex -->1 : Good dataset  (e.g. Run 18 27 GeV data ~288, Run 16 200 GeV ~727)'
    printf "\e[31m The estimated SanityIndex for this dataset using the quantity <${HISTNAME}> is '${SANITYINDEX}'\n"
    printf "\e[39m "
    echo '--------------------------------------------------------------------------'
fi    

    if [ -f jumpcheck_${HISTNAME}_log ]; then
	rm jumpcheck_${HISTNAME}_log
    fi

    # Here we start to check the fist (${NR_qalist}) qualities with the default logic
    echo 'DOING STEP-1: Prechecking the jumps without any cuts ('${NSIGMA_JUMP}' weighted error)'
    bash jumpcheck.sh ${ROOTFILE} ${HISTNAME} 0 ${NSIGMA_JUMP}>jumpcheck_${HISTNAME}_log
    echo 'STEP-1 DONE.'
    echo 'DOING STEP-2: Doing 1st round of BADRUNS checking, '${NSIGMA_BADRUNS_0}'-RMS'
    bash badrunfinder.sh ${ROOTFILE} ${HISTNAME} Cut_for_${HISTNAME}.dat 0 ${NSIGMA_BADRUNS_0}>>jumpcheck_${HISTNAME}_log
    echo 'STEP-2 DONE.'
    echo 'DOING STEP-3: Checking the jumps for 1st round of BADRUNS rejection'
    if [ -f badrunlist_${HISTNAME}_round0.list ]; then
	bash jumpcheck.sh ${ROOTFILE} ${HISTNAME} badrunlist_${HISTNAME}_round0.list ${NSIGMA_JUMP} >>jumpcheck_${HISTNAME}_log
    else
	bash jumpcheck.sh ${ROOTFILE} ${HISTNAME} 0 ${NSIGMA_JUMP}>>jumpcheck_${HISTNAME}_log
    fi
    echo 'STEP-3 DONE.'

    echo 'DOING STEP-4: Rechecking every region '
    #Here we recheck every region to make sure we don't miss any jumps
    mv Cut_for_${HISTNAME}.dat Run_region_${HISTNAME}.list
    NR_jumps=$(($(cat "Run_region_${HISTNAME}.list" | wc -l)-1))

    #    echo 'test for NR_jumps '${NR_jumps}''
    NR_qalist=$((${NR_qalist}+1))    
    if [ $NR_jumps -lt 1 ]; then
	#Here if the first quantity couldn't find the jump, we should visit the next quantity
	#Case 1
	echo 'There is no jump associated with this qantity <'${HISTNAME}'>, we will search in the next one'
	echo '######################'
	echo '  '
	if [ -f badrunlist_${HISTNAME}_round0.list ]; then
	    mv badrunlist_${HISTNAME}_round0.list badrunlist_${HISTNAME}_10SIGMA.list
	fi

	#Here we also need to check if it's the last quantity in the qa.list, we need to exit this code
	if [ ${NR_qalist} -gt ${NR_HISTNAME} ]; then
	    echo 'We went through all the quantities that you have, but find no jumps, will exit'
	    exit
	fi
	
	#Case 1 finished
    else
	#But if we find some jumps in this step, we should do the following search
	#Case 2
	root -l "readoutroot.C("'"'"${ROOTFILE}"'"'","'"'"${HISTNAME}"'"'")" -q | awk '{if ($1 == ($1+0)) print $0}'> recheck_${HISTNAME}_raw.dat
	
	if [ -f newcuts_in_recheck ]; then
	    rm newcuts_in_recheck
	fi
	for i in `seq 1 $NR_jumps`
	do
	    tempnumber=$i
	    lowlimit=$(awk 'NR=='${i}' {print($1)}' Run_region_${HISTNAME}.list)
	    highlimit=$(awk 'NR=='$((${i}+1))' {print($1)}' Run_region_${HISTNAME}.list)
	    awk '{if (($1>='${lowlimit}')&&($1<='${highlimit}')) print $0}' recheck_${HISTNAME}_raw.dat > recheck_region_${tempnumber}.dat

	    NR_eachregion=$(cat recheck_region_${tempnumber}.dat | wc -l)
	    #We need at least 3 runs to do the jumpcheck
	    if [ ${NR_eachregion} -lt 3 ]; then
		rm recheck_region_${tempnumber}.dat
		continue
	    fi
	
	    bash jumpcheck.sh ${ROOTFILE} ${HISTNAME} badrunlist_${HISTNAME}_round0.list ${NSIGMA_JUMP} recheck_region_${tempnumber}.dat >>jumpcheck_${HISTNAME}_log
	
	    if [ -f "Cut_for_${HISTNAME}.dat" ]; then
		NR_newcut=$(($(cat "Cut_for_${HISTNAME}.dat" | wc -l)-1))
		if [ $NR_newcut -gt 1 ]; then
		    for ii in `seq 1 $NR_newcut`
		    do
			awk '{if(NR=='$((${ii}+1))'){if($1>1)print ($1)}}' Cut_for_${HISTNAME}.dat >> newcuts_in_recheck
		    done
		fi
		rm Cut_for_${HISTNAME}.dat
	    fi
	    rm recheck_region_${tempnumber}.dat
	done # recheck for every region for the jumps

	rm recheck_${HISTNAME}_raw.dat

	if [ -f "newcuts_in_recheck" ]; then
	    NR_newrecheck=$(cat newcuts_in_recheck | wc -l )
	    #echo 'Found '$NR_newrecheck' more jumps in this step'
	    cat Run_region_${HISTNAME}.list newcuts_in_recheck | awk '{if($1>1)print($0)}' | sort -u -k1,1 -n > Final_region_${HISTNAME}.list
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
	echo 'Here is the final jumps for '${HISTNAME}', saved in Final_region_'${HISTNAME}'.list'
	cat Final_region_${HISTNAME}.list
	echo '######################'
	echo 'Here is the final BADRUNS for '${HISTNAME}', saved in badrunlist_'${HISTNAME}'_10SIGMA.list and badrunlist_'${HISTNAME}'_5SIGMA.list '
	if [ -f badrunlist_${HISTNAME}_round0.list ]; then
	    mv badrunlist_${HISTNAME}_round0.list badrunlist_${HISTNAME}_10SIGMA.list
	fi
	if [ -f badrunlist_${HISTNAME}_round1.list ]; then
	    mv badrunlist_${HISTNAME}_round1.list badrunlist_${HISTNAME}_5SIGMA.list
	fi
	echo 'Runing details please find in jumpcheck_'${HISTNAME}'_log'

	#Need to track the number of jumps to go out of this loop
	NR_jumps=$(($(cat "Final_region_${HISTNAME}.list" | wc -l)-2))
    fi #Case 2 finished

done
##########  PART II ########

if [ -f badrunlist_${HISTNAME}_10SIGMA.list  -a  -f badrunlist_${HISTNAME}_5SIGMA.list ]; then
    cat badrunlist_${HISTNAME}_5SIGMA.list badrunlist_${HISTNAME}_10SIGMA.list | awk '{if($1>1)print($0)}' | sort -u -k1,1 -n >Final_badruns.list
    cat Final_badruns.list
else
    " ">Final_badruns.list
fi

cat Final_region_${HISTNAME}.list | awk '{print($1, "#'${HISTNAME}'")}' > Final_regions.list


### From here we start to loop all the other quantities

for i in `seq $NR_qalist $NR_HISTNAME`
do
    # Get the next quatity name
    HISTNAME=$(awk 'NR=='${i}' {print($2)}' ${NAMELIST})
    echo ' '
    echo '######################'
    echo 'Working on quantity '${i}', <'${HISTNAME}'>'

    # Check if we have set the n-weighted-error, if not get the corrected n 
    if [ -z "$3" ]
    then
	NSIGMA_JUMP=$(root -l "readoutnrms.C("'"'"${ROOTFILE}"'"'","'"'"${HISTNAME}"'"'","'"'"${ROOTFILE}"'"'")" -q | awk '{if ($1 == ($1+0)) print $1}' )
	echo 'BY AUTO CHECK, WILL USE N-SIGMA =' ${NSIGMA_JUMP} 'AND 1% DIFFERENCE TO DETERMINE JUMPS FOR '${HISTNAME}' '
    else
	echo 'Will use the given N-weighted-error '${NSIGMA_JUMP}''
    fi
    
    
    if [ -f jumpcheck_${HISTNAME}_log ]; then
	rm jumpcheck_${HISTNAME}_log
    fi

    # Using the region list generated in the last step, recheck every region
    NR_jumps=$(($(cat "Final_regions.list" | wc -l)-1))
        
    # First you need to generate the txt file for the new quantity by using readoutroot.C
    root -l "readoutroot.C("'"'"${ROOTFILE}"'"'","'"'"${HISTNAME}"'"'")" -q | awk '{if ($1 == ($1+0)) print $0}'> recheck_${HISTNAME}_raw.dat

    if [ -f newcuts_in_recheck ]; then
	rm newcuts_in_recheck
    fi
    # Check every region we have in the last step
    for ii in `seq 1 $NR_jumps`
    do	
	tempnumber=$ii
	lowlimit=$(awk 'NR=='${ii}' {print($1)}' Final_regions.list)
	highlimit=$(awk 'NR=='$((${ii}+1))' {print($1)}' Final_regions.list)
	awk '{if (($1>='${lowlimit}')&&($1<='${highlimit}')) print $0}' recheck_${HISTNAME}_raw.dat > recheck_region_${tempnumber}.dat

	NR_eachregion=$(cat recheck_region_${tempnumber}.dat | wc -l)
	#We need at least 3 runs to do the jumpcheck, but if one region has less than 10 runs, the statistics are not enough to get the weight, we set the low limit here
	if [ ${NR_eachregion} -lt 10 ]; then
	    rm recheck_region_${tempnumber}.dat
	    continue
	fi

	#	cp recheck_region_${tempnumber}.dat recheck_lasttest_${tempnumber}.dat
	#       debug_NR=$(cat recheck_region_${tempnumber}.dat | wc -l )
	# 	echo "working on region ${ii}, have ${debug_NR} runs"
	
	bash jumpcheck.sh ${ROOTFILE} ${HISTNAME} Final_badruns.list ${NSIGMA_JUMP} recheck_region_${tempnumber}.dat >>jumpcheck_${HISTNAME}_log
	
	if [ -f "Cut_for_${HISTNAME}.dat" ]; then
	    NR_newcut=$(($(cat "Cut_for_${HISTNAME}.dat" | wc -l)-1))
	    if [ $NR_newcut -gt 1 ]; then
		for iii in `seq 1 $NR_newcut`
		do
		    awk '{if(NR=='$((${iii}+1))') print ($1, "#'${HISTNAME}'")}' Cut_for_${HISTNAME}.dat >> newcuts_in_recheck
		done
	    fi
	    rm Cut_for_${HISTNAME}.dat
	fi
	rm recheck_region_${tempnumber}.dat
    done   # Finished the loop for every region
    
    rm recheck_${HISTNAME}_raw.dat # clean the raw data file

    if [ -f "newcuts_in_recheck" ]; then
	NR_newrecheck=$(cat newcuts_in_recheck | wc -l )

	cat Final_regions.list newcuts_in_recheck | awk '{if($1>1)print($0)}' | sort -u -k1,1 -n > Final_region_${HISTNAME}.list
	rm Final_regions.list newcuts_in_recheck
    else
	mv Final_regions.list Final_region_${HISTNAME}.list
    fi
    ###Get the badruns based on the final region
    bash badrunfinder.sh ${ROOTFILE} ${HISTNAME} Final_region_${HISTNAME}.list 1 ${NSIGMA_BADRUNS_1} >>jumpcheck_${HISTNAME}_log

    echo '######################'
    printf "\e[39m Here is the list of jumps after checking ${i} quantities new regions added in Final_regions.list \n"
    printf "\e[39m "
    cat Final_region_${HISTNAME}.list
    echo '######################'
    printf "\e[34m BADRUNS associated with quantity: <$HISTNAME> are saved in badrunlist_${HISTNAME}_5SIGMA.list \n"
    printf "\e[39m "
    if [ -f badrunlist_${HISTNAME}_round1.list ]; then
	mv badrunlist_${HISTNAME}_round1.list badrunlist_${HISTNAME}_5SIGMA.list
    fi
    
    # Initialize for the next quantity
    cp Final_region_${HISTNAME}.list Final_regions.list
    # Merge all the badruns
    if [ -f "badrunlist_${HISTNAME}_5SIGMA.list" ]; then
	NR_new_badruns=$(cat badrunlist_${HISTNAME}_5SIGMA.list | wc -l )
	#echo 'Found '$NR_newrecheck' more jumps in this step'
	cat Final_badruns.list badrunlist_${HISTNAME}_5SIGMA.list | awk '{if($1>1)print($0)}' | sort -u -k1,1 -n > Final_temp_badruns.list
	rm Final_badruns.list
	mv Final_temp_badruns.list Final_badruns.list
    fi
   
done

    printf "\e[39m "
    printf "\e[39m "
    printf "\e[34m Here is the final list of jumps saved in Final_regions.list \n"
    printf "\e[39m "
    cat Final_region_${HISTNAME}.list  
    printf "\e[34m Here is the final list of BADRUNS saved in Final_badruns.list \n"
    printf "\e[39m "
    cat  Final_badruns.list
 
