#!/bin/bash
# Run-by-run QA script developed by Fudan-BNL group for Isobar data analysis
# author (c) Y. Hu, P. Tribedy, S. Choudhury Feb 12, 2020
# last edit on Mar 01, 2020

FILE=$1
input=$2
badrunlist=$3
nsigma=$4
datafile=$5


if [ -z "$1" ]
  then
    printf "\e[31m ERROR: (Argument 1) No ROOT file supplied, will exit \n"
     printf "\e[39m "
     printf "\e[34m Try: bash jumpcheck.sh FILENAME [e.g qahist.root ] HISTNAME [e.g runidvsrefmult ] BADRUNLIST [allbadruns.list] NSIGMA [e.g 5 ] \n" 
     printf "\e[39m "
    exit
fi

if [ -z "$2" ]
  then
    printf "\e[31m ERROR: (Argument 2) No HISTNAME supplied, will exit \n"
     printf "\e[39m "
     printf "\e[34m Try: bash jumpcheck.sh FILENAME [e.g qahist.root ] HISTNAME [e.g runidvsrefmult ] BADRUNLIST [allbadruns.list] NSIGMA [e.g 5 ] \n" 
     printf "\e[39m "
    exit
fi

if [ -z "$3" ]
then
    printf "\e[33m WARNING: (Argument 3) No BADRUNLIST supplied, will not reject the badruns in jumpcheck \n"
    printf "\e[39m "
    nobadruns=1
else
    if [ ${badrunlist} == 0 ]
    then
	printf "\e[33m WARNING: (Argument 3) No BADRUNLIST supplied, will not reject the badruns in jumpcheck \n"
	printf "\e[39m "
	nobadruns=1
    else
	nobadruns=0
    fi
fi

if [ -z "$4" ]
then
    printf "\e[33m WARNING: (Argument 4) No NSIGMA supplied, will try to use the default NSIGMA=5 \n"
    printf "\e[39m "
    nsigma=5 

fi

if [ -z "$5" ]
then
    datamode=0
else
    printf "\e[33m Will use the given data "${datafile}" to run the jumpcheck\n"
    printf "\e[39m "
    datamode=1
fi


if [ ! -f "readoutroot.C" ]
then
    printf "\e[31m ERROR: NEED readoutroot.C, will exit \n"
    printf "\e[39m "
    exit
fi
if [ -f ${input}_raw.dat ]; then
    rm ${input}_raw.dat
fi

#Here we generate the data file from the given root file:
if [ ${datamode} -eq 0 ]
then
    root -l "readoutroot.C("'"'"${FILE}"'"'","'"'"${input}"'"'")" -q > ${input}_raw.dat
    NR_rawdata=$(($(cat "${input}_raw.dat" | wc -l )-4))
else
    cp ${datafile} ${input}_raw.dat
    NR_rawdata=$(cat "${input}_raw.dat" | awk '{if ($1 == ($1+0)) print $1}'| wc -l )
fi


#Here we reject the badruns if we have the badrunlist input
if [ ${nobadruns} -eq 1 ]
then
    cat ${input}_raw.dat | awk '{if ($1 == ($1+0)) print $0}' | awk '{print $0, NR}' | sort -n -k 8,8 >"file_Goodruncut_sort"
else
    cat "${badrunlist}" | awk '{if ($1 == ($1+0)) print $1}' > badruns_runid
    grep -vf badruns_runid ${input}_raw.dat | awk '{if ($1 == ($1+0)) print $0}' | awk '{print $0, NR}' | sort -n -k 8,8 >"file_Goodruncut_sort_old"
    rm badruns_runid
    #After the badruns rejection, the mean need to be updatad
    New_gbmean=$(awk '{if ($3!=0){ u += $2*(1/($3*$3)); w += (1/($3*$3)); m = u/w; s += (1/($3*$3)); sig = (1/s)^0.5 }} END {print m}' "file_Goodruncut_sort_old")
    awk '{print ($1, $2, $3, '$New_gbmean', $5, $6, $7, $8)}' file_Goodruncut_sort_old | sort -n -k 8,8 > "file_Goodruncut_sort"
    rm file_Goodruncut_sort_old
fi

jumpnline=$(cat "file_Goodruncut_sort" | wc -l )
firstrun=$(cat file_Goodruncut_sort | sort -n -k1,1 | awk 'NR==1 {print ($1)}')
lastrun=$(cat file_Goodruncut_sort | sort -n -r -k1,1 | awk 'NR==1 {print ($1)}')

#Redraw the data with normalizing and smoothing.
gnuplot <<EOF	
	set table "Jump_smooth_${input}.dat"
	set samples $jumpnline 
	plot "file_Goodruncut_sort" u (\$8):((\$2)/(\$4)) smooth bezier	
EOF


# #To calculate the 1st derivative, 2nd derivative of the smoothed data
# #To calculate and find out where the smoothed data has the biggest slope or peak 
# #combine all the generated data file

        NR_jump_total=$(($(cat "Jump_smooth_${input}.dat" | wc -l )-6))
	cat Jump_smooth_${input}.dat | awk '{if(NR>4) {print($1, $2)}}' |awk '{ {y1=y2; y2=$2; dy1=dy; d2y=dy2; dy=(y2-y1); dy2=(dy-dy1); slope=dy2*d2y; peak=dy*dy1} if(NR>2&&NR<'${NR_jump_total}') {print($1, $2, dy, dy2, slope, peak, 1.)} }' > "Jump_total_${input}.dat"
	
	#find the biggest slope($5) and peak($6) points
	 awk '{if(($5<0.)||($6<0.)||($3==0.)||($4==0.)) { print int($1+0.5),$2,$3,$4,$5,$6}}' "Jump_total_${input}.dat" | sort -n -k1>"file_temp_peak_slope_${input}.dat"
	join -1 8 -2 1 <(sort -k 8,8 "file_Goodruncut_sort") <(sort -k 1,1 "file_temp_peak_slope_${input}.dat") | sort -n -k1 > "Peak_Slope_Points_jump_${input}.dat"
	#if no peak and slope find will exit
	NR_temp_peakslope=$(cat "file_temp_peak_slope_${input}.dat" | wc -l )
	if [ ${NR_temp_peakslope} -eq 0 ]; then
	    printf "\e[34m We found no jump for ${input} with the NSIGMA=${nsigma} from Runid = ${firstrun} to Runid = ${lastrun}, will exit \e[39m \n"
	    rm  "Jump_total_${input}.dat" "file_temp_peak_slope_${input}.dat" 
	    rm "Peak_Slope_Points_jump_${input}.dat"
	    rm ${input}_raw.dat
	    rm "file_Goodruncut_sort"
	    rm Jump_smooth_${input}.dat
	    exit
	fi
	#check if the first line is the peak
	NR_firstpk=$(awk 'NR==1{if($6<=0){print(1)}else{print(0)}} ' "file_temp_peak_slope_${input}.dat")
	if [ $NR_firstpk -gt 0 ]; then
	    awk 'BEGIN {print(1, 2)} ' "file_temp_peak_slope_${input}.dat" > "temp_slope_peak_slope_NR.txt" 
	fi

	#find the BiggestSlope-Peak-BiggestSlope regions
	awk 'BEGIN {slope1=1;slope2=1;slope3=1;peak2=1;peak3=1;} {if (slope1<=0 && peak2<=0 && slope3<=0 ) print(NR-3, NR-1); slope1=slope2; slope2=slope3; peak2=peak3; slope3=$5; peak3=$6}' "file_temp_peak_slope_${input}.dat" >> "temp_slope_peak_slope_NR.txt"
	#check if the last line is the peak
	NR_lastpk=$(awk 'END {if($6<=0) {print(1)} else {print(0)}} ' "file_temp_peak_slope_${input}.dat")
	if [ $NR_lastpk -gt 0 ]; then
	    awk 'END {print(NR-1, NR)} ' "file_temp_peak_slope_${input}.dat" >> "temp_slope_peak_slope_NR.txt" 
	fi
#exit;

	if [ -f local_mean_sigma.txt ]; then
	    rm local_mean_sigma.txt
	fi
	if [ -f peakline.txt ]; then
	    rm peakline.txt
	fi
	
	NR_temp_slope_peak_slope=$(cat "temp_slope_peak_slope_NR.txt" | wc -l )
	
	if [ ${NR_temp_slope_peak_slope} -eq 0 ]; then
	    printf "\e[34m We found no jump for ${input} with the NSIGMA=${nsigma} from Runid = ${firstrun} to Runid = ${lastrun} \e[39m \n"
	    rm "file_Goodruncut_sort" "file_temp_peak_slope_${input}.dat" "temp_slope_peak_slope_NR.txt" "Peak_Slope_Points_jump_${input}.dat"
	    rm ${input}_raw.dat
	    rm Jump_smooth_${input}.dat Jump_total_${input}.dat
	    exit
	fi	    

	for i in `seq 1 $NR_temp_slope_peak_slope`
	do
	    tempnumber=$i
	    	  
	    lowlimit=$(awk 'NR=='${i}' {print($1)}' temp_slope_peak_slope_NR.txt)
	    highlimit=$(awk 'NR=='${i}' {print($2)}' temp_slope_peak_slope_NR.txt)
	    if [ $lowlimit -eq 1 ] && [ $highlimit -eq 2 ]; then
		    peakline=1
	    else
		peakline=$(($lowlimit+1))
	    fi

	    lowlimit_1=$(awk '(NR=='$lowlimit') {print($1)}' "file_temp_peak_slope_${input}.dat")
	    peakline_1=$(awk '(NR=='$peakline') {print($1)}' "file_temp_peak_slope_${input}.dat")
	    highlimit_1=$(awk '(NR=='$highlimit') {print($1)}' "file_temp_peak_slope_${input}.dat")

	    #match with the original data

	    if [ -f rawpeak_tempregion_${tempnumber}.dat ]; then
		rm rawpeak_tempregion_${tempnumber}.dat
	    fi
	    dpeaklow=$((${peakline_1}-${lowlimit_1}))
	    dpeakhigh=$((${highlimit_1}-${peakline_1}))
	    dpeakn10=$((${peakline_1}-5))
	    dpeakp10=$((${peakline_1}+5))

	    if [ ${dpeaklow} -gt 5 ]
	    then
		awk '($8>='$dpeakn10')&&($8<='${peakline_1}') {print($0)}' "file_Goodruncut_sort" > rawpeak_tempregion_${tempnumber}.dat
	    else
		if [ ${dpeaklow} -gt 0 ]
		then
		    awk '($8>='${lowlimit_1}')&&($8<='${peakline_1}') {print($0)}' "file_Goodruncut_sort" > rawpeak_tempregion_${tempnumber}.dat
		else
		    awk '($8>='${dpeakn10}')&&($8<='${peakline_1}') {print($0)}' "file_Goodruncut_sort" > rawpeak_tempregion_${tempnumber}.dat
		fi
	    fi
	    if [ ${dpeakhigh} -gt 5 ]
	    then
		awk '($8>'${peakline_1}')&&($8<='$dpeakp10') {print($0)}' "file_Goodruncut_sort" >> rawpeak_tempregion_${tempnumber}.dat
	    else
		if [ ${dpeakhigh} -gt 0 ]
		then
		awk '($8>'$peakline_1')&&($8<='$highlimit_1') {print($0)}' "file_Goodruncut_sort" >> rawpeak_tempregion_${tempnumber}.dat
		else
		awk '($8>'$peakline_1')&&($8<='$dpeakp10') {print($0)}' "file_Goodruncut_sort" >> rawpeak_tempregion_${tempnumber}.dat
		fi
	    fi
	    runnb=$(cat rawpeak_tempregion_${tempnumber}.dat | wc -l )	    
	    #calculate the local mean and sigma
	    #if we calculate by error of every run
	    #reference https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
	    
	    awk '{if ($3!=0){ u += $2*(1/($3*$3)); w += (1/($3*$3)); m = u/w; s += (1/($3*$3)); sig = (1/s)^0.5 }} END {print '$tempnumber',m, sig}' "rawpeak_tempregion_${tempnumber}.dat" >> local_mean_sigma.txt 
#	    awk '{if ($3!=0){ u += $2*(1/($3*$3)); w += (1/($3*$3)); m = u/w; s += (1/($3*$3)); sig = (1/s)^0.5 }} END {print '$tempnumber',m, sig*'$runnb'^0.5}' "rawpeak_tempregion_${tempnumber}.dat" >> local_mean_sigma.txt 

	    awk '(NR=='$peakline') {print($0)}' "file_temp_peak_slope_${input}.dat" >> peakline.txt
	    
	done

	#comebackhere
	#Gang suggests add |(u1-u2)/u1| > 1%
	awk '{ n1=n2; u1=u2; s1=s2; u2=$2; s2=$3; n2=$1 } {if( n1!=n2 ) {if( (NR>1)&&((u1-u2)^2-'$nsigma'*'$nsigma'*(s1^2+s2^2) >= 0.0 )&&(((u1-u2)/u1)^2 >= 0.0001)) print( n1, n2, (u1-u2)^2, (s1^2+s2^2), (u1-u2)^2-'$nsigma'*'$nsigma'*(s1^2+s2^2), ((u1-u2)/u1)^2-0.0001)}}' "local_mean_sigma.txt" > jump_in_tempregion.txt #give the jump place between two region numbers


	if [ -f Cut_for_${input}.dat ]; then
	    rm Cut_for_${input}.dat
	fi

	
	NR_peaknumber=$(cat "jump_in_tempregion.txt" | wc -l )
	echo ${firstrun} '#FirstRun ' >> Cut_for_${input}.dat

	if [ ${NR_peaknumber} -eq 0 ]; then
	    printf "\e[34m We found no jump for ${input} with the NSIGMA=${nsigma} from Runid = ${firstrun} to Runid = ${lastrun} \e[39m \n"
	    rm "file_temp_peak_slope_${input}.dat" "temp_slope_peak_slope_NR.txt" "jump_in_tempregion.txt"
	    rm rawpeak_tempregion_*.dat 	
	    rm ${input}_raw.dat
	    rm file_Goodruncut_sort
	    rm Jump_smooth_${input}.dat Jump_total_${input}.dat Peak_Slope_Points_jump_${input}.dat
	    rm peakline.txt
	    rm local_mean_sigma.txt
	    exit
	fi

	for i in `seq 1 $NR_peaknumber`
	do
	    if [ -f cutpoint_${i}.txt ]; then
		rm cutpoint_${i}.txt
	    fi
	    tempnumber=$i
	    lowlimit=$(awk 'NR=='${i}' {print($1)}' jump_in_tempregion.txt)
	    highlimit=$(($lowlimit+1))

	    lowlimit_1=$(awk 'NR=='${lowlimit}' {print($1+1)}' temp_slope_peak_slope_NR.txt)
	    highlimit_1=$(awk 'NR=='${highlimit}' {print($1+1)}' temp_slope_peak_slope_NR.txt)

	    awk ' {if((NR>'${lowlimit_1}')&&(NR<'${highlimit_1}')){ cutpoint1=cutpoint2; cutid1=cutid2; cutpoint2=($10)^2; cutid2=$2; if(cutpoint1<cutpoint2){cutpoint1=cutpoint2; cutid1=cutid2 } else {cutpoint2=cutpoint1; cutid2=cutid1}; if(NR-1=='${lowlimit_1}') print(cutpoint2^0.5, cutid2); else print(cutpoint1^0.5, cutid1) } }' "Peak_Slope_Points_jump_${input}.dat" > cutpoint_${i}.txt

	    awk 'END{print $2}' cutpoint_${i}.txt >> Cut_for_${input}.dat
	done	    
	echo ${lastrun} '#LastRun' >> Cut_for_${input}.dat
	
	printf "\e[34m We found ${NR_peaknumber} jump for ${input} with the NSIGMA=${nsigma} from Runid = ${firstrun} to Runid = ${lastrun} \e[39m \n"
	awk '{print($0)}' Cut_for_${input}.dat

	#cache file clean up
	rm "file_temp_peak_slope_${input}.dat" "temp_slope_peak_slope_NR.txt" "jump_in_tempregion.txt"
	rm rawpeak_tempregion_*.dat cutpoint_*.txt	
	rm ${input}_raw.dat
	rm file_Goodruncut_sort
	rm Jump_smooth_${input}.dat Jump_total_${input}.dat Peak_Slope_Points_jump_${input}.dat
	rm peakline.txt
	rm local_mean_sigma.txt
exit;
