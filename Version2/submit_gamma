#!/bin/tcsh
set daynum=$1
set cent=$2
set first=$3
set endnum=$4
set max=$5
set OUTPUT="/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/output/combine"

#%@first=$3*5
#endnum= 'expr $beginnum \* 5 - 1'

#set up output directory
echo ${OUTPUT}
if ( -e ${OUTPUT}/ ) then
#echo "dir: ${OUTPUT}/ already exists"
else
mkdir ${OUTPUT}/
endif

#set up output directory
echo ${OUTPUT}
if ( -e ${OUTPUT}/Data_${daynum}_${cent}_debug ) then
#echo "dir: ${OUTPUT}/Data_${daynum}_${cent}_debug already exists"
#rm -r ${OUTPUT}/Data_${daynum}_${cent}_debug
else
mkdir ${OUTPUT}/Data_${daynum}_${cent}_debug
endif

#set up data directory
if ( -e ${OUTPUT}/Data_${daynum}_${cent} ) then
#echo "dir: ${OUTPUT}/Data_${daynum}_${cent} already exists"
#rm -r ${OUTPUT}/Data_${daynum}_${cent}
else
mkdir ${OUTPUT}/Data_${daynum}_${cent}
endif

#submit the Scheduler template
star-submit-template -template ./Combine_g.xml -entities daynum=${daynum},cent=${cent},first=${first},endnum=${endnum},max=${max}

echo "submitting job day number->${daynum} cent->${cent} first->${first} done"
