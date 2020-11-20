#!/bin/tcsh
#daynum = today's date

set daynum=$1
set cent = $2
set opt_weight = $3
set lam = $4
set lum = $5
set OUTPUT="/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/output/${lam}_${lum}"

#set up output directory
echo ${OUTPUT}
if ( -e ${OUTPUT}/ ) then
echo "dir: ${OUTPUT}/ already exists"
else
mkdir ${OUTPUT}/
endif

#set up output directory
echo ${OUTPUT}
if ( -e ${OUTPUT}/Data_${daynum}_${cent}_debug ) then
echo "dir: ${OUTPUT}/Data_${daynum}_${cent}_debug already exists"
rm -r ${OUTPUT}/Data_${daynum}_${cent}_debug
endif
mkdir ${OUTPUT}/Data_${daynum}_${cent}_debug

#set up data directory
if ( -e ${OUTPUT}/Data_${daynum}_${cent} ) then
echo "dir: ${OUTPUT}/Data_${daynum}_${cent} already exists"
rm -r ${OUTPUT}/Data_${daynum}_${cent}
endif
mkdir ${OUTPUT}/Data_${daynum}_${cent}

#submit the Scheduler template
star-submit-template -template ./Analysis.xml -entities daynum=${daynum},cent=${cent},opt_weight=${opt_weight},lam=${lam},lum=${lum}

echo "submitting job day number->${daynum} cent->${cent} lum->${lum} lam->${lam} done"
