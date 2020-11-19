#!/bin/tcsh
#daynum = today's date

set cent = $1
set opt_weight = $2
set lam = $3
set lum = $4
set OUTPUT="/star/u/brian40/PicoDst/PicoDst_NoMaker/Gamma112/output/"

#set up output directory
echo ${OUTPUT}
if ( -e ${OUTPUT}/${lam}_${lum} ) then
echo "dir: ${OUTPUT}/${lam}_${lum} already exists"
else
mkdir ${OUTPUT}/${lam}_${lum}
endif

#set up output directory
echo ${OUTPUT}
if ( -e ${OUTPUT}/${lam}_${lum}/Gamma112_${cent}_debug ) then
echo "dir: ${OUTPUT}/${lam}_${lum}/Gamma112_${cent}_debug already exists"
rm -r ${OUTPUT}/${lam}_${lum}/Gamma112_${cent}_debug
endif
mkdir ${OUTPUT}/${lam}_${lum}/Gamma112_${cent}_debug

#set up data directory
if ( -e ${OUTPUT}/${lam}_${lum}/Gamma112_${cent} ) then
echo "dir: ${OUTPUT}/${lam}_${lum}/Gamma112_${cent} already exists"
rm -r ${OUTPUT}/${lam}_${lum}/Gamma112_${cent}
endif
mkdir ${OUTPUT}/${lam}_${lum}/Gamma112_${cent}

#submit the Scheduler template
star-submit-template -template ./Analysis_v2.xml -entities cent=${cent},opt_weight=${opt_weight},lam=${lam},lum=${lum}

echo "submitting job day cent->${cent} lum->${lum} lam->${lam} done"
