#!/bin/bash

number=11
number2=0

while [ $number -le 11 ]
do
    while [ $number2 -le 4 ]
    do
        echo $number
        hadd ./temp/${number}cen${number2}.weight_112_module_new.root ./temp/Data_${number}_${number2}/t*cen${number2}.weight_112_module_new.root
        hadd ./temp/${number}cen${number2}.v2_fullEP_eff_pT02_module.root ./temp/Data_${number}_${number2}/t*cen${number2}.v2_fullEP_eff_pT02_module.root
        hadd ./temp/${number}cen${number2}plam.root ./temp/Data_${number}_${number2}/t*plam.root
        ((number2++))
    done
	((number2=0))
    ((number++))
done

