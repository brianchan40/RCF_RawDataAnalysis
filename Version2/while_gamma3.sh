#!/bin/bash

number=11
number2=1

#while [ $number -le 09 ]
#do
    while [ $number2 -le 8 ]
    do
    echo $number
    hadd ./temp/${number}cen${number2}.gamma112_fullEP_eff_pT02_module.root ./temp/Data_${number}_${number2}/t*cen${number2}.gamma112_fullEP_eff_pT02_module.root
    hadd ./temp/${number}cen${number2}plam.root ./temp/Data_${number}_${number2}/t*plam.root
    ((number2++))
done
#((number2=0))
#((number++))
#done
