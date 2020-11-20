#!/bin/bash

number=0
number2=$1
number3=0

while [ $number -le 16 ]
do
    while [[ $number2 -le $1 ]]
    do
        echo $number
        ((number3=16-number))
        echo $number3
        hadd ./temp/Gamma112_${number2}/t${number3}cen${number2}.weight_112_module_new.root ./temp/Gamma112_${number2}/${number3}*cen${number2}.weight_112_module_new.root
        hadd ./temp/Gamma112_${number2}/t${number3}cen${number2}.v2_fullEP_eff_pT02_module.root ./temp/Gamma112_${number2}/${number3}*cen${number2}.v2_fullEP_eff_pT02_module.root
        rm ./temp/Gamma112_${number2}/${number3}*cen${number2}.weight_112_module_new.root
        rm ./temp/Gamma112_${number2}/${number3}*cen${number2}.v2_fullEP_eff_pT02_module.root
        ((number2++))
    done
    number2=$1
    ((number++))
done

