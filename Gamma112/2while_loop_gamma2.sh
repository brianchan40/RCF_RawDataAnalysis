#!/bin/bash

number=0
number2=$1
number3=0

while [ $number -le 19 ]
do
    while [[ $number2 -le $1 ]]
    do
        echo $number
        ((number3=19-number))
        echo $number3
        hadd ./temp2/Gamma112_${number2}/t${number3}cen${number2}.gamma112_fullEP_eff_pT02_module.root ./temp2/Gamma112_${number2}/${number3}*cen${number2}.gamma112_fullEP_eff_pT02_module.root
        rm ./temp2/Gamma112_${number2}/${number3}*cen${number2}.gamma112_fullEP_eff_pT02_module.root
        ((number2++))
    done
    number2=$1
    ((number++))
done

