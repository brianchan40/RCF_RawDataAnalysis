#!/bin/bash

number=$3
number2=$5
number3=0
number4=$1

while [ $number -le $4 ]
do
    echo $number
    ((number3=$2-number))
    echo $number3

    hadd ./temp/Data_${number4}_${number2}/Data_${number3}/t_${number3}cen${number2}.gamma112_fullEP_eff_pT02_module.root ./temp/Data_${number4}_${number2}/Data_${number3}/sched*_${number3}*cen${number2}.gamma112_fullEP_eff_pT02_module.root
    hadd ./temp/Data_${number4}_${number2}/Data_${number3}/t_${number3}plam.root ./temp/Data_${number4}_${number2}/Data_${number3}/sched*_${number3}*plam.root
    ((number++))
done

