#!/bin/bash

number=$3
number2=$5
number3=0
number4=$1

mkdir ./temp/Data_${number4}_${number2}/done_folders/

while [ $number -le $4 ]
do
    echo $number
    ((number3=$2-number))
    echo $number3
    hadd ./temp/Data_${number4}_${number2}/t${number3}cen${number2}.weight_112_module_new.root ./temp/Data_${number4}_${number2}/Data_${number3}*/t_${number3}*cen${number2}.weight_112_module_new.root
    hadd ./temp/Data_${number4}_${number2}/t${number3}cen${number2}.v2_fullEP_eff_pT02_module.root ./temp/Data_${number4}_${number2}/Data_${number3}*/t_${number3}*cen${number2}.v2_fullEP_eff_pT02_module.root
    hadd ./temp/Data_${number4}_${number2}/t${number3}cen${number2}plam.root ./temp/Data_${number4}_${number2}/Data_${number3}*/t_${number3}*plam.root
    mv ./temp/Data_${number4}_${number2}/Data_${number3}*/ ./temp/Data_${number4}_${number2}/done_folders/
    ((number++))
done

