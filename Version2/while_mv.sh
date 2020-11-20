#!/bin/bash

number=0
number2=5
number3=0
number4=$1

while [ $number2 -le 8 ]
do
    while [ $number -le $2 ]
    do
        echo $number
        ((number3=$2-number))
        echo $number3

        mkdir ./temp/Data_${number4}_${number2}/Data_${number3}
        mv ./temp/Data_${number4}_${number2}/sched*_${number3}*cen* ./temp/Data_${number4}_${number2}/Data_${number3}
        mv ./temp/Data_${number4}_${number2}/sched*_${number3}*plam.root ./temp/Data_${number4}_${number2}/Data_${number3}

        ((number++))
    done
    ((number2++))
    ((number=0))
#echo $number2
done

