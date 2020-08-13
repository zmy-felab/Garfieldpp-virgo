#!/bin/bash

if [ ! -d "./result" ];then
    mkdir result
fi
if [ ! -d "./runlog" ];then
    mkdir runlog
fi

voltage=760
pressure=760.0
temperature=293.15

for i in `seq 1 5`
do
    let voltage=voltage+20
    let presure+=1
    presure=`expr $presure + 1`
    temperature=`echo "$temperature+0.5"|bc`
    echo $i $voltage $presure $temperature
done




