#!/bin/bash

if [ ! -d "./result" ];then
    mkdir result
fi
if [ ! -d "./runlog" ];then
    mkdir runlog
fi

nEvents=10
voltage=700
pressure=760.0
temperature=293.15

for i in `seq 0 1`
do
    # let voltage=voltage+20
    pressure=`echo "$pressure-1."|bc`
    # temperature=`echo "$temperature+2."|bc`
    echo $i $nEvents $voltage $pressure $temperature

    ./CeramicGEM -n $nEvents -v $voltage -p pressure -t $temperature > "./runlog/${nEvents}_${voltage}_${pressure}_${temperature}.log" 
done
