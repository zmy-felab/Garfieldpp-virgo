#!/bin/bash
make
if [ $? -ne 0 ]; then
    echo "make error!"
    exit
fi

if [ ! -d "./result" ]; then
    mkdir result
fi
if [ ! -d "./runlog" ]; then
    mkdir runlog
fi

nEvents=1
voltage=740

for i in `seq 0 3`
do
    echo $i $nEvents ${voltage}V
    ./x-ray -n $nEvents -v $voltage > "./runlog/${nEvents}_${voltage}.log" &

    nEvents=$(($nEvents+1))
    # voltage=$(($voltage+20))
done
