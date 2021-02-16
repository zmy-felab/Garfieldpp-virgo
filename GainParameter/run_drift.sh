#!/bin/bash

# complie code
make CeramicGEM
if [ $? -ne 0 ]
then
    echo "make error!"
    exit
fi

# create directory
if [ ! -d "./result" ]
then
    mkdir result
fi
if [ ! -d "./runlog" ]
then
    mkdir runlog
fi

# number of events
nEvents=10000
# mixture gas
gas=ar_90_co2_10

# electric field in induction region
induce=(2.0 4.0 8.0)
for i in ${induce[@]}
do
    # electric field in drift region
    drift=(0.1 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 5.0 6.0 7.0 8.0)
    for j in ${drift[@]}
    do
        # voltage of GEM
        for k in 700 800 900 1000
        do
            echo D_${j}kV_I_${i}_GEM_${k}V
            # run and save log
            ./CeramicGEM -n $nEvents -v $k -d $j -i $i > "./runlog/${nEvents}_${k}_${j}_${i}_${gas}.log" &
        done
    done
done