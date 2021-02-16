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

# electric field in drift region
drift=(1.0 2.0 4.0)
for i in ${drift[@]}
do
    # electric field in induction region
    induce=(0.5 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0)
    for j in ${induce[@]}
    do
        # exclude D_1.0kV_I_3.0kV 
        if [ `echo $i!=1.0|bc` -eq 1 -o `echo $j!=3.0|bc` -eq 1 ]
        then
            # voltage of GEM
            for k in 700 800 900 1000
            do
                echo D_${i}kV_I_${j}kV_GEM_${k}V
                # run and save log
                ./CeramicGEM -n $nEvents -v $k -d $i -i $j > "./runlog/${nEvents}_${k}_${i}_${j}_${gas}.log" &
            done
        fi
    done
done