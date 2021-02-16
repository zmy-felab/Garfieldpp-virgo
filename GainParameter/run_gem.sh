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
# start voltage of GEM
voltage=700
# electric field in drift region
driftE=1.0
# electric field in induction region
inductionE=3.0
# mixture gas
gas=ar_90_co2_10

for i in `seq 0 20`
do
    echo $i $nEvents ${voltage}V ${driftE}kV/cm ${inductionE}kV/cm
    # run and save log
    ./CeramicGEM -n $nEvents -v $voltage > "./runlog/${nEvents}_${voltage}_${driftE}_${inductionE}_${gas}.log" &

    voltage=$(($voltage+20))
done