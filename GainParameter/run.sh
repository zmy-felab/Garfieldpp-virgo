#!/bin/bash
make CeramicGEM
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

nEvents=10
pressure=760.0
temperature=293.15
voltage=700
driftE=1.0
inductionE=3.0
rim=80

for i in `seq 0 5`
do
    echo $i $nEvents ${pressure}Torr ${temperature}K ${voltage}V ${driftE}kV/cm ${inductionE}kV/cm ${rim}Rim
    ./CeramicGEM -n $nEvents -v $voltage > "./runlog/${nEvents}_${pressure}_${temperature}_${voltage}_${driftE}_${inductionE}_${rim}.log" &

    voltage=$(($voltage+20))
    # pressure=`echo "$pressure-1."|bc`
    # temperature=`echo "$temperature+2."|bc`
    # rim=`expr $rim + 10`
done
