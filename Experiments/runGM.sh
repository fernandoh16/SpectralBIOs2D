#!/bin/bash

# Generate Measurements
# ../build/GenerateMeasurements --NumberOfParameters=100 --Gamma=1.0 --NDirections=4 --NControlPoints=200 --MeasFile=... --TrueBoundary=...

rm -r TrueBoundary
rm -r MeasurementsFolder

mkdir -p TrueBoundary 
mkdir -p MeasurementsFolder 

G=(0.1)
ND=(16)
NCP=200;
NP=50;

for g in "${G[@]}"; do
	for nd in "${ND[@]}"; do
		MF=MeasurementsFolder/Measurements_${g}_${nd}.txt ;
		TB=TrueBoundary/TrueBoundary_${g}_${nd}.txt ;
        ../build/GenerateMeasurements --NumberOfParameters=${NP} --Gamma=${g} --NDirections=${nd} --NControlPoints=${NCP} --MeasFile=${MF} --TrueBoundary=${TB}
    done
done