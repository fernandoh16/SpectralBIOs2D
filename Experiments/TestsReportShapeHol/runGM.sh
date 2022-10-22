#!/bin/bash

# Generate Measurements
# ../build/GenerateMeasurements --NumberOfParameters=100 --Gamma=1.0 --NDirections=4 --NControlPoints=200 --MeasFile=... --TrueBoundary=...

rm -r TrueBoundary
rm -r MeasurementsFolder

mkdir -p TrueBoundary 
mkdir -p MeasurementsFolder 

#G=(1.0 0.5 0.1)
G=(0.0)
#ND=(4 8 16)
ND=(32)
NCP=200 ;
NP=48 ;
M=5.0;
for g in "${G[@]}"; do
	for nd in "${ND[@]}"; do
		MF=MeasurementsFolder/Measurements_${g}_${nd}.txt ;
		TB=TrueBoundary/TrueBoundary_${g}_${nd}.txt ;
        ../../build/GenerateMeasurements --NumberOfParameters=${NP} --Gamma=${g} --NDirections=${nd} --NControlPoints=${NCP} --MeasFile=${MF} --TrueBoundary=${TB} --Magnitude=${M}
    done
done