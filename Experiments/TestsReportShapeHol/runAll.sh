#!/bin/bash

# Generate Measurements
# ../build/GenerateMeasurements --NumberOfParameters=100 --Gamma=1.0 --NDirections=4 --NControlPoints=200 --MeasFile=... --TrueBoundary=...

rm -r TrueBoundary
rm -r MeasurementsFolder

mkdir -p TrueBoundary 
mkdir -p MeasurementsFolder 

G=(1.0 0.5 0.1)
ND=(4 8 16)
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

NM=25;
Lmax=16 ;
rm -r Prior
rm -r Posterior
mkdir -p Prior
mkdir -p Posterior
for L in `seq 0 ${Lmax}`; do
	for g in "${G[@]}"; do
		for nd in "${ND[@]}"; do
			Pri=Prior/prior_${L}_${g}_${nd}.txt
			Pos=Posterior/posterior_${L}_${g}_${nd}.txt
			MF=MeasurementsFolder/Measurements_${g}_${nd}.txt ;
			# echo ${MF}
			mpirun -np 2 ../../build/TestBI_Halton_Full --Level=${L} --NumberOfParameters=${NP} --Gamma=${g} --NDirections=${nd} --NControlPoints=${NCP} --NumberOfModes=${NM} --Prior=${Pri} --Posterior=${Pos} --Measurements=${MF} --Magnitude=${M}
		done
	done
done