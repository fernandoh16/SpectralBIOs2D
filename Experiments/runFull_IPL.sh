#!/bin/bash
#
# Run All: Measurements + BI
#
Root=$1;
Exec=$2;
#
G=(1.0 0.5 0.1)
#G=(1);
ND=(4 8 16)
#ND=(1);
NCP=200 ;
NP=48 ;
M=1.0;
NM=25;
#
rm -r ${Root}
mkdir ${Root}
mkdir ${Root}/Measurements 
mkdir ${Root}/TrueBoundary 
mkdir ${Root}/Prior
mkdir ${Root}/Posterior
#
for g in "${G[@]}"; do
	for nd in "${ND[@]}"; do
		MF=${Root}/Measurements/Measurements_${g}_${nd}.txt ;
		TB=${Root}/TrueBoundary/TrueBoundary_${g}_${nd}.txt ;
        ${Exec}/GenerateMeasurements --NumberOfParameters=${NP} --Gamma=${g} --NDirections=${nd} --NControlPoints=${NCP} --MeasFile=${MF} --TrueBoundary=${TB} --Magnitude=${M} --NumberOfModes=${NM}
    done
done
#
Lmax=14 ;
#Lmax=0;
#
for L in `seq 0 ${Lmax}`; do
	for g in "${G[@]}"; do
		for nd in "${ND[@]}"; do
			Pri=${Root}/Prior/prior_${L}_${g}_${nd}.txt
			Pos=${Root}/Posterior/posterior_${L}_${g}_${nd}.txt
			MF=${Root}/Measurements/Measurements_${g}_${nd}.txt ;
			mpirun -np 2 ${Exec}/TestBI_IPL_Full --Level=${L} --NumberOfParameters=${NP} --Gamma=${g} --NDirections=${nd} --NControlPoints=${NCP}  --Prior=${Pri} --Posterior=${Pos} --Measurements=${MF} --Magnitude=${M}
		done
	done
done