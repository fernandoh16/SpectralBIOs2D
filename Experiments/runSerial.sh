#!/bin/bash
G=(0.1)
ND=(8)
NCP=100;
Lmax=3 ;
NP=50;
rm -r Prior
rm -r Posterior
mkdir -p Prior
mkdir -p Posterior
for L in `seq 0 ${Lmax}`; do
	for g in "${G[@]}"; do
		for nd in "${ND[@]}"; do
			Pri=Prior/prior_${L}_${g}_${nd}.txt
			Pos=Posterior/posterior_${L}_${g}_${nd}.txt
			MF=MeasurementsFolder/Measurements_${g}_${nd}.txt
			echo ${MF}
			../build/TestBI_IPL_Serial --Level=${L} --NumberOfParameters=${NP} --Gamma=${g} --NDirections=${nd} --NControlPoints=${NCP} --NumberOfModes=50 --Prior=${Pri} --Posterior=${Pos} --Measurements=${MF}
		done
	done
done