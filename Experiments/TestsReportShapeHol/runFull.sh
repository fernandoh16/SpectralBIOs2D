#!/bin/bash
G=(0.1);
ND=(16);
NCP=200;
Lmax=(1) ;
NP=100;
NM=25;

rm -r Prior
rm -r Posterior
mkdir -p Prior
mkdir -p Posterior
for L in "${Lmax[@]}"; do
	for g in "${G[@]}"; do
		for nd in "${ND[@]}"; do
			Pri=Prior/prior_${L}_${g}_${nd}.txt
			Pos=Posterior/posterior_${L}_${g}_${nd}.txt
			MF=MeasurementsFolder/Measurements_${g}_${nd}.txt ;
			# echo ${MF}
			mpirun -np 2 ../build/TestBI_Halton_Full --Level=${L} --NumberOfParameters=${NP} --Gamma=${g} --NDirections=${nd} --NControlPoints=${NCP} --NumberOfModes=${NM} --Prior=${Pri} --Posterior=${Pos} --Measurements=${MF} --Magnitude=${M}
		done
	done
done