rm -r TrueBoundary
rm -r MeasurementsFolder

mkdir -p TrueBoundary 
mkdir -p MeasurementsFolder 

G=(0.1)
ND=(32)
NCP=200;
NP=100;

for g in "${G[@]}"; do
	for nd in "${ND[@]}"; do
		MF=MeasurementsFolder/Measurements_${g}_${nd}.txt ;
		TB=TrueBoundary/TrueBoundary_${g}_${nd}.txt ;
        ../build/GenerateMeasurements --NumberOfParameters=${NP} --Gamma=${g} --NDirections=${nd} --NControlPoints=${NCP} --MeasFile=${MF} --TrueBoundary=${TB}
    done
done

Lmax=14 ;
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
			mpirun -np 100 ../build/TestBI_IPL_Full --Level=${L} --NumberOfParameters=${NP} --Gamma=${g} --NDirections=${nd} --NControlPoints=${NCP} --NumberOfModes=100 --Prior=${Pri} --Posterior=${Pos} --Measurements=${MF}
		done
	done
done