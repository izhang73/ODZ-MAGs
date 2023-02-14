#!/bin/bash

echo Which ocean is this?

read ocean

echo This is $ocean

echo Which set is this?

read setname

echo This is $setname

for value in napA narG nirK nirS norB nosZ nxrA amoA_AOA amoA_AOB nrfA nifH; 
do
	hmmsearch -o ${ocean}_PROKKA/HMMout/all_${setname}_${ocean}_${value}.out ../HMM_profiles/${value}.hmm ${ocean}_PROKKA/all_${setname}_${ocean}_PROKKA.faa
done


