#!/bin/bash

#run HMM searching using new nor HMM profiles on MAGs from any ocean
echo Which ocean is this?

read ocean

echo This is $ocean

for value in bnor cnor enor gnor nnor qnor snor nod; 
do
	hmmsearch -o new_nors_HMMrun/all_${ocean}_${value}.hits ../HMM_profiles/nor_HMMs/${value}.hmm all_${ocean}_PROKKA.faa
done


