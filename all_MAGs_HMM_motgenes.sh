#!/bin/bash

echo Which ocean is this?

read ocean

echo This is $ocean

for value in aer cheA cheB cheR cheV cheY fliG fliM fliNY motAC motBD; 
do
	hmmsearch -o che_mot_genes_HMMrun/all_${ocean}_${value}.out ../HMM_profiles/chem_mot_HMMs/${value}.hmm all_${ocean}_PROKKA.faa
done


