#!/bin/bash

echo Which ocean is this?

read ocean

echo This is $ocean

for value in bnor cnor enor gnor nnor qnor snor nod; 
do
	awk -F ' ' '{print $9}' new_nors_HMMrun/HMMout/all_${ocean}_${value}.hits > new_nors_HMMrun/HMMout/all_${ocean}_${value}_IDs.txt
	seqtk subseq all_${ocean}_PROKKA.faa new_nors_HMMrun/HMMout/all_${ocean}_${value}_IDs.txt > new_nors_HMMrun/HMMout/all_${ocean}_${value}_seqs.faa
done
