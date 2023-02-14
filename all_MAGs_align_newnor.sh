#!/bin/bash

echo Which ocean is this?

read ocean

echo This is $ocean

for value in bnor cnor enor gnor nnor qnor snor nod; 
do
	mafft-linsi --leavegappyregion new_nors_HMMrun/HMMout/all_${ocean}_${value}_seqs.faa > new_nors_HMMrun/HMMout/all_${ocean}_${value}_aligned.faa
done