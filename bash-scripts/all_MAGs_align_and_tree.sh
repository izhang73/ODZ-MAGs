#!/bin/bash

#From HMM hits for denitrification genes, get protein sequences, align them with MAFFT, and make trees if needed for each metagenome
echo Which metagenome is this?

read varname

echo This is $varname


for value in napA narG nirK nirS norB nosZ; 
do
	awk -F ' ' '{print $9}' HMMout/all_${varname}_MAGs_${value}.hits > HMMout/all_${varname}_MAGs_${value}_IDs.txt
	seqtk subseq ../all_${varname}_MAGs_over70_cat.faa all_${varname}_MAGs_${value}_IDs.txt > Sequences/all_${varname}_MAGs_${value}_seqs.faa
	mafft-linsi --leavegappyregion Sequences/all_${varname}_MAGs_${value}_seqs.faa > Sequences/Alignments/all_${varname}_MAGs_${value}_aligned.faa
done

for value in napA narG nirK nirS norB nosZ; 
do
	/Users/irene/Desktop/tools/trimAl/source/trimal -in Sequences/Alignments/all_${varname}_MAGs_${value}_aligned.faa -automated1 -htmlout Sequences/Alignments/trimmed/auto_${value}_from_${varname}.html -out Sequences/Alignments/trimmed/all_${varname}_MAGs_${value}_MAFFT_trimal_auto.faa
	/Users/irene/Desktop/tools/iqtree-1.6.12-MacOSX/bin/iqtree -s Sequences/Alignments/trimmed/all_${varname}_MAGs_${value}_MAFFT_trimal_auto.faa -bb 1000 -mset Dayhoff,JTT,LG,WAG
done
