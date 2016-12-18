#!/bin/bash

input=~/documents/school/grad/townsend/ntetncras/input/*

for f in $input;
do
	( ./splicer.py -csv $f -gf /Users/jaystanley/Documents/School/Grad/townsend/ntetncras/ncr1_genes.gff -fas /Users/jaystanley/Documents/School/Grad/townsend/ntetncras/ntet_cds.fasta & ); 
done    