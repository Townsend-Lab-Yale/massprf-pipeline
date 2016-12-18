#!/bin/bash
targetscript="extractGenePolymorphism.py"

arg1="-gn"

arg2="-pm"
param2="20161201/ellisonSNP.txt"

arg3="-gm"
param3="20161201/neurospora_crassa_OR74A_10.fasta"


outfiles="*_polymorphs.txt"
outdir="20161201/output/"
for f in $@
do
	echo "Processing $f.."
	python "$targetscript" "$arg1" "$f" "$arg2" "$param2" "$arg3" "$param3"
done

mkdir $outdir
mv $outfiles $outdir