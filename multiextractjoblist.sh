#!/bin/bash

input=/ysm-gpfs/home/jss245/scratch/20161207extractPoly/ortho/output/*_polymap.txt

for f in $input;
do
	echo "cd /ysm-gpfs/home/jss245/scratch/20161207extractPoly"";"" ./extractPoly.py -gn" $f "-pm ./ellisonSNP.txt -gm ./ncr1_assembly.fasta -o ./pmOut/ \n"> jobs.list 
done    