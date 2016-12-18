#!/bin/bash
cd /ysm-gpfs/home/jss245/scratch/20161215aligntrim/ATout;
for f in Pol_*.txt
do
	div=`echo "Div_${f:4:9}"*`
	echo "source ~/.bashrc; cd /ysm-gpfs/home/jss245/scratch/20161215aligntrim/ATout; MASSPRF_preprocess -p $f -d $div -o 1 -ci_m 1 -s 1 -t 1 > /ysm-gpfs/home/jss245/scratch/20161215consensus/${f:4:9}_MACPRF_prep.txt" >> /ysm-gpfs/home/jss245/scratch/20161215consensus/jobs.list
done

