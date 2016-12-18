#!/bin/bash



for f in *_polymorphs_*
do
	string=*_${f:0:8}_*
	mv $string Div_${f:0:8}_`expr $string : '\(.[NtLA]*.[0-9]*\)'`${f:19};
	mv $f Pol_${f:0:8}${f:19};
done;