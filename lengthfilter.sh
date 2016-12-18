#!/bin/bash

for f in *.txt
do
	len=`grep -Po ".*?\\d+.*?\\d+" <<< $f | tail -n1`;
	len=${len:1};
	if [ $len < 20 ]; then
		echo $f;
	fi
done