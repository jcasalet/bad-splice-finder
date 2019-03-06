#!/bin/bash 

if [ $# -ne 1 ]
then
	echo "usage: $0 <data-with-sequences.csv> "
	exit 1
fi

INPUT=$1

../PYTHON/MyMaxentScan.py $INPUT 
