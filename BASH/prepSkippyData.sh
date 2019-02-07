#!/bin/bash

if [ $# -ne 1 ]
then
	echo "usage: $0 <input>" 
	exit 1
fi

INPUT_FILE=$1

if [ ! -f $INPUT_FILE ]
then
	echo "file $INPUT_FILE doesn't exist"
	exit 2
fi

# final-data.csv
# "","lib_hg19","sequence","mut_coord","allele_w","allele_m","wt5score","mu5score","wt3score","mu3score"

# skippy-input.txt
# chr3 37083775 C T

for line in $(sed -n '1,$ p' $INPUT_FILE)
do
	#echo $line
	chr=$(echo $line | awk -F, '{print $2}' | cut -d: -f1 | cut -d\" -f2)
	pos=$(echo $line | awk -F, '{print $4}' | cut -d: -f2 | cut -d- -f1)
	w=$(echo $line | awk -F, '{print $5}' | cut -d\" -f2)
	m=$(echo $line | awk -F, '{print $6}' | cut -d\" -f2)
	echo $chr $pos $w $m
done


