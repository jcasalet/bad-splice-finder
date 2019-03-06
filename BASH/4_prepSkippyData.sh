#!/bin/bash

#POSITION=0_NUMBERING
POSITION=1_NUMBERING

if [ $# -ne 1 ]
then
	echo "usage: $0 <data-with-maxentscan.csv>" 
	exit 1
fi

INPUT_FILE=$1

if [ ! -f $INPUT_FILE ]
then
	echo "file $INPUT_FILE doesn't exist"
	exit 2
fi

# input 
# "","lib_hg19","sequence","mut_coord","allele_w","allele_m","wt5score","mu5score","wt3score","mu3score"
# lib_hg19,mut_coord,allele_w,allele_m,in_vivo_wu,in_vivo_mu,in_vivo_ws,in_vivo_ms,vit_wu,vit_mu,vit_ws,vit_ms,sequence,wt5start,wt5sequence,wt5score,mu5start,mu5sequence,mu5score,wt3start,wt3sequence,wt3score,mu3start,mu3sequence,mu3score
#chr3:37083672-37083842,chr3:37083774-37083775,C,T,5515,780,385,10,5807,792,2042,115,ggattacttctcccattttgtcccaactggttgtatctcaagcatgaattcagcttttccttaaagtcacttcatttttattttcagtgaagaactgttctaccagatactcatttatgattttgccaattttggtgttctcaggttatcggtaagtttagatccttttca,37083768,gttctacca,-14.7144839217625,37083766,ctgttctat,-12.660073394558811,37083758,gtgaagaactgttctaccagata,2.4974727419772793,37083758,gtgaagaactgttctatcagata,2.5474478513497334

# skippy-input.txt
# chr3 37083775 C T

for line in $(sed -n '1,$ p' $INPUT_FILE)
do
	#echo $line
	chr=$(echo $line | awk -F, '{print $2}' | cut -d: -f1 | cut -d\" -f2)
	if [ $POSITION == "1_NUMBERING" ]
	then
		pos=$(echo $line | awk -F, '{print $2}' | cut -d: -f2 | cut -d- -f2)
	else
		pos=$(echo $line | awk -F, '{print $2}' | cut -d: -f2 | cut -d- -f1)
	fi
		
	w=$(echo $line | awk -F, '{print $3}' | cut -d\" -f2)
	m=$(echo $line | awk -F, '{print $4}' | cut -d\" -f2)
	echo $chr $pos $w $m
done


