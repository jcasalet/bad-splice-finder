#!/bin/bash

if [ $# -ne 1 ]
then
	echo "usage $0 <data-with-skippy.csv>"
	exit 1
fi

INPUT_FILE=$1

echo lib_hg19,mut_coord,allele_w,allele_m,in_vivo_wu,in_vivo_mu,in_vivo_ws,in_vivo_ms,vit_wu,vit_mu,vit_ws,vit_ms,sequence,wt5start,wt5sequence,wt5score,mu5start,mu5sequence,mu5score,wt3start,wt3sequence,wt3score,mu3start,mu3sequence,mu3score,lor,exonlen,fivePrimeSSscore,threePrimeSSscore,deltaSS3,deltaSS5,n_ese,n_ess

# lib_hg19,mut_coord,allele_w,allele_m,in_vivo_wu,in_vivo_mu,in_vivo_ws,in_vivo_ms,vit_wu,vit_mu,vit_ws,vit_ms,sequence,wt5start,wt5sequence,wt5score,mu5start,mu5sequence,mu5score,wt3start,wt3sequence,wt3score,mu3start,mu3sequence,mu3score,lor,exonlen,fivePrimeSSscore,threePrimeSSscore,deltaSS3,deltaSS5
#chr3:37083672-37083842,chr3:37083774-37083775,C,T,5515,780,385,10,5807,792,2042,115,ggattacttctcccattttgtcccaactggttgtatctcaagcatgaattcagcttttccttaaagtcacttcatttttattttcagtgaagaactgttctaccagatactcatttatgattttgccaattttggtgttctcaggttatcggtaagtttagatccttttca,37083775,Cagatactc,-1.1415193178520904,37083775,Tagatactc,-3.506761467131295,37083758,gtgaagaactgttctacCagata,2.4974727419772793,37083758,gtgaagaactgttctacTagata,1.4422214557905244,0,64,0.682,11.78,8.17,0.254

for line in $(sed -n '1,$ p' $INPUT_FILE) 
do
	if [ $(echo $line | grep -c ^lib) != 0 ]
	then
		continue
	fi
	sequence=$(echo $line | awk -F, '{print $13}')
	ese_count=$(curl http://genes.mit.edu/cgi-bin/rescue-ese_new.pl -d Human=1 -d sequence=$sequence -d process=true 2>/dev/null | grep "total matches"  | cut -d\< -f1 | awk '{print $3}' | awk '{$1=$1};1')
	ess_count=$(curl http://genes.mit.edu/cgi-bin/fas-ess.pl -d set=FAS-hex2 -d sequence=$sequence 2>/dev/null | grep -o 'font color="red"' | wc -l | awk '{$1=$1};1')
	echo ${line},${ese_count},${ess_count}
done



# get the ESE score
#curl http://genes.mit.edu/cgi-bin/rescue-ese_new.pl -d Human=1 -d sequence=ggattacttctcccattttgtcccaactggttgtatctcaagcatgaattcagcttttccttaaagtcacttcatttttattttcagtgaagaactgttctaccagatactcatttatgattttgccaattttggtgttctcaggttatcggtaagtttagatccttttca -d process=true 2>/dev/null | grep "total matches"  | cut -d\< -f1 | awk '{print $3}' | awk '{$1=$1};1'

# get the ESS score
#curl http://genes.mit.edu/cgi-bin/fas-ess.pl -d set=FAS-hex2 -d sequence=ggattacttctcccattttgtcccaactggttgtatctcaagcatgaattcagcttttccttaaagtcacttcatttttattttcagtgaagaactgttctaccagatactcatttatgattttgccaattttggtgttctcaggttatcggtaagtttagatccttttca 2>/dev/null | grep -o 'font color="red"' | wc -l | awk '{$1=$1};1'

