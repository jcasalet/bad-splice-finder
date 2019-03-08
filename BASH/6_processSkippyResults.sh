#!/bin/bash
  
if [ $# -ne 3 ]
then
        echo "usage: $0 <skippy-output> <data-with-maxentscan.csv> <problematic-results-out>"
        exit 1
fi

# print header
echo "lib_hg19,mut_coord,allele_w,allele_m,in_vivo_wu,in_vivo_mu,in_vivo_ws,in_vivo_ms,vit_wu,vit_mu,vit_ws,vit_ms,sequence,wt5start,wt5sequence,wt5score,mu5start,mu5sequence,mu5score,wt3start,wt3sequence,wt3score,mu3start,mu3sequence,mu3score,lor,exonlen,fivePrimeSSscore,threePrimeSSscore,deltaSS3,deltaSS5"

while IFS='' read -r line || [[ -n "$line" ]]; do
        lor=$(echo $line | awk '{print $17}')
	exonlen=$(echo $line | awk '{print $10}')
	fivePrimeSSscore=$(echo $line | awk '{print $23}')
	threePrimeSSscore=$(echo $line | awk '{print $24}')
	deltaSS3=$(echo $line | awk '{print $25}')
	deltaSS5=$(echo $line | awk '{print $28}')
	chr=$(echo $line | awk '{print $1}')
	loc=$(echo $line | awk '{print $2}')
	((locminus1=$loc - 1))
	mut=$(echo $line | awk '{print $3}')
	x=$(echo $mut | cut -d- -f1)
	y=$(echo $mut | cut -d\> -f2)
	searchMaxEntString="$chr:$locminus1-$loc,$x,$y"
	# chr20   49565187        C->G 
	searchSkippyString="$chr\t$loc\t$x->$y"	
	if [ $(grep -c $searchSkippyString $1) -eq 0 ]
	then
		echo "not found in $1: $searchSkippyString" >> $3
	elif [ $(grep -c $searchSkippyString $1) -gt 1 ]
	then
		echo "found_multiple_times_in:${1}:${searchSkippyString}" >> $3
	else
		thisLine=$(grep $searchMaxEntString $2)
		echo ${thisLine},$lor,$exonlen,$fivePrimeSSscore,$threePrimeSSscore,$deltaSS3,$deltaSS5

	fi	 
done < "$1"

# maxentscan header 
#lib_hg19,mut_coord,allele_w,allele_m,in_vivo_wu,in_vivo_mu,in_vivo_ws,in_vivo_ms,vit_wu,vit_mu,vit_ws,vit_ms,sequence,wt5start,wt5sequence,wt5score,mu5start,mu5sequence,mu5score,wt3start,wt3sequence,wt3score,mu3start,mu3sequence,mu3score

# maxentscan sample line
#chr3:37083672-37083842,chr3:37083774-37083775,C,T,5515,780,385,10,5807,792,2042,115,ttttgtcccaactggttgtatctcaagcatgaattcagcttttccttaaagtcacttcatttttattttcagtgaagaactgttctaccagatactcattt,37083672,None,-99.0,37083672,None,-99.0,37083672,None,-99.0,37083672,None,-99.0

# skippy header

# skippy sample line
#chr3    37083775        C->T    internal_exon_protein_coding    MLH1    ENSE00001747618 37083759        37083822        1       64      2       0       1       0       2       1       0       1.73    17      0.5312  0.2656  1.546   0.682   11.78   8.17    0       0       0       0       0       0       0.136   -1.818  0.254   2.037   0.253   0.023   0.400   1.271

# output header from this script (30 fields)
#lib_hg19,mut_coord,allele_w,allele_m,in_vivo_wu,in_vivo_mu,in_vivo_ws,in_vivo_ms,vit_wu,vit_mu,vit_ws,vit_ms,sequence,wt5start,wt5sequence,wt5score,mu5start,mu5sequence,mu5score,wt3start,wt3sequence,wt3score,mu3start,mu3sequence,mu3score,lor,exonlen,fivePrimeSSscore,threePrimeSSscore,deltaSS3,deltaSS5

# sample output from this script (30 fields)
#chr2:189854936-189855106,chr2:189855033-189855034,G,A,2345,322,3520,484,2535,258,9624,1388,acttaaaaacagaaagtgttttactactagattgtgattctatttgaaggttcattaatattttttcattcattatttttagggtatcaaaggtccagctg,189855025,aaggtccag,-1.3012874202034164,189855025,aaggtccaa,-1.9967689952451928,189855014,ttagggtatcaaaggtccagctg,-6.214633523023855,189855014,ttagggtatcaaaggtccaactg,-14.965012826743692,0,54,0.591,9.35,10.55,1
