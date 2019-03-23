#!/bin/bash
  
if [ $# -ne 2 ]
then
        echo "usage: $0 <data.csv> <path-to-jar>" 
        exit 1
fi

DATA_FILE=$1
JAR_FILE=$2


# print header
echo "lib_hg19,mut_coord,allele_w,allele_m,in_vivo_wu,in_vivo_mu,in_vivo_ws,in_vivo_ms,vit_wu,vit_mu,vit_ws,vit_ms,sequence,wt5start,wt5sequence,wt5score,mu5start,mu5sequence,mu5score,wt3start,wt3sequence,wt3score,mu3start,mu3sequence,mu3score,lor,exonlen,fivePrimeSSscore,threePrimeSSscore,deltaSS3,deltaSS5,n_ese,n_ess,esm,fisher" 

while IFS='' read -r line || [[ -n "$line" ]]
do
	if [ $(echo $line | grep -c ^lib) != 0 ]
	then
		continue
	fi
        a=$(echo $line | awk -F, '{print $5}')
        c=$(echo $line | awk -F, '{print $6}')
        b=$(echo $line | awk -F, '{print $7}')
        d=$(echo $line | awk -F, '{print $8}')
	#echo $a $b $c $d
	fisher=$(java -cp $JAR_FILE spliceutil.FisherExact $a $b $c $d)
	echo ${line},$fisher
	
done < "$DATA_FILE"
