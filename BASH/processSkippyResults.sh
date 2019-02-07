#!/bin/bash
  
if [ $# -ne 3 ]
then
        echo "usage: $0 <skippy-in> <soemedi-in> <multi-out>"
        exit 1
fi

while IFS='' read -r line || [[ -n "$line" ]]; do
        lor=$(echo $line | awk '{print $17}')
	exonlen=$(echo $line | awk '{print $10}')
	fivePrimeSSscore=$(echo $line | awk '{print $23}')
	threePrimeSSscore=$(echo $line | awk '{print $24}')
	deltaSS3=$(echo $line | awk '{print $25}')
	deltaSS5=$(echo $line | awk '{print $28}')
	chr=$(echo $line | awk '{print $1}')
	loc=$(echo $line | awk '{print $2}')
	((locplus1=$loc + 1))
	mut=$(echo $line | awk '{print $3}')
	x=$(echo $mut | cut -d- -f1)
	y=$(echo $mut | cut -d\> -f2)
	searchString="$chr:$loc-$locplus1,$x,$y"
	if [ $(grep -c $searchString $2) -gt 1 ]
	then
		echo "$searchString" >> $3
	else
		#sed -ne "/'"$searchString"'/ s/$/\,'"$lor"'/" $2	
		thisLine=$(grep $searchString $2)
		echo ${thisLine},$lor,$exonlen,$fivePrimeSSscore,$threePrimeSSscore,$deltaSS3,$deltaSS5

	fi	 
done < "$1"


#chr3    37083774        C->T    internal_exon_protein_coding    LRRFIP2 ENSE00000760648 37083636        37083806        -1      171     4       0       2       0       2       0       2       3.934   33      0.3860  0.1930  0       0       8.99    7.94    0       0       0       0       0       0       0.530   1.952   0.190   1.741   0.347   0.787   0.337   0.991
#chr3    37083774        C->T    first_last_exon_protein_coding  LRRFIP2 ENSE00001856842 37083645        37083806        -1      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA
#chr20   49565186        C->G    non_exonic      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA



# lib_hg19,mut_coord,allele_w,allele_m,in_vivo_wu,in_vivo_mu,in_vivo_ws,in_vivo_ms,vit_wu,vit_mu,vit_ws,vit_ms,sequence,wt5start,wt5sequence,wt5score,mu5start,mu5sequence,mu5score,wt3start,wt3sequence,wt3score,mu3start,mu3sequence,mu3score
#chr3:37083672-37083842,chr3:37083774-37083775,C,T,5515,780,385,10,5807,792,2042,115,ggattacttctcccattttgtcccaactggttgtatctcaagcatgaattcagcttttccttaaagtcacttcatttttattttcagtgaagaactgttctaccagatactcatttatgattttgccaattttggtgttctcaggttatcggtaagtttagatccttttca,37083768,gttctacca,-14.7144839217625,37083766,ctgttctat,-12.660073394558811,37083758,gtgaagaactgttctaccagata,2.4974727419772793,37083758,gtgaagaactgttctatcagata,2.5474478513497334

# 17. lortotal - a log-odds ratio score that intergrates the three types of change (ESE loss, ESE gain and ESS gain) into a score that reflects how likely the combination of changes is for a splice-affecting variant. i.e. high positive scores are indicative of splice-affecting variants and negative scores are more indicative of splice-neutral SNPs.
