#!/bin/bash

if [ $# -ne 2 ]
then
	echo "usage $0 <soemedi-label-file> <data-with-skippy-file>"
	exit 1
fi

SOEMEDI_FILE=$1
INPUT_FILE=$2

# ID,ESM,mwSS_scoredif,HGMD_SSV_scoredifsum,ESE_sum_mexon,intronlength,dG_wexon,exonln,ESS_sum_mexon,chasin_wmdif,distance_closest_SS,scoresum_wexon,HI_imp_score,mean_exonphastCons,N_introns,exonpos,PPTscore,ExAC_SSV_scoredifsum,nucleotide_change,chrom,start,ref,alt,hg19_libpos
# chr1-12254011-12254089W,0,0,0,13,2.83681268176447,-252.7,78,29,0.0182407000000002,4,18.21,0.592,0.94,9,0.7,221,7.26,T>C,chr1,12254015,T,C,chr1-12253934-12254104

# lib_hg19,mut_coord,allele_w,allele_m,in_vivo_wu,in_vivo_mu,in_vivo_ws,in_vivo_ms,vit_wu,vit_mu,vit_ws,vit_ms,sequence,wt5start,wt5sequence,wt5score,mu5start,mu5sequence,mu5score,wt3start,wt3sequence,wt3score,mu3start,mu3sequence,mu3score,lor,exonlen,fivePrimeSSscore,threePrimeSSscore,deltaSS3,deltaSS5,n_ese,n_ess
#chr3:37083672-37083842,chr3:37083774-37083775,C,T,5515,780,385,10,5807,792,2042,115,ggattacttctcccattttgtcccaactggttgtatctcaagcatgaattcagcttttccttaaagtcacttcatttttattttcagtgaagaactgttctaccagatactcatttatgattttgccaattttggtgttctcaggttatcggtaagtttagatccttttca,37083775,Cagatactc,-1.1415193178520904,37083775,Tagatactc,-3.506761467131295,37083758,gtgaagaactgttctacCagata,2.4974727419772793,37083758,gtgaagaactgttctacTagata,1.4422214557905244,0,64,0.682,11.78,8.17,0.254,10,12

for line in $(sed -n '1,$ p' $INPUT_FILE)
do
        if [ $(echo $line | grep -c ^lib_hg19) != 0 ]
        then
		echo ${line},esm
                continue
        fi
	chrom=$(echo $line | awk -F: '{print $1}')
	startPos=$(echo $line | awk -F: '{print $2}' | cut -d, -f1)	
	endPos=$(echo $line | awk -F: '{print $2}' | cut -d- -f2 | cut -d, -f1)
	sline=$(grep $chrom $SOEMEDI_FILE | grep $startPos | grep $endPos)
	if [ -n "$sline" ]
	then
        	ESM=$(echo $sline | awk -F, '{print $2}')
        	echo ${line},${ESM}
	fi
done
