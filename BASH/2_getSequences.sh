#!/bin/bash

# "upstream" = 5 prime end
# "downstream" = 3 prime end
# Mutations were mapped to internal exons ≤ 100 nucleotides in length and selected for those that fit into 170 nucleotide genomic windows, which include 15 nucleotides of downstream intronic sequence and ≥ 55 nucleotides of upstream intronic sequence
PAD_5PRIME=55
PAD_3PRIME=15
#PAD_5PRIME=0
#PAD_3PRIME=0



INPUTFILE=../DATA/data.csv

OUTPUTFILE=../DATA/data-with-sequences.csv

if [ ! -f $INPUTFILE ]
then
	echo "input file $INPUTFILE does not exist - exiting!"
	exit 1
elif [ -f $OUTPUTFILE ]
then
	echo "output file $OUTPUTFILE already exists - exiting!"
	exit 2
fi


# chrX-54837931-54838101,chrX:54838076-54838077,C,T,398,218,9,0,1180,623,24,4

# wget -O /tmp/bob http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr3:37083672,37083842

# <?xml version="1.0" standalone="no"?>
#<!DOCTYPE DASDNA SYSTEM "http://www.biodas.org/dtd/dasdna.dtd">
#<DASDNA>
#<SEQUENCE id="chr2" start="189855633" stop="189855803" version="1.00">
#<DNA length="171">
#tatgaagaccaattagaaaaataccatgtaaactttatcaatcattctag
#attattaacagattttaataattttgctggttttatacatttcctagggc
#ttcgatggacgaaatggagaaaagggtgaaacaggtgctcctggattaaa
#ggtaaatcacaacaaaaatca
#</DNA>
#</SEQUENCE>
#</DASDNA>


for line in $(sed -n '1,$ p' $INPUTFILE)
do
	orig_coordinate=$(echo $line | sed 's/-/:/' | cut -d, -f1 | sed 's/-/,/')
	# chr2:189864471,189864641
	start_coord=$(echo $orig_coordinate | cut -d: -f2 | cut -d, -f1)
	# add PAD_5PRIME to get 5-prime site 
	((start_coord=$start_coord + $PAD_5PRIME))
	end_coord=$(echo $orig_coordinate | cut -d: -f2 | cut -d, -f2)
	# subtract PAD_3PRIME to get 3-prime site
	((end_coord=$end_coord - $PAD_3PRIME))
	chr=$(echo $orig_coordinate | cut -d: -f1) 
	coordinate=${chr}:${start_coord},${end_coord}
	wget -q -o /dev/null -O /tmp/bob http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=$coordinate > /dev/null 2>&1
	sequence=$(awk '/\<DNA/,/\/DNA/ {print}' /tmp/bob | grep -v DNA | tr -d '[:space:]') 
	echo ${line},${sequence} | sed 's/\-/\:/' >> $OUTPUTFILE
	echo -n '#'
done

sed -i '' '1,1 s/$/sequence/' $OUTPUTFILE
