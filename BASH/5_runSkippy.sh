#!/bin/bash

TEMPFILE=/tmp/melissa75
SKIPPY_HOME=/Users/jcasaletto/Desktop/GRAD_SCHOOL/UCSC/WINTER_2019/ROTATION/SKIPPY/Skippy_v1.3

if [ $# -ne 1 ]
then
	echo "usage: $0 <skippy-input> "
	exit 1
fi

while IFS='' read -r line || [[ -n "$line" ]]; do
	echo $line > $TEMPFILE 
	#${SKIPPY_HOME}/score_Skippy_75.pl -build hg19 -score_ectopic -positive_strand $TEMPFILE
	${SKIPPY_HOME}/score_Skippy_75.pl -build hg19 -positive_strand $TEMPFILE
	rm $TEMPFILE
done < "$1"
