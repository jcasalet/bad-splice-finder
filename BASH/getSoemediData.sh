#!/bin/bash

# get the data in TSV format
wget http://fairbrother.biomed.brown.edu/data/bulk_download.txt -O data.csv

# convert TSV to CSV
sed -i '' -e $'s/\t/,/g'  data.csv 

