BAD SPLICE FINDER
=================

1. install maxentpy in your virtual environment

`pip install git+https://github.com/kepbod/maxentpy`

2. download the Soemedi data set

`BASH/getSoemediData.sh`

3. download the sequences associated with the loci specified in data.csv

`BASH/getSequences.sh`

4. run MyMaxEntScan.py
`python PYTHON/MyMaxEntScan.py data.csv > enriched-data.csv` 
