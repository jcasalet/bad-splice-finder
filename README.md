BAD SPLICE FINDER
=================

1. install maxentpy in your virtual environment

`$ pip install git+https://github.com/kepbod/maxentpy`

2. download the Soemedi data set

`$ BASH/getSoemediData.sh`

3. download the sequences associated with the loci specified in data.csv

`$ BASH/getSequences.sh`

4. get 5' and 3' splice site scores in neighborhood of mutation 

`$ python PYTHON/MyMaxEntScan.py enriched-data.csv > final-data.csv` 

5. build predictive model using logistic regression and random forest

`$ R R/binary_model_selection.R` 
