BAD SPLICE FINDER
=================

1. install maxentpy in your virtual environment

`$ pip install git+https://github.com/kepbod/maxentpy`

2. download the Soemedi data set

`$ BASH/1_getSoemediData.sh`

3. download the sequences associated with the loci specified in data.csv

`$ BASH/2_getSequences.sh`

4. get 5' and 3' splice site scores in neighborhood of mutation 

`$ PYTHON/MyMaxEntScan.py data-with-sequences.csv > data-with-maxentscan.csv` 

5. build predictive model using logistic regression and random forest

`$ Rscript R/binary_model_selection.R data-with-maxentscan.csv` 

OR run the script in RStudio to see the graphs.
