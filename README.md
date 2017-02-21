# massprf-pipeline
pipeline for genome wide MASSPRF analysis

## Installation

1) Install conda package manager with `pip install conda`

2) Make sure your conda is up to date by running `conda update conda`

3) Update python to 3.5, or (optionally) create a python 3.5 virtual environment
````conda update python````

4) Get package dependencies
    - Add Bioconda channel to Conda:
     
     ````conda config --add channels bioconda````
     
    - Biopython: ````conda install biopython````
    
    - Gffutils: ````conda install gffutils````
    
    - pyvcf: ````conda install pyvcf````
    
# Usage:

To run, type python massprf-pipeline.py <cli>
