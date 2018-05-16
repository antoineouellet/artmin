# artmin
Mineral identification routine. From chemical analyses, find the most likely mineral phase(s).

## Getting started
### Prerequisistes
Python 3.x (64-bit recommended) with additional libraries 'numpy', 'pandas', 'scipy' (essential) and 'hdbscan', 'sklearn', 'matplotlib', 'seaborn' (for clustering and plotting).
### Installing
Unzip contents in a folder. Essential for running are run.py (or debug.py), which is the actual routine with user defined inputs, artmin.py and artmin_database.py (modules containing methods called in main routine), a pure phase database (e.g. DBartmin_vXX) and a solid solution database (e.g. DBsolidsolutions_vXX). Other files may be kept in another folder, separate from the essential run files.
##Running the tests
