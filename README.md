# artmin
Mineral identification routine. From chemical analyses, find the most likely mineral phase(s).

## Getting started
### Prerequisistes
Python 3.x (64-bit recommended) with additional libraries 'numpy', 'pandas', 'scipy' (essential) and 'hdbscan', 'sklearn', 'matplotlib', 'seaborn' (for clustering and plotting).
### Installing
Unzip contents in a folder. Essential for running are run.py (or debug.py), which is the actual routine with user defined inputs, artmin.py and artmin_database.py (modules containing methods called in main routine), a pure phase database (e.g. DBartmin_vXX) and a solid solution database (e.g. DBsolidsolutions_vXX). Other files may be kept in another folder, separate from the essential run files.
## Deployment
## Contributing
All coding was done from a geologist/physicist background with very little coding experience. As such many things are likely to be sub-optimal at best, both on a usability/application design and efficiency of numerical calculations.
## Built with
Python 3.6.5 (64-bit)
coded on Sublime Text 3.0 and Atom Editor 1.25.0
## Authors
* **Antoine Ouellet** - *Initial work* - (antoine.ouellet@gmail.com)
## License
Developped internally at IOS Service Géoscientifiques - (www.iosgeo.com/)
## Acknowledgments
* Alexandre Néron for coaching and first version of mineral identification
