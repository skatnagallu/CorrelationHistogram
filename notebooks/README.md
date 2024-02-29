The example notebooks show the analysis of the correlation histograms of the three oxides.

The epos files are not provided. The epos files used in the publication can be downloaded from here [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10677562.svg)](https://doi.org/10.5281/zenodo.10677562) after embargo.

The notebooks can be run after the installation of the correlationHistogramAnalysis module.

All the notebooks have the following structure.
 1. import necessary modules.
 2. Processing DFT energies to identify possible reactions. (The results from the relevant successful DFT calculations were extracted via command line scripting and compiled manually into the final_energy_new.csv .)
 3. Processing APT data (epos files) to overlay identified reactions in the correlation histograms.

 Notes:
 The current version of the code uses ad-hoc readepos routines.

 The future version of the code will use the reading routines provided by IFES Technical Committee. And also expand the use case of the current module.
