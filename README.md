# CorrelationHistogram
For in-depth analysis of correlation histograms in atom probe tomography. 

A unique advantage of straight flight path atom probes is the rich information in correlation histograms of the collected data. The APT datasets contain information regarding the number of ion hits detected in each pulsing window, referred to as multiplicity.For correlation histograms, we focus on events with a multiplicity greater than 1. We preprocess the data to arrange multiple hits with a value greater than two, creating all possible sorted pairs of mass-to-charge values. So, a multiplicity value of 3, corresponding to a detection of triplet $(m_1, m_2, m_3)$ in a pulsing window is processed to 3 sorted pairs { $(m_1, m_2)$, $(m_2, m_3)$, $(m_1, m_3)$ } such that the first mass-to-charge value is always smaller than its partner. The resulting mass-to-charge pairs are plotted on a 2D histogram containing two mass-to-charge axes. As discussed by Saxey, molecular dissociation often leaves a track on the correlation histogram, allowing us to investigate the type of reactions occurring in the analysis of each oxide.


# Installation
Install the latest verion using 
```
pip install correlationHistogramAnalysis
```
# Publication
In this  pre-print[1] the module is used to analyse the oxygen loses in iron oxides analysed from APT. Based on possible reactions predicted from ground state density functional theory (DFT).

# Reference
[1] Kim, Se-Ho, Shalini Bhatt, Daniel K. Schreiber, Jörg Neugebauer, Christoph Freysoldt, Baptiste Gault, and Shyam Katnagallu. "Understanding atom probe's analytical performance for iron oxides using correlation histograms and ab initio calculations." arXiv preprint arXiv:2401.12784 (2024).