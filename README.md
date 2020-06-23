# Ibogaine_analysis_2020

This repository includes the custom programs employed in the work entitled  "Gamma band alterations and REM-like traits underpin the acute effect of the atypical psychedelic ibogaine. Gonzalez, et al, 2020c"

Two main functions are included: Cluster_Permutation_Correction_xfreq and TimeSeries2OPentropy, which perform respectively the  and the permutation entropy analysis. All the other analysis employed in the paper is currently available under the standard Matlab toolboxes or the Chronux and EEGLAB toolboxes. 

Cluster_Permutation_Correction_xfreq performs the cluster-based permutation multiple comparisons correction on the power spectrum data between 2 conditions (ibogaine and control), allowing us to quantify the spectral differences between conditions in a data-driven way. In addition to the main function, the program Example_cluster_correction provides an example of how to use the function with simulated data. 

TimeSeries2OPentropy calculates the permutation entropy of a general time-series, in our case the intracraneal electroencephalogram. For details regarding this analysis consult: "Decreased electrocortical temporal complexity distinguishes sleep from wakefulness. Gonzalez, et al, 2019" available at https://www.nature.com/articles/s41598-019-54788-6
