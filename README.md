# Transcriptomics

MATLAB code for analysis of single cell RNA sequencing data.

The main wrapper class for transcriptomic data is GeneSet. It contains a lot of visualization routines, plus some old clustering methods (Harris et al bioRxiv 2015).

To cluster using the ProMMT algorithm (Harris et al bioRxiv 2017) use the MixNB class.

To use the nbtSNE algorithm use the ComputetSNE function in GeneSet

To run latent factor analysis use NBpca

NOTE: it is quite possible some of this code calls helper functions that are not yet in the github repository. Email me if you get an error, I will upload them.
