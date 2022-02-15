# TNBC-RPPA

This repository contains the MATLAB code for the network and BMRA for the following paper:

A systems biology approach to investigate kinase signal trans-duction networks that are involved in Triple Negative Breast Cancer resistance to Cisplatin  
Nupur Mukherjee, Alacoque Browne, Laura Ivers,  Tapesh Santra, Mattia Cremona, Bryan T Hennessy, Norma O’Donovan, John Crown, Walter Kolch, Dirk Fey**, Alex J. Eustace1**

./rawdata/ contains the raw RPPA data  
./processeddata/ contains the processed data: mean aggregates over replicates, fold-changes, and infferred network changes  

load_RPPA_data: Load the RPPA data into matlab and pre-processes the data.  
analysis_\*.m: As indicated by their filename, the analysis_\*.m files are scripts performing the analysis tasks. The workflow is as follows:

The workflow is:  
load_RPPA_data  
analysis_sen_vs_res  
analysis_cluster_combined  
analysis_bmra  
