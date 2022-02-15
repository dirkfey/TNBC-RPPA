%% analysis_bmra ... BMRA analysis on all cell lines
%
%  Requires: ./BMRA/bmra_run.m, ./BMRA/fig_graph_extended.m, 

%% Run BMRA analysis (carefull, this might take some time):

% To run uncomment this code:
%cd BMRA
%bmra_run.m
%cd ..

%% Generate heatmaps of cell networks and calculate conensus networks:
cd BMRA
bmra_analysis_heatmaps
cd ..

%% Visualise consensus network as graph:
cd BMRA
fig_graph_extended.m
cd ..
