# Figures 2 and S1

This folder contains scripts to reproduce graphs from Figures 2 and S1. 
The two `RData` files contain the results of all pairwise analyses for the Complexes and Modules benchmarks. 
`Plot2x2Heatmap.R` is the main workhorse file; this creates the layout and plots all the lines.
The other two R files load the relevant datafiles, format them to a consistent way, and then call `Plot2x2Heatmap.R` to make the plots.
