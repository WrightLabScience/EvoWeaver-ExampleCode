# Figures 2 and S1

This folder contains scripts to reproduce graphs from Figures 2 and S1. 
The two `RData` files contain the results of all pairwise analyses for the Complexes and Modules benchmarks. 
`Plot2x2Heatmap.R` is the main workhorse file; this creates the layout and plots all the lines.
The other two R files load the relevant datafiles, format them to a consistent way, and then call `Plot2x2Heatmap.R` to make the plots.

## Note

Internal algorithm names differ from those in the published article; see the [README](https://github.com/WrightLabScience/EvoWeaver-ExampleCode/blob/1afeb5936470833a74cc2104994f54b363cc1103/README.md) for a description of which algorithm names in data correspond to algorithms in published material.
