## Figure 4 generation

This folder contains the necessary materials to generate Figure 4 from the manuscript, comparing STRING to EvoWeaver. 
The `readingSTRINGdata.R` script shows how data from STRING were parsed; the output of this file is `StringPredictions50.RData`.
Note that this script requires text files from STRING not included here--these files are available from STRING or on the Zenodo
repository. 

Use of `StringPredictions50.RData` rather than the output of `readingSTRINGdata.R` is preferred since the former 
contains an additional line of metadata for EvoWeaver scores not recomputed in `readingSTRINGdata.R`. This column is the minimum
number of leaves between each pair of trees and can easily be computed by taking `length(labels(tree))` for each tree in the set, 
then using `min(lab_tree1, lab_tree2)` to get the minimum of the number of leaves between each tree. This computation is straightforward but time consuming, so the results are provided precomputed within `StringPredictions50.RData`. These
results can also be found in `Fig3_and_Stats/ModuleTreeStatistics.RData`.

`PlotStringEW.R` reproduces Figure 4. 
