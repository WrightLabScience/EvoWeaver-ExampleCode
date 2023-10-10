# EvoWeaver-ExampleCode

This Github repository details how to reproduce analyses and figures from our manuscript, "EvoWeaver: Large-scale prediction of gene functional associations from coevolutionary signals". 

EvoWeaver is distributed through the Bioconductor software platform ([release version](https://bioconductor.org/packages/release/bioc/html/SynExtend.html), [development version](https://bioconductor.org/packages/devel/bioc/html/SynExtend.html)).

Raw sequence data was downloaded from [KEGG](https://www.kegg.jp/). The total amount of data is too large to include on GitHub or Zenodo (>100GB). These sequences were used to construct phylogenetic trees using DECIPHER, which are available on Zenodo (see below). Files `ProteinComplexTrees.RData` and `ModulesEvoWeaverObject.RData` correspond to the trees for the Complexes benchmark and the Modules benchmark, respectively. The Modules benchmark data are in a premade `EvoWeaver` object, whereas the Complexes benchmark data is a list of trees. Species trees for the Complexes and Modules benchmarks are also available on Zenodo. Gene index information for the Modules benchmark is available on Zenodo, and gene index information for the complexes benchmark is available in the `Fig3` folder.

These trees were used to calculate coevolutionary signal using `EvoWeaver`. An example of running these analyses for input data is shown in the `RunningEvoWeaver` folder. These files demonstrate the use of EvoWeaver on a pair of phylogenetic trees with associated metadata. Actual analyses used this framework on a compute cluster, with each node running analyses for a single pair of trees. 

The resulting pairwise scores for all algorithms were then used for downstream tasks. Scores for all pairs for all algorithms on the Complexes benchmark are given in `Fig3/ComplexPredsAllPairs.RData`, and for the Modules benchmark in the `ModulePredsAllPairs.RData` file on Zenodo. These pairwise scores were used to generate ROC and PRC curves, using the `Fig3/GenComplexStatistics.R` and `Fig3/GenBlockStatistics.R` scripts. The results of these scorings for the Complexes and Modules benchmark are provided in `Fig2_FigS1/ComplexStatistics.RData` and `Fig2_FigS1/ModuleStatistics.RData`, respectively. Each of these files contain a list, where each entry is the results for a particular algorithm. 

The pairwise scores were also used for multiclass classification. This classification was done with `Fig3/CVMulticlassClassif.R`. This R script generates the results used in Figure 3, as well as Supplemental Table 2. The data saved from this script is available in `ModulesMulticlassData.RData` on Zenodo. Figure 4 used specific scores from the multiclass data, which are available in the `Fig4/B3GNT5_ST6GAL1_plottingdata.RData`. Generating the figure also requires the `ModulePredsAllPairs.RData` file from Zenodo.

All associated data is stored on [Zenodo](https://zenodo.org/record/8423025) (DOI: 10.5281/zenodo.8423025).
