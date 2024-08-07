# Online Information for "EvoWeaver: Large-scale prediction of gene functional associations from coevolutionary signals"

This repository contains scripts to reproduce the figures in "EvoWeaver: Large-scale prediction of gene functional associations from coevolutionary signals", by Aidan H. Lakshman and Erik S. Wright.

This repository has the following structure (shown in the order they're discussed below):
```
--MasterScript.R
--Figure2
--Figure3
--Figure4
--Figure5
--Figure6
--SupplementalFigures
  |--ComplexesWithValidation
  |--CORUM
  |--HomologyTest
  |--Misclasses
  |--ModulesValidation
  |--PPAlgorithmsExample
  |--Runtime
  |--ReferenceTreeRobustness
--OutputFigures
  |--MainFigures
  |--SupplFigures
--Data
  |--HelperScripts
  |--Modules
  |--Multiclass
  |--SupplementalData
     |--Complexes
     |--CORUM
     |--Homology
     |--ModulesValidation
     |--Runtime
  |--OtherData
--OtherDataGeneration
```

Each folder that generates figures (`FigureX`, and all folders in `SupplementalFigures`) contains a file of the form `GenerateDataAndFig_XXX.R` to generate the associated figure. These scripts depend on a variable `SourceDir`, which is set to the base path to this folder.

All figures can be generated at once by setting the value of `SourceDir` in `MasterScript.R`, then calling:

```
R -f MasterScript.R
```

A detailed explanation of each folder follows.

## Main Figure Folders

Each of these folders contains scripts to generate the figures appearing in the main text. These scripts output figures to `OutputFigures/MainFigures`.

**Note: Some figures were generated using Biorender**. As such, the associated scripts will generate the portions of the figure that are generated in R. The following figures were finalized in Biorender:
- Fig. 1: all panels were created in Biorender.
- Fig. 3: panels (c,d) were created in Biorender, plus legends for (a,b) and the overall arrangement of the figure.
- Fig. 5: panel (a) was created in Biorender, plus legends for (b) and the overall arrangement of the figure.

<!----><a name="figure_2"></a>
### Figure2

`GenerateDataAndFig_Fig2.R` generates the following:

**Data:**
- `Data/Modules/ModuleStatistics.RData`: AUROCs for component predictors and ensemble methods on the Modules benchmark.
- `Data/Modules/ModuleStatisticsKEGG.RData`: Same as above, but using the KEGG Taxonomy tree
- `Data/Modules/EnsembleModels.RData`: Trained ensemble models (Random Forest, Neural Network, Logistic Regression) used in each fold of the benchmark.
- `Data/Modules/EnsembleModelsKEGG.RData`: Same as above, but using the KEGG Taxonomy tree
- `Data/Modules/ExtendedModuleEnsembleMethods.RData`: Extra ensemble methods evaluated on the Modules benchmark, used for supplemental figures.
- `Data/Modules/ExtendedModuleEnsembleMethodsKEGG.RData`: Same as above, but using the KEGG Taxonomy tree.

**Figures:**
- `OutputFigures/MainFigures/2_FigModule.pdf`: Figure 2 from the manuscript.

<!----><a name="figure_3"></a>
### Figure3

`GenerateDataAndFig_Fig3.R` generates the following:

**Data:**
- `Data/Modules/ModulePredsAllPairs.RData`: Updates this file with the multiclass predictions for each pair.
- `Data/Multiclass/MulticlassModuleData.RData`: Multiclass data (see below).
- `Data/Multiclass/MulticlassModuleDataKEGG.RData`: Same as above, but using the KEGG Taxonomy tree
- `Data/Multiclass/MisclassesByEnsembleProbability.csv`: All mispredictions from the multiclass benchmark along with component algorithm scores, original KO identifiers for each module block pair, and ensemble prediction probabilities. This file is used to create `Top100Misclasses.xlsx`, which contains manually determined data on separation in KEGG for each misprediction where EvoWeaver predicted "Direct Connection" but the pair was "Same Pathway".

**Figures:**
- `OutputFigures/MainFigures/3a_Heatmap.png`: Panel (a) from Figure 3 in the manuscript.
- `OutputFigures/MainFigures/3b_featimportance.pdf`: Panel (b) from Figure 3 in the manuscript.
- `OutputFigures/SupplFigures/0ConfHeatmap.png`: Supplemental heatmap showing a confusion matrix of EvoWeaver's predictions at a cutoff of 0% confidence.
- `OutputFigures/SupplFigures/0ConfKEGGHeatmap.png`: Same as above, but using the KEGG reference tree.
- `OutputFigures/SupplFigures/50ConfKEGGHeatmap.png`: Same as above, but using the KEGG reference tree at a cutoff of 50% confidence.

`Data/Multiclass/MulticlassModuleData.RData` contains the following objects:
- `testSets`: train/test sets for each fold of training
- `subpreds`: data used for training/testing (subset of `AllPairs` object, which is contained in `Data/Modules/ModulePredsAllPairs.RData`)
- `Pairings`: module block pairs, where row `i` corresponds to the row `i` in `subpreds`
- `allconf`: feature importances for each component algorithm in the random forest models across all folds of testing
- `allpredictions`: 5-column matrix containing ensemble model predictions. Columns 1-5 correspond to "Direct Connection", "Same Module", "Same Pathway", "Same Global Pathway", "Unrelated", respectively. Row `i` corresponds to row `i` of `Pairings` and `subpreds`.
- `HeatmapsList`: List containing confusion matrices at 0% and 50% confidence. This matrix contains raw counts--the heatmaps shown in the manuscript are normalized by dividing each row by its row sum.
- `Fig5Row`: row of `AllPairs`/`subpreds`/`Pairings`/`allpredictions` used in Figure 5.

<!----><a name="figure_4"></a>
### Figure4

`GenerateDataAndFig_Fig4.R` generates the following:

**Figures:**
- `OutputFigures/MainFigures/4_FigStringEW.pdf`: Figure 4 from the manuscript.
- `OutputFigures/SupplFigures/StringComparisonExtended.pdf`: Supplemental figure containing more information on the STRING vs. EvoWeaver benchmark.

<!----><a name="figure_5"></a>
### Figure5

`GenerateDataAndFig_Fig5.R` generates the following:

**Figures:**
- `OutputFigures/MainFigures/5b_FigDiscoveryViolin.pdf`: Figure 5(b) from the manuscript.
- `OutputFigures/MainFigures/5c-f_FigDiscoveryGraphs.pdf`: Figure 5(c-f) from the manuscript.

This plot depends on `Data/OtherData/B3GNT5_ST6GAL1_plottingdata.RData`, which is a manually created datafile containing the following relevant information:
- `M1,M2`: Gene tree, positional data, original module block name, and gene name for each gene group
- `subspec`: Reference tree subset to the lowest clade containing all organisms in which at least one of the genes is present
- `tree1,tree2`: `subspec` tree with branches and leaves colored by presence/absence and gain/loss patterns of each gene
- `colocdiffs`: Differences in gene indices among the gene groups
- `colocdirs`: Relative gene orientation among the gene groups (1=same, 0=different)
- `seqsM1, seqsM2`: Sequences for each gene group (`character`)
- `nvdtCorr60vecs`: Gene vectors for each gene group
- `RPCVecs`: unnormalized RP-reduced cophenetic vectors for each gene group
- `ofile`: `phylo` object containing only the organisms that have both genes on the same chromosome (used for plotting gene organization data)

<!----><a name="figure_6"></a>
### Figure6

`GenerateDataAndFig_Fig6.R` generates the following:

**Figures:**
- `OutputFigures/MainFigures/6_FigCaseStudy.pdf`: Figure 6 from the manuscript.

<!----><a name="figure_suppl"></a>
## SupplementalFigures Folders

Each of these folders generates at least one component of a supplemental figure.

<!----><a name="figure_scomp"></a>
### ComplexesWithValidation

This folder tests the performance of EvoWeaver on the KEGG Complexes benchmark using gene pair holdouts, then retests using gene holdouts and complex holdouts. The script `GenrateDataAndFigures_Complexes.R` generates the following:

**Data:**
- `Data/SupplementalData/Complexes/ComplexStatistics.RData`: AUROCs for component predictors and ensemble methods on the Complexes benchmark.
- `Data/SupplementalData/Complexes/ComplexStatisticsKEGG.RData`: Same as above, but using the KEGG Taxonomy tree
- `Data/SupplementalData/Complexes/ComplexEnsembleModules.RData`: Trained ensemble models (Random Forest, Neural Network, Logistic Regression) used in each fold of the benchmark.
- `Data/SupplementalData/Complexes/ComplexEnsembleModulesKEGG.RData`: Same as above, but using the KEGG Taxonomy tree
- `Data/SupplementalData/Complexes/ExtendedComplexEnsembleMethods.RData`: Extra ensemble methods evaluated on the Complexes benchmark.
- `Data/SupplementalData/Complexes/ExtendedComplexEnsembleMethodsKEGG.RData`: Same as above, but using the KEGG Taxonomy tree.
- `Data/SupplementalData/Complexes/ComplexStatistics_complexgeneholdouts.RData`: AUROCs for component predictors and ensemble methods on the Complexes benchmark using gene group holdouts.
- `Data/SupplementalData/Complexes/ComplexStatistics_fullcomplexholdouts.RData`: AUROCs for component predictors and ensemble methods on the Complexes benchmark using complex holdouts.


**Figures**:
- `OutputFigures/SupplFigures/SXX_FigComplex.pdf`: Results of EvoWeaver on the Complexes benchmark.
- `OutputFigures/SupplFigures/SXX_ComplexesValidation.pdf`: Comparison of additional ensemble methods, gene group holdouts, and complex holdouts on the Complexes benchmark.

<!----><a name="figure_smod"></a>
### ModulesValidation

This folder performs a similar analysis to `ComplexesWithValidation`. The script `GenerateDataAndFig_ModulesValidation.R` generates the following:

**Data:**
- `Data/SupplementalData/ModulesValidation/ModuleStatistics_modulegeneholdouts.RData`: AUROCs for component predictors and ensemble methods on the Modules benchmark using gene group holdouts.
- `Data/SupplementalData/ModulesValidation/ModuleStatistics_fullmoduleholdouts.RData`: AUROCs for component predictors and ensemble methods on the Modules benchmark using module holdouts.
- `Data/SupplementalData/ModulesValidation/MulticlassModuleData_ModuleHoldout.RData`: Data for the Multiclass benchmark using module holdouts. See "Figure3" for an explanation of this `.RData` file.

**Figures:**
- `OutputFigures/SupplFigures/SXX_ModulesVerification.pdf`: Results of extra ensemble models, gene group holdouts, and module holdouts on the Modules benchmark.
- `OutputFigures/SupplFigures/50ConfModuleHoldoutHeatmap.png`: Multiclass confusion matrix (as in Fig. 3(a)) at a cutoff of 50% confidence using module holdouts.
- `OutputFigures/SupplFigures/0ConfModuleHoldoutHeatmap.png`: Multiclass confusion matrix (as in Fig. 3(a)) at a cutoff of 0% confidence using module holdouts.

<!----><a name="figure_scorum"></a>
### CORUM

This folder analyzes EvoWeaver's predictions on the CORUM dataset. It also benchmarks transfer learning by comparing models trained on KEGG and tested on CORUM, and vice versa. The script `GenerateDataAndFig_CORUM.R` generates the following:

**Figures:**
- `OutputFigures/SupplFigures/SXX_CORUMrocs.pdf`: ROCs of EvoWeaver on KEGG, plus transfer learning results of ensemble methods between KEGG and CORUM.

<!----><a name="figure_shom"></a>
### HomologyTest

This folder analyzes the extent to which different gene groups have a high degree of homology. Ideally, sequences in different gene groups will not have high homology to minimize leakage introduced by comparing two purportedly different groups that actually contain identical sequences. The file `GenerateDataAndFig_Homology.R` generates the following:

**Figures:**
- `OutputFigures/SXX_ECDFHomology.pdf`: empirical cumulative distribution functions (ECDFs) for each dataset. Each plot contains two ECDFs, corresponding to the PID of pairs of sequences taken from either the same gene group or different gene groups. Cutoffs at 40% (often used as a lower cutoff for homology) and 20% PID are marked on each plot.

Additionally, if `REGENERATE_PIDS` is set to `TRUE` in `GenerateDataAndFig_Homology.R`, the following will be generated:

**Data:**
- `Data/SupplementalData/Homology/ComplexPairwise_PID.RData`: Pairwise PIDs for sequences in the Complexes benchmark. `PID_positive` contains PIDs for pairs of sequences drawn from the same gene group, whereas `PID_negative` contains PIDs for pairs of sequences drawn from different gene groups.
- `Data/SupplementalData/Homology/ModulePairwise_PID.RData`: Same as above, but for the Modules benchmark.
- `Data/SupplementalData/Homology/CorumPairwise_PID.RData`: Same as above, but for the CORUM benchmark.

Note that this is set to `FALSE` by default, since this calculation is relatively slow.

<!----><a name="figure_smis"></a>
### Misclasses

This folder generates a plot of the top 100 misclassifications wherein a pair that was "Same Pathway" was classified as "Direct Connection" by EvoWeaver's ensemble predictor. This analysis relies on `Data/Multiclass/Top100Misclasses.xlsx`, which was created manually by examining KEGG pathways for each pair. The script `GenerateDataAndFig_Misclasses.R` generates the following:

**Figures:**
- `OutputFigures/SupplFigures/SXX_Misclasses.pdf`: Ensemble confidence vs. degrees of separation in KEGG for top 100 pairs classified as "Direct Connection" despite actually being "Same Pathway" in KEGG.

<!----><a name="figure_sppex"></a>
### PPAlgorithmsExample

This folder generates example scenarios for the Phylogenetic Profiling algorithms, as well as their expected behavior. `GenerateDataAndFig_PP.R` generates the following:

**Figures:**
- `OutputFigures/SupplFigures/SXX_ExampleTreesPP.pdf`: example scenarios, expected behavior, and actual result for the Phylogenetic Profiling algorithms.

<!----><a name="figure_sruntime"></a>
### Runtime

This folder generates runtime data for each algorithm in EvoWeaver. The file `GenerateDataAndFig_Runtime.R` generates the following:

**Figures:**
- `OutputFigures/SupplFigures/SXX_Runtime.pdf`: Time complexity for each algorithm in EvoWeaver, along with example datapoints.

It is possible to regenerate the benchmark data by setting `REGENERATE_RUNTIMES` to `TRUE` in `GenerateDataAndFig_Runtime.R`. Note that this will take a very long time. If set to `TRUE`, the following file is generated:

**Data:**
- `Data/SupplementalData/Runtime/EvoWeaver_RuntimeBenchmark.RData`: Runtime (in seconds) for each algorithm for different input sizes.

<!----><a name="figure_sreftree"></a>
### ReferenceTreeRobustness

This folder empirically measures the impact of an improved reference tree on EvoWeaver's results. The datafiles for this analysis are generated in `Figure2`, `ModulesValidation`, and `ComplexesWithValidation`. The file `GenerateDataAndFig_ReferenceTree.R` generates the following:

**Figures:**
- `OutputFigures/SupplFigures/SXX_ReferenceTreeComparison.pdf`: Empirical quantification of the impact of an improved reference tree on EvoWeaver's predictions.


## OutputFigures

`OutputFigures` contains two folders, `MainFigures` and `SupplFigures`. The former contains the main figures in the manuscript, and the latter contains the figures from the manuscript's supplemental information.

## Data

### HelperScripts

`Data/HelperScripts` contains some helper scripts used in other analyses. The following scripts are included:
- `ColorPalettes.R`: Color palettes used throughout the EvoWeaver manuscript. The initial colors were generated using David Nichols's website, [Coloring for Colorblindness](https://davidmathlogic.com/colorblind/). Subpalettes were created by darkening/lightening the initial colors.
- `Plot2x2Heatmap.R`: Script to generate the layout used in Figure 2 and the Supplemental Complexes benchmark figure.
- `PredictionCheck.R`: Script to calculate statistics for a predictive algorithm. The main function used is `vcheckans`, which takes as input a vector of predictions and a vector of actuals, and then returns a list containing the original predictions, true positive rate at each threshold, false positive rate at each threshold, area under the receiver-operating curve (AUROC), and area under the precision-recall curve (AUPRC).

### Modules

Data relating to the Modules benchmark. See [`Figure2`](#figure_2) and [`ModulesValidation`](#figure_smod).

### Multiclass

Data relating to the Multiclass benchmark. See [`Figure3`](#figure_3) and [`ModulesValidation`](#figure_smod).

### SupplementalData

Data relating to supplemental figures. See ["`SupplementalFigures` Folders"](#figure_suppl) for more information.

### OtherData

Miscellaneous data used for Supplemental Tables or other analyses.

## OtherDataGeneration

Most of the files in `Data` are datafiles used to generate the final figures. Creating these datafiles is relatively intensive and not easily replicable. Scripts used to generate these files are contained in `OtherGeneration`, though most of these were performed using a distributed compute cluster. Additional larger datafiles are stored on Zenodo (DOI: 10.5281/zenodo.13255153).