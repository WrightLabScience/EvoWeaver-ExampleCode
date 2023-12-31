# EvoWeaver-ExampleCode

This Github repository details how to reproduce analyses and figures from our manuscript, "EvoWeaver: Large-scale prediction of gene functional associations from coevolutionary signals". 

EvoWeaver is a part of the SynExtend package, distributed through the Bioconductor software platform ([release version](https://bioconductor.org/packages/release/bioc/html/SynExtend.html), [development version](https://bioconductor.org/packages/devel/bioc/html/SynExtend.html)). SynExtend v1.15.1, available from Bioconductor release v3.19, was used for this manuscript. The most up-to-date version of source code for SynExtend can be found on its [GitHub repository](https://github.com/npcooley/SynExtend).

Raw sequence data was downloaded from [KEGG](https://www.kegg.jp/). The total amount of data is too large to include on GitHub or Zenodo (>100GB). These sequences were used to construct phylogenetic trees using DECIPHER, which are available on Zenodo (see below). Files `ProteinComplexTrees.RData` and `ModulesEvoWeaverObject.RData` correspond to the trees for the Complexes benchmark and the Modules benchmark, respectively. The Modules benchmark data are in a premade `EvoWeaver` object, whereas the Complexes benchmark data is a list of trees. Species trees for the Complexes and Modules benchmarks are also available on Zenodo. Gene index information for the Modules benchmark is available on Zenodo, and gene index information for the complexes benchmark is available in the `Fig3` folder.

These trees were used to calculate coevolutionary signal using `EvoWeaver`. An example of running these analyses for input data is shown in the `RunningEvoWeaver` folder. These files demonstrate the use of EvoWeaver on a pair of phylogenetic trees with associated metadata. Actual analyses used this framework on a compute cluster, with each node running analyses for a single pair of trees. 

The resulting pairwise scores for all algorithms were then used for downstream tasks. Scores for all pairs for all algorithms on the Complexes benchmark are given in `Fig3/ComplexPredsAllPairs.RData`, and for the Modules benchmark in the `ModulePredsAllPairs.RData` file on Zenodo. These pairwise scores were used to generate ROC and PRC curves, using the `Fig3/GenComplexStatistics.R` and `Fig3/GenBlockStatistics.R` scripts. The results of these scorings for the Complexes and Modules benchmark are provided in `Fig2_FigS1/ComplexStatistics.RData` and `Fig2_FigS1/ModuleStatistics.RData`, respectively. Each of these files contain a list, where each entry is the results for a particular algorithm. 

Pairwise scores are also used in comparison with STRING. Details on how STRING scores were parsed are available in the helper script in the `Fig4` folder.

The pairwise scores were also used for multiclass classification. This classification was done with `Fig3/CVMulticlassClassif.R`. This R script generates the results used in Figure 3, as well as Supplemental Table 2. The data saved from this script is available in `ModulesMulticlassData.RData` on Zenodo. Figure 4 used specific scores from the multiclass data, which are available in the `Fig5/B3GNT5_ST6GAL1_plottingdata.RData`. Generating the figure also requires the `ModulePredsAllPairs.RData` file from Zenodo.

We also investigated misclassifications from EvoWeaver. The file `EvoWeaverMisclassifications.xlsx` is an Excel spreadsheet containing information on all of EvoWeaver's misclassifications, wherein EvoWeaver predicted a pair of genes to be either "Direct Connection" or "Same Module" when they ostensibly belonged to either "Same Global Pathway" or "Unrelated" in KEGG. The first 15 lines are manually annotated with additional information on the pairing. Each module block is of the form `ID_SUBMODULE_BLOCK`. For example, `M00917_1_2` refers to the second block of the first submodule of `M00917`, which is `K00559`. The second block of the second submodule is `K00227`, and would be written `M00917_2_2`. "Actual Category" is a number corresponding to the classification schema of Fig. 3--1 represents "Direct Connection", 2 is "Same Module", 3 is "Same Pathway", 4 is "Same Global Pathway", and 5 is "Unrelated". The first 15 lines also include the gene's common name, its KO groups, the primary host organism, and biosynthetic pathways they are involved in. As there were 1,946 misclassifications with at least 50% confidence in either "Direct Connection" or "Same Module", we did not have time to investigate all of them in detail.

Scripts in each folder contain a `basepath` and `localpath` variable if they read/write data. `basepath` refers to the root folder in this repository, and `localpath` is the folder containing the script. If a file is not found, it is likely found on Zenodo.

Figure 1 lacks data and was produced in Biorender--as such this repository does not contain a script to reproduce it.

All associated data is stored on [Zenodo](https://zenodo.org/doi/10.5281/zenodo.8423024) (DOI: 10.5281/zenodo.10266140).

Algorithms sometimes have different names in internal data than they do in the journal article. This is due to difficulties with special characters, a desire to keep function names concise, and other reason. The below table lists the internal and published names:

| Name in Publication | Name in Datafile | Description |
| :----: | :----: | :----: |
| G/L Correlation | CorrGL  | Correlation of gain/loss events |
| G/L Distance | GainLoss | Distance between gain/loss events |
| P/A Jaccard | Jaccard | Jaccard index of presence/absence profiles |
| P/A MI | PAMI | Bidirectional mutual information of presence/absence profiles |
| RP ContextTree | RPCT | ContextTree using cophenetic distance randomly projection |
| RP MirrorTree | RPMT | MirrorTree using cophenetic distance with random projection |
| Tree Distance | TreeDistance | Normalized Robinson-Foulds distance |
| Gene Distance | Coloc | Number of coding regions separating genes |
| Moran's I | ColocMoran | Moran's I applied to gene distance |
| Orientation MI | TranscripMI | Bidirectional mutual information of gene orientation |
| Gene Vector | NVDT | Correlation of Gene Natural Vectors in amino acid space |
| Sequence Info | SequenceMI | Mutual information of sites in concatenated sequence alignment |

