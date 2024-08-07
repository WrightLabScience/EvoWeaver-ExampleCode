## Main script to regenerate all figures used in the text
SourceDir <- NULL

if(is.null(SourceDir))
  stop("Please set the value of SourceDir prior to running this script.")

## Main text figures
source(file.path(SourceDir, "Figure2", "GenerateDataAndFig_Fig2.R"))
source(file.path(SourceDir, "Figure3", "GenerateDataAndFig_Fig3.R"))
source(file.path(SourceDir, "Figure4", "GenerateDataAndFig_Fig4.R"))
source(file.path(SourceDir, "Figure5", "GenerateDataAndFig_Fig5.R"))
source(file.path(SourceDir, "Figure6", "GenerateDataAndFig_Fig6.R"))

## Supplemental Figures
SupplDir <- file.path(SourceDir, "SupplementalFigures")
source(file.path(SupplDir, "ComplexesWithValidation", "GenerateDataAndFig_Complexes.R"))
source(file.path(SupplDir, "ModulesValidation", "GenerateDataAndFig_ModuleValidation.R"))
source(file.path(SupplDir, "CORUM", "GenerateDataAndFig_CORUM.R"))
source(file.path(SupplDir, "HomologyTest", "GenerateDataAndFig_Homology.R"))
source(file.path(SupplDir, "PPAlgorithmsExample", "GenerateDataAndFig_PP.R"))
source(file.path(SupplDir, "Runtime", "GenerateDataAndFig_Runtime.R"))
source(file.path(SupplDir, "ReferenceTreeRobustness", "GenerateDataAndFig_ReferenceTree.R"))
source(file.path(SupplDir, "Misclasses", "GenerateDataAndFig_Misclasses.R"))
