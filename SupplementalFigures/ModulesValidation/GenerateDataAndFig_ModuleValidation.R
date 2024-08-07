library(e1071)
library(randomForest)
library(neuralnet)

PLOT_FEAT_IMPORTANCE_MULTICLASS <- FALSE

datadir <- file.path(SourceDir, "Data")
outdir <- file.path(SourceDir, "Data", "SupplementalData", "ModulesValidation")
figdir <- file.path(SourceDir, "OutputFigures", "SupplFigures")
curdir <- file.path(SourceDir, "SupplementalFigures", "ModulesValidation")

## source helper scripts
source(file.path(SourceDir, "Data", "HelperScripts", "PredictionCheck.R"))

## generate data
cat("Generating data with gene group holdouts...\n")
source(file.path(curdir, "Binary_Revalidation_ModuleGeneGroup_holdout.R"))
cat("Generating data with module holdouts...\n")
source(file.path(curdir, "Binary_Revalidation_ModuleHoldout.R"))

cat("Plotting binary validation plots...\n")
source(file.path(curdir, "ModuleValidationPlot.R"))

cat("Generating and plotting multiclass data...\n")
source(file.path(curdir, "CVMulticlassClassif_ModuleHoldouts.R"))
