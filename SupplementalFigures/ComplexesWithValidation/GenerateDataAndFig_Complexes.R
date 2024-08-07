library(e1071)
library(randomForest)
library(neuralnet)

datadir <- file.path(SourceDir, "Data", "SupplementalData", "Complexes")
outdir <- file.path(SourceDir, "Data", "SupplementalData", "Complexes")

figdir <- file.path(SourceDir, "OutputFigures", "SupplFigures")

curdir <- file.path(SourceDir, "SupplementalFigures", "ComplexesWithValidation")

## source helper scripts
source(file.path(SourceDir, "Data", "HelperScripts", "PredictionCheck.R"))
source(file.path(SourceDir, "Data", "HelperScripts", "Plot2x2Heatmap.R"))

## generate data
cat("Generating main data...\n")
USE_KEGG <- FALSE
source(file.path(curdir, "GenComplexStatistics.R"))
cat("Rerunning with KEGG species tree...\n")
USE_KEGG <- TRUE
source(file.path(curdir, "GenComplexStatistics.R"))
cat("Generating data with gene group holdouts...\n")
source(file.path(curdir, "Binary_Revalidation_ComplexGeneGroup_holdout.R"))
cat("Generating data with complex holdouts...\n")
source(file.path(curdir, "Binary_Revalidation_Complexholdout.R"))

cat("Plotting Complexes benchmark...\n")
source(file.path(curdir, "ComplexesPrepAndPlot.R"))
cat("Plotting validation plots...\n")
source(file.path(curdir, "ComplexesValidationPlot.R"))
