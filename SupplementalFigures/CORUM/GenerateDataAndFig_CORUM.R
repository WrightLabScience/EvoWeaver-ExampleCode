library(randomForest)
library(neuralnet)

## Will take a while
REGENERATE_CORUM_DATA_AND_PREDICTIONS <- FALSE

datadir <- file.path(SourceDir, "Data")
outdir <- file.path(SourceDir, "Data", "SupplementalData", "CORUM")
figdir <- file.path(SourceDir, "OutputFigures", "SupplFigures")
curdir <- file.path(SourceDir, "SupplementalFigures", "CORUM")


source(file.path(SourceDir, "Data", "HelperScripts", "PredictionCheck.R"))

source(file.path(curdir, "CalcEnsembleModelsCORUM.R"))
source(file.path(curdir, "PrepPlotCORUMData.R"))
