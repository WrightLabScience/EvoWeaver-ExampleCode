## Directory to source data from
datadir <- file.path(SourceDir, "Data", "Multiclass")

figdir <- file.path(SourceDir, "OutputFigures")

curdir <- file.path(SourceDir, "Figure4")

source(file.path(SourceDir, "Data", "HelperScripts", "PredictionCheck.R"))
source(file.path(SourceDir, "Data", "HelperScripts", "ColorPalettes.R"))
source(file.path(curdir, "PlotStringEW.R"))
