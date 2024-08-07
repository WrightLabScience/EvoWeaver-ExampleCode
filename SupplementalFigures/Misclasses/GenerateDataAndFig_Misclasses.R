library(xlsx)

datadir <- file.path(SourceDir, "Data")
figdir <- file.path(SourceDir, "OutputFigures", "SupplFigures")
curdir <- file.path(SourceDir, "SupplementalFigures", "Misclasses")

source(file.path(curdir, "PlotMisclasses.R"))
