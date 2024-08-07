library(phytools)
library(SynExtend)
library(circlize)
library(dendextend)
library(sm)

## Directory to source data from
datadir <- file.path(SourceDir, "Data")

figdir <- file.path(SourceDir, "OutputFigures")

curdir <- file.path(SourceDir, "Figure5")

palette("default")
source(file.path(curdir, 'PrepPlotGainLoss.R'))
source(file.path(curdir, 'Fig5PrepPlot.R'))
