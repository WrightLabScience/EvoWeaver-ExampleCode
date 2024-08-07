library(igraph)

## Directory to source data from
datadir <- file.path(SourceDir, "Data")

figdir <- file.path(SourceDir, "OutputFigures")

curdir <- file.path(SourceDir, "Figure6")

source(file.path(curdir, "PlotCaseStudy.R"))
source(file.path(curdir, "BuildAllCaseStudies.R"))
