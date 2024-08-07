## Directory to source data from
datadir <- file.path(SourceDir, "Data")

figdir <- file.path(SourceDir, "OutputFigures")
curdir <- file.path(SourceDir, "SupplementalFigures", "B3GNT5Connections")

source(file.path(SourceDir, "Figure6", "PlotCaseStudy.R"))
source(file.path(curdir, "BuildExtendedFig5Plot.R"))
