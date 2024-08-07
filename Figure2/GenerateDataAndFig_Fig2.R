## Directory to source data from
datadir <- file.path(SourceDir, "Data")

## Directory to save generated datafiles to
## by default, overwrites the existing files
outdir <- file.path(datadir, "Modules")

figdir <- file.path(SourceDir, "OutputFigures")

curdir <- file.path(SourceDir, "Figure2")

## load function to calculate AUROCs
source(file.path(datadir, "HelperScripts", "PredictionCheck.R"))

# Generate the data for original tree
IS_KEGG <- FALSE
source(file.path(curdir, "GenBlockStatistics.R"))

# Generate the data for KEGG tree
IS_KEGG <- TRUE
source(file.path(curdir, "GenBlockStatistics.R"))

# Plot the tree
source(file.path(curdir, "ModulesPrepAndPlot.R"))
