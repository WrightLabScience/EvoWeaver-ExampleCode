library(randomForest)
library(pheatmap)
library(igraph)

## Directory to source data from
datadir <- file.path(SourceDir, "Data", "Modules")

## Directory to save generated datafiles to
outdir <- file.path(SourceDir, "Data", "Multiclass")

figdir <- file.path(SourceDir, "OutputFigures")

curdir <- file.path(SourceDir, "Figure3")

message("Note: Fig 3 was created using Biorender.
Only panels a,b are generated in R.")
# Generate Multiclass data for original tree
USE_KEGG <- FALSE
source(file.path(curdir, "CVMulticlassClassif.R"))

# Generate the data for KEGG tree
USE_KEGG <- TRUE
source(file.path(curdir, "CVMulticlassClassif.R"))
