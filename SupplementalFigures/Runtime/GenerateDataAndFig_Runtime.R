library(DECIPHER)
library(SynExtend)
library(minpack.lm)

## This will take a long, long time
REGENERATE_RUNTIMES <- FALSE

outdir <- file.path(SourceDir, "Data", "SupplementalData", "Runtime")
datadir <- outdir
figdir <- file.path(SourceDir, "OutputFigures", "SupplFigures")
curdir <- file.path(SourceDir, "SupplementalFigures", "Runtime")

if(REGENERATE_RUNTIMES){
  cat("Measuring runtimes [THIS MAY TAKE A WHILE!]\n")
  source(file.path(curdir, "RuntimeServerScript.R"))
}

cat("Plotting runtimes...\n")
source(file.path(curdir, "RuntimeDualScaling.R"))
