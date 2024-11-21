library(DECIPHER)

## Set this to TRUE to regenerate the XXX_PID.RData files
## this may take a while if NUM_SAMPLES is big
REGENERATE_PIDS <- FALSE
NUM_SAMPLES <- 10000

## Directory to source data from
datadir <- file.path(SourceDir, "Data", "SupplementalData", "Homology")
outdir <- file.path(SourceDir, "Data", "SupplementalData", "Homology")

figdir <- file.path(SourceDir, "OutputFigures", "SupplFigures")

curdir <- file.path(SourceDir, "SupplementalFigures", "HomologyTest")

if(REGENERATE_PIDS){
  cat("Generating Complexes pairwise PIDs...\n")
  USE_COMPLEXES <- TRUE
  USE_CORUM <- FALSE
  source(file.path(curdir, "FindHomologyForSequences.R"))
  source(file.path(curdir, "FindHomologyForSequencesWithBlast.R"))

  cat("Generating Modules pairwise PIDs...\n")
  USE_COMPLEXES <- FALSE
  USE_CORUM <- FALSE
  source(file.path(curdir, "FindHomologyForSequences.R"))
  source(file.path(curdir, "FindHomologyForSequencesWithBlast.R"))

  cat("Generating Modules pairwise PIDs for False Positives...\n")
  USE_BLAST <- FALSE
  source(file.path(curdir, "FindHomologyFalsePositives.R"))
  USE_BLAST <- TRUE
  source(file.path(curdir, "FindHomologyFalsePositives.R"))

  #cat("Generating CORUM pairwise PIDs...\n")
  #USE_COMPLEXES <- FALSE
  #USE_CORUM <- TRUE
  #source(file.path(curdir, "FindHomologyForSequences.R"))

  source(file.path(curdir, "ReconcileDatafiles.R"))
}

cat("plotting results...\n")
source(file.path(curdir, "PlotHomologyECDF.R"))
