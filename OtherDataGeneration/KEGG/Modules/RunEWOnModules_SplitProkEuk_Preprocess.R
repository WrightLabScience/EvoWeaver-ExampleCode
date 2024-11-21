## This script was used to generate AllPairs data for Modules on distributed compute

## We will need the following:
##	- Pairings.RData: to get the job id
##	- EvoWeaver object
##	- inds object
##  - original species tree
##	- KEGG species tree
NUM_PAIRS <- 100L

PKG_NAME <- "SynExtend_1.15.4.tar.gz"
LIB <- "./libraries/"
dir.create(LIB)
Sys.setenv(TMPDIR=getwd())
install.packages(PKG_NAME, lib=LIB)
.libPaths(c(LIB, .libPaths()))
library(SynExtend)

Args <- commandArgs(trailingOnly=TRUE)
subval <- as.integer(Args[1L])
start_row <- subval*NUM_PAIRS + 1L
end_row <- (subval+1L) * NUM_PAIRS

load('Pairings.RData') # 0MB, Pairings
load("ModulesEvoWeaverObject.RData") # 650MB, pw -
load("ModsWithPositions.RData") # 15MB, inds
load('KEGGSpeciesTree.RData') # tree -
load("ModulesSpeciesTree.RData") # SpecTree -

## load taxonomy information
library(xlsx)
KEGG_taxon <- read.xlsx('KEGGOrganismBreakdown.xlsx', sheetIndex = 1)
KEGG_euks <- KEGG_taxon$KEGG.orgcode[KEGG_taxon$Taxa.Level.1=="Eukaryotes"]
KEGG_proks <- KEGG_taxon$KEGG.orgcode[KEGG_taxon$Taxa.Level.1=="Prokaryotes"]

#end_row <- ifelse(end_row > nrow(Pairings), nrow(Pairings), end_row)

## can we just process everything?
start_row <- 1L
end_row <- nrow(Pairings)

## normalize species tree heights
h1 <- attr(SpecTree, 'height')
h2 <- attr(tree, 'height')

SpecTree <- dendrapply(SpecTree, \(x){attr(x, "height") <- attr(x,'height')/h1;x})
tree <- dendrapply(tree, \(x){attr(x, "height") <- attr(x,'height')/h2;x})

## first create subpw objects to make processing faster
Subpairings <- Pairings[start_row:end_row,]
modules_used <- unique(unlist(Subpairings))
coloc_modules <- coloc_modules[modules_used]
subpw <- unclass(pw[modules_used])

pw_euk <- lapply(subpw, \(x) subset(x, KEGG_euks))
pw_prok <- lapply(subpw, \(x) subset(x, KEGG_proks))
rm(subpw)

coloc_proks <- lapply(coloc_modules, \(x){
  subv <- strsplit(x,'_') |> vapply(.subset, character(1L), 1L)
  pos <- which(subv %in% KEGG_proks)
  x[pos]
})
coloc_euks <- lapply(coloc_modules, \(x){
  subv <- strsplit(x,'_') |> vapply(.subset, character(1L), 1L)
  pos <- which(subv %in% KEGG_euks)
  x[pos]
})

SpecTreeProk <- subset(SpecTree, KEGG_proks)
SpecTreeEuk <- subset(SpecTree, KEGG_euks)
treeProk <- subset(tree, KEGG_proks)
treeEuk <- subset(tree, KEGG_euks)

save(pw_prok, coloc_proks, SpecTreeProk, treeProk, file="ProkaryoteEWData.RData")
save(pw_euk, coloc_euks, SpecTreeEuk, treeEuk, file="EukaryoteEWData.RData")
