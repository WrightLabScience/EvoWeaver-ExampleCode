library(SynExtend)
basepath <- './'
localpath <- file.path(basepath, "RunningEvoWeaver")
# This contains two trees, two sets of location data, and a species tree
load(file.path(localpath, "ExampleData.RData"))

algos <- c("GainLoss", "Jaccard", "Hamming", "CorrGL", "MutualInformation", "PAPV",
            "MirrorTree", "ContextTree", "TreeDistance",
            "NVDT", "ResidueMI")
colocAlgos <- c("Coloc", "ColocMoran", "TranscripMI")

ew1 <- EvoWeaver(TwoComplexTrees, MySpeciesTree=ComplexSpecTree)
ew2 <- EvoWeaver(TwoLocs, MySpeciesTree=ComplexSpecTree)

predsNoColoc <- predict(ew1, Method=algos, TreeMethods='JRF', useDNA=FALSE, Processors=NULL)
predsColoc <- predict(ew1, Method=colocAlgos, Processors=NULL)

preds <- c(vapply(predsNoColoc, \(x) x[[2]], numeric(1L)),
           vapply(predsColoc, \(x) x[[2]], numeric(1L)))

names(preds) <- c(algos, colocAlgos)

preds
