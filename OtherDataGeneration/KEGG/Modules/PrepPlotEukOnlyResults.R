basedir <- getwd()
datadir <- file.path(basedir, "Data", "Modules")
source(file.path(basedir, "Data", "HelperScripts", "PredictionCheck.R"))
source(file.path(basedir, "OtherDataGeneration", "KEGG", "Modules", "Plot2x2HeatmapEuk.R"))
load(file.path(datadir, 'EukOnlyModuleStatistics.RData'))
load(file.path(datadir, "ModuleStatistics.RData"))

RawScores[is.na(RawScores)] <- 0

outdir <- file.path(basedir, "OutputFigures", "SupplFigures")
outfile <- file.path(outdir, "SXX_EukOnlyModules.pdf")

BlockStatistics <- Results$Original[c(3,2,4,1,7,6,5,10,11,12,9,8)]
#RawScores <- df
names(BlockStatistics) <- c("GLDistance", "PAOverlap", "GLMI", "PAJaccard",
                            "RPContextTree", "RPMirrorTree", "TreeDistance",
                            "GeneDistance", "MoransI", "OrientationMI",
                            "GeneVector", "SequenceInfo")
names(BlockStatistics) <- c("G/L Distance", "P/A Overlap",
                            "G/L MI", "P/A Jaccard",
                            "RP ContextTree", "RP MirrorTree", "Tree Distance",
                            "Gene Distance", "Moran's I", "Orientation MI",
                            "Gene Vector", "Sequence Info")
colnames(RawScores)[seq_along(BlockStatistics)] <- names(BlockStatistics)
# Reordering all the algorithms by AUROC

rnames <- names(BlockStatistics)

rnames[1:4] <- names(sort(vapply(BlockStatistics[1:4], \(x) x$AUROC, numeric(1L)),decreasing=TRUE))
rnames[5:7] <- names(sort(vapply(BlockStatistics[5:7], \(x) x$AUROC, numeric(1L)),decreasing=TRUE))
rnames[8:10] <- names(sort(vapply(BlockStatistics[8:10], \(x) x$AUROC, numeric(1L)),decreasing=TRUE))
rnames[11:12] <- names(sort(vapply(BlockStatistics[11:12], \(x) x$AUROC, numeric(1L)),decreasing=TRUE))
BlockStatistics <- BlockStatistics[rnames]
RawScores <- RawScores[,c(names(BlockStatistics), 'isTP')]

colors <- rep(1:4, times=c(4,3,3,2))
for(i in seq_along(BlockStatistics)){
  BlockStatistics[[i]]$Color <- colors[i]
}

plot_fig_2x2(BlockStatistics, RawScores, outfile)





