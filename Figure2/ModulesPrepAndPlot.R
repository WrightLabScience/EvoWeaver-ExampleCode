outfile <- file.path(figdir, "MainFigures", '2_FigModule.pdf')
load(file.path(outdir, 'ModuleStatistics.RData'))

names(BlockStatistics) <- c("GLDistance", "PAOverlap", "GLMI", "PAJaccard",
                       "RPContextTree", "RPMirrorTree", "TreeDistance",
                       "GeneDistance", "MoransI", "OrientationMI",
                       "GeneVector", "SequenceInfo")
names(BlockStatistics) <- c("G/L Distance", "P/A Overlap",
                            "G/L MI", "P/A Jaccard",
                       "RP ContextTree", "RP MirrorTree", "Tree Distance",
                       "Gene Distance", "Moran's I", "Orientation MI",
                       "Gene Vector", "Sequence Info")
names(EnsembleBlockStatistics) <- c("Logistic Regression", "Random Forest", "Neural Network")
colnames(RawScores)[seq_along(BlockStatistics)] <- names(BlockStatistics)
# Reordering all the algorithms by AUROC

rnames <- names(BlockStatistics)

rnames[1:4] <- names(sort(vapply(BlockStatistics[1:4], \(x) x$AUROC, numeric(1L)),decreasing=TRUE))
rnames[5:7] <- names(sort(vapply(BlockStatistics[5:7], \(x) x$AUROC, numeric(1L)),decreasing=TRUE))
rnames[8:10] <- names(sort(vapply(BlockStatistics[8:10], \(x) x$AUROC, numeric(1L)),decreasing=TRUE))
rnames[11:12] <- names(sort(vapply(BlockStatistics[11:12], \(x) x$AUROC, numeric(1L)),decreasing=TRUE))
BlockStatistics <- BlockStatistics[rnames]
EnsembleBlockStatistics <- EnsembleBlockStatistics[names(sort(vapply(EnsembleBlockStatistics, \(x) x$AUROC, numeric(1L)),decreasing=TRUE))]
RawScores <- RawScores[,c(names(BlockStatistics), 'isTP')]

source(file.path(datadir, "HelperScripts", 'Plot2x2Heatmap.R'))
plot_fig_2x2(BlockStatistics, EnsembleBlockStatistics, RawScores, outfile)
