basepath <- "./"

outfile <- '2_FigModule.pdf'
outfile <- file.path(basepath, "Fig2_FigS1", outfile)
load(file.path(basepath, "Fig2_FigS1", 'ModuleStatistics.RData'))
names(BlockStatistics) <- c("CorrGL", "Jaccard", "GainLoss", "PAMI",
                       "RPCT", "RPMT", "TreeDistance",
                       "Coloc", "ColocMoran", "TranscripMI",
                       "NVDT", "SequenceMI")
names(BlockStatistics) <- c("G/L Correlation", "P/A Jaccard", "G/L Distance",
                            "P/A MI",
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

source(file.path(basepath, 'Fig2_FigS1', 'Plot2x2Heatmap.R'))
plot_fig_2x2(BlockStatistics, EnsembleBlockStatistics, RawScores, outfile)
