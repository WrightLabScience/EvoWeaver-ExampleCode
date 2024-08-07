outfile <- file.path(figdir, 'SXX_FigComplex.pdf')
load(file.path(datadir, 'ComplexStatistics.RData'))

reordervec <- c("PAJaccard", "GLDistance", "PAOverlap",
                "GLMI",
                "RPContextTree", "RPMirrorTree",
                "TreeDistance",
                "GeneDistance", "OrientationMI", "MoransI",
                "SequenceInfo", "GeneVector")
toPlotMain <- ComplexStatistics[reordervec]
RawDataComplex <- RawScores[,c(reordervec, 'isTP')]
names(toPlotMain)
toPlotEnsemble <- EnsembleComplexStatistics
names(toPlotMain) <- c("P/A Jaccard", "G/L Distance", "P/A Overlap",
                       "G/L MI",
                       "RP ContextTree", "RP MirrorTree", "Tree Distance",
                       "Gene Distance", "Orientation MI", "Moran's I",
                       "Sequence Info", "Gene Vector")
names(toPlotEnsemble) <- c("Logistic Regression", "Random Forest", "Neural Network")
colnames(RawDataComplex) <- c(names(toPlotMain), 'isTP')
# Reordering all the algorithms by AUROC

rnames <- names(toPlotMain)
toPlotEnsemble <- toPlotEnsemble[c(2,3,1)]

plot_fig_2x2(toPlotMain, toPlotEnsemble, RawDataComplex, outfile, isComplex=TRUE)
