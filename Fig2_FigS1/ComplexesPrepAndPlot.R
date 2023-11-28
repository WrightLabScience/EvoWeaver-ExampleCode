basepath <- './'

outfile <- 'S1_FigComplex.pdf'
outfile <- file.path(basepath, "Fig2_FigS1", outfile)
load(file.path(basepath, "Fig2_FigS1", 'ComplexStatistics.RData'))
reordervec <- c("CorrGL", "Jaccard", "GainLoss",
                "MutualInformation",
                "ContextTree", "MirrorTree",
                "TreeDistance",
                "Coloc", "ColocMoran", "TranscripMI",
                "NVDT", "SequenceMI")
toPlotMain <- ComplexStatistics[reordervec]
RawDataComplex <- RawScores[,c(reordervec, 'isTP')]
names(toPlotMain)
toPlotEnsemble <- EnsembleComplexStatistics
names(toPlotMain) <- c("G/L Correlation", "P/A Jaccard", "G/L Distance",
                       "P/A MI",
                       "RP ContextTree", "RP MirrorTree", "Tree Distance",
                       "Gene Distance", "Moran's I", "Transcription MI",
                       "Gene Vector", "Sequence Info")
names(toPlotEnsemble) <- c("Logistic Regression", "Random Forest", "Neural Network")
colnames(RawDataComplex) <- c(names(toPlotMain), 'isTP')
# Reordering all the algorithms by AUROC

rnames <- names(toPlotMain)

# rnames[1:4] <- names(sort(vapply(toPlotMain[1:4], \(x) x$AUROC, numeric(1L)),decreasing=TRUE))
# rnames[5:7] <- names(sort(vapply(toPlotMain[5:7], \(x) x$AUROC, numeric(1L)),decreasing=TRUE))
# rnames[8:10] <- names(sort(vapply(toPlotMain[8:10], \(x) x$AUROC, numeric(1L)),decreasing=TRUE))
# rnames[11:12] <- names(sort(vapply(toPlotMain[11:12], \(x) x$AUROC, numeric(1L)),decreasing=TRUE))
# toPlotMain <- toPlotMain[rnames]
# toPlotEnsemble <- toPlotEnsemble[names(sort(vapply(toPlotEnsemble, \(x) x$AUROC, numeric(1L)),decreasing=TRUE))]

# Reordering to match Modules
toPlotMain <- toPlotMain[c(1,3,2,4,5:8,10,9,12,11)]
toPlotEnsemble <- toPlotEnsemble[c(2,3,1)]

source(file.path(basepath, 'Fig2_FigS1', 'Plot2x2Heatmap.R'))
plot_fig_2x2(toPlotMain, toPlotEnsemble, RawDataComplex, outfile, isComplex=TRUE)
