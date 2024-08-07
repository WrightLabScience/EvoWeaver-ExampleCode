LOWER_BOUND <- 0.65
MAX_PAUC <- 0.02
LEGEND_CEX <- 1.0
TEXT_CEX <- 1.0
NORMALIZE_PAUC <- TRUE

algs <- c("PAJaccard", "GLMI", "PAOverlap", "GLDistance",
          "RPContextTree", "MoransI") #
ensalgs <- c("Logit", "RandomForest", "NeuralNetwork")

oldnames <- c(algs, ensalgs)
newnames <- c("P/A Jaccard", "G/L MI", "P/A Overlap", "G/L Distance",
              "RP ContextTree", "Moran's I",
              "Logistic Regression", "Random Forest", "Neural Network")

load(file.path(datadir, "SupplementalData", "Complexes", "ComplexStatisticsKEGG.RData"))
ComplexKEGGStat <- c(ComplexStatistics[algs], EnsembleComplexStatistics[ensalgs])
names(ComplexKEGGStat) <- newnames
ComplexKEGGStat_low <- vapply(ComplexKEGGStat, \(x) suppressWarnings(AUC(x=x$FPR, y=x$TPR, to=MAX_PAUC)), numeric(1L))
ComplexKEGGStat <- vapply(ComplexKEGGStat, \(x) x$AUROC, numeric(1L))

load(file.path(datadir, "SupplementalData", "Complexes", "ComplexStatistics.RData"))
ComplexOrigStat <- c(ComplexStatistics[algs], EnsembleComplexStatistics[ensalgs])
names(ComplexOrigStat) <- newnames
ComplexOrigStat_low <- vapply(ComplexOrigStat, \(x) suppressWarnings(AUC(x=x$FPR, y=x$TPR, to=MAX_PAUC)), numeric(1L))
ComplexOrigStat <- vapply(ComplexOrigStat, \(x) x$AUROC, numeric(1L))

load(file.path(datadir, "Modules", "ModuleStatistics.RData"))
BlockKEGGStat <- c(BlockStatistics[algs], EnsembleBlockStatistics[ensalgs])
names(BlockKEGGStat) <- newnames
BlockKEGGStat_low <- vapply(BlockKEGGStat, \(x) suppressWarnings(AUC(x=x$FPR, y=x$TPR, to=MAX_PAUC)), numeric(1L))
BlockKEGGStat <- vapply(BlockKEGGStat, \(x) x$AUROC, numeric(1L))

load(file.path(datadir, "Modules", "ModuleStatisticsKEGG.RData"))
BlockOrigStat <- c(BlockStatistics[algs], EnsembleBlockStatistics[ensalgs])
names(BlockOrigStat) <- newnames
BlockOrigStat_low <- vapply(BlockOrigStat, \(x) suppressWarnings(AUC(x=x$FPR, y=x$TPR, to=MAX_PAUC)), numeric(1L))
BlockOrigStat <- vapply(BlockOrigStat, \(x) x$AUROC, numeric(1L))

pdf(file.path(figdir, "SXX_ReferenceTreeComparison.pdf"),
    onefile=TRUE, width=4.3*2, height=4.3*2)
par(mar=c(2.75,3,1.25,1)+0.1, mgp=c(1.75,0.5,0))
layout(matrix(1:4, nrow=2, byrow=TRUE))

alldata <- rbind(ComplexOrigStat, ComplexKEGGStat,
                 BlockOrigStat, BlockKEGGStat,
                 ComplexOrigStat_low, ComplexKEGGStat_low,
                 BlockOrigStat_low, BlockKEGGStat_low)
## rearrange to same ordering as in Figure 2
alldata <- alldata[,c(1,4,3,2,5,6,8,9,7)]
colList <- list(red=c('#E82B70','#C80B50','#A80030','#780000'),
                blue=c('#1E88E5','#0E68C5','#0048A5','#0028C5'),
                yellow=c('#FFC107','#DFA100','#BF8100'),
                other=c('#5CDC63','#3CBC83','#7CFC83'),
                black=c('black', 'gray60', 'gray80'))

cols <- c(colList$red, colList$blue[1], colList$yellow[3], 'black', 'gray40', 'gray60')
pchs <- c(16,15,17,8,16,16,16,15,8)
dist_arr <- 0.15
o_arr <- 0.005
yl <- "KEGG Reference Tree (AUROC)"
xl <- "Original Reference Tree (AUROC)"
titles <- c("Complexes Benchmark", "Modules Benchmark", "", "")
axticks <- seq(0.65,1,0.05)
axlabs <- sprintf("%.02f",axticks)
for(i in 1:4){
  if(i == 3){
    dist_arr <- dist_arr * (MAX_PAUC/(1-LOWER_BOUND))
    o_arr <- o_arr * (MAX_PAUC/(1-LOWER_BOUND))
    LOWER_BOUND <- 0
    yl <- "KEGG Reference Tree (Partial AUROC)"
    xl <- "Original Reference Tree (Partial AUROC)"
    axticks <- seq(0,0.02,0.005)
    if(NORMALIZE_PAUC)
      axlabs <- sprintf("%.02f",axticks/MAX_PAUC)
    else
      axlabs <- sprintf("%.03f",axticks)
  }

  ## complexes, blocks
  plot(c(LOWER_BOUND,1), c(LOWER_BOUND,1), type='l', lty=2,
       ylim=c(LOWER_BOUND,ifelse(i<3,1,MAX_PAUC)),
       xlim=c(LOWER_BOUND,ifelse(i<3,1,MAX_PAUC)),
       ylab=yl,
       xlab=xl,
       xaxs='i', yaxs='i',
       axes=FALSE, frame.plot=TRUE,
       main=titles[i])
  axis(1, at=axticks, labels=c('', axlabs[-1]))
  axis(2, at=axticks, labels=c('', axlabs[-1]))
  p <- axticks[1] - 0.03*(max(axticks)-min(axticks))
  text(p, p, labels = axlabs[1], srt=-45, xpd=NA)
  points(x=alldata[i*2-1,], alldata[i*2,], col=cols, pch=pchs)
  # up arrow
  arrows(x0=LOWER_BOUND+o_arr, y0=LOWER_BOUND+o_arr,
         x1=LOWER_BOUND+o_arr, y1=LOWER_BOUND+dist_arr,
         lwd=2, length=0.05, col='#BF8100')
  text(adj=c(0.5,0.5),
       x=LOWER_BOUND+2.5*o_arr,
       y=LOWER_BOUND+o_arr+0.5*dist_arr,
       labels="KEGG tree better",
       cex=TEXT_CEX, col='#BF8100', srt=90, font=2)
  arrows(x0=LOWER_BOUND+o_arr, y0=LOWER_BOUND+o_arr,
         x1=LOWER_BOUND+dist_arr, y1=LOWER_BOUND+o_arr,
         lwd=2, length=0.05, col='#0028C5')
  text(adj=c(0.5,0.5),
       x=LOWER_BOUND+o_arr+0.55*dist_arr,
       y=LOWER_BOUND+2.5*o_arr,
       labels="Orig. tree better",
       cex=TEXT_CEX, col='#0028C5', font=2)
  legend('bottomright', inset=c(0.015,0.015),
         pch=pchs, col=cols, legend=colnames(alldata), cex=LEGEND_CEX)
}

dev.off(dev.list()['pdf'])
