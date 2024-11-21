basedir <- getwd()
AllResultFileName <- file.path(basedir, "Data", "SupplementalData", "CORUM", "ExternalAlgorithmResultsCORUM.RData")
CORUMEWResultFileName <- file.path(basedir, "Data","SupplementalData","CORUM","CORUMNuclearPredictions.RData")
figdir <- file.path(basedir, 'OutputFigures','SupplFigures')

source(file.path(basedir, "Data","HelperScripts","ColorPalettes.R"))
load(AllResultFileName)
load(CORUMEWResultFileName)

LEGEND_CEX <- 0.7
CUTOFF <- 0.1
CUTOFF_Y <- 0.6
roc_list <- roc_list[c(3:4,2,1,5:7,10:12,8:9)]
cols <- rep(1:4, times=c(4,3,3,2))
subcols <- c(1:4, 1:3, 1:3, 1:2)
all_cols <- character(length(roc_list))
for(i in seq_along(roc_list)){
  all_cols[i] <- EW_shades[[cols[i]]][subcols[i]]
}

ltys <- c(1,5,3,4)
pdf(file=file.path(figdir, "SXX_CORUM_Comparison.pdf"), width=4.3*2, height=4.3*2)
par(mar=c(3,2.75,1,1.25)+0.1, mgp=c(1.5,0.5,0))
layout(matrix(1:4, nrow=2, byrow=TRUE))

## Plot EvoWeaver ROCs
aurocs <- vapply(roc_list, \(x) x$AUROC, numeric(1L))
#o <- order(aurocs, decreasing=TRUE)
#roc_list <- roc_list[o]
#aurocs <- aurocs[o]
plot(c(0,1), c(0,1), type='l', lty=2, xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i',
     xlab="False positive rate", ylab='True positive rate', main="CORUM Complexes")
for(i in seq_along(roc_list)){
  col <- roc_list[[i]]$Color
  lines(x=c(0,roc_list[[i]]$FPR,1), y=c(0,roc_list[[i]]$TPR,1),
        col=all_cols[i], lty=ltys[subcols[i]])
}
old_cnames <- c("GLMI","GLDistance","PAOverlap","PAJaccard",
                "RPMirrorTree","RPContextTree","TreeDistance",
                "GeneDistance","MoransI","OrientationMI",
                "GeneVector", "SequenceInfo")
new_cnames <- c("G/L MI", "G/L Distance","P/A Conservation", "P/A Jaccard",
                "RP MirrorTree","RP ContextTree", "Tree Distance",
                "Gene Distance", "Moran's I", "Orientation MI",
                "Gene Vector", "Sequence Info")
names(new_cnames) <- old_cnames
bolding <- rep(1, length(roc_list))
bolding[which.max(aurocs)] <- 2
legend('bottomright', title="AUROCs", title.font=2, bty='n', text.font=bolding,
        legend=paste0(new_cnames[names(roc_list)], " (", sprintf("%.03f", aurocs), ")"),
       col=all_cols, lty=ltys[subcols], cex=LEGEND_CEX)

## Plot partial ROCs (low FPR)
calc_partial_roc <- function(tpr, fpr, cutoff){
  require(DescTools)
  p <- fpr <= cutoff
  suppressWarnings(round(AUC(fpr[p], tpr[p]) / cutoff, 3))
}

aurocs <- vapply(roc_list, \(x) calc_partial_roc(x$TPR, x$FPR, CUTOFF), numeric(1L))
## use previous ordering, don't need to reorder here
plot(c(0,1), c(0,1), type='l', lty=2, xlim=c(0,CUTOFF), ylim=c(0,CUTOFF_Y), xaxs='i', yaxs='i',
     xlab="False positive rate", ylab='True positive rate', main="CORUM Complexes (Low FPR)")
for(i in seq_along(roc_list)){
  lines(x=c(0,roc_list[[i]]$FPR,1), y=c(0,roc_list[[i]]$TPR,1),
        col=all_cols[i], lty=ltys[subcols[i]])
}

bolding <- rep(1, length(roc_list))
bolding[which.max(aurocs)] <- 2
legend('topleft', inset=c(0,0.01), title="Partial AUROCs", title.font=2, bty='n', text.font=bolding,
       legend=paste0(new_cnames[names(roc_list)], " (", sprintf("%.03f", aurocs), ")"),
       col=all_cols, lty=ltys[subcols], cex=LEGEND_CEX)

## Plot other algorithms ##

## Main AUROCs
algnames <- c("Hamming", "Jaccard", "MutInf", "PPP", "SVDPhy", "CladeOScopeAll")
newalgnames <- c("Hamming BPP", "Jaccard BPP", "MI BPP", "PrePhyloPro", "SVD-phy", "CladeOScope")
OtherResults <- AllResults[algnames]
OtherResults$EvoWeaverLogit <- list(Result=ens_rocs$Logit)
OtherResults$EvoWeaverTransfer <- list(Result=TransferRes$Logit)
algnames <- c(algnames, "EvoWeaverLogit", "EvoWeaverTransfer")
newalgnames <- c(newalgnames, "EvoWeaver CORUM Ensemble", "EvoWeaver KEGG Ensemble")
names(newalgnames) <- algnames
aurocs <- vapply(OtherResults, \(x) x$Result$AUROC, numeric(1L))
o <- order(aurocs, decreasing=TRUE)
OtherResults <- OtherResults[o]
aurocs <- aurocs[o]
plot(c(0,1), c(0,1), type='l', lty=2, xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i',
     xlab="False positive rate", ylab='True positive rate', main="")
for(i in seq_along(OtherResults)){
  lines(x=c(0,OtherResults[[i]]$Result$FPR,1), y=c(0,OtherResults[[i]]$Result$TPR,1), col=i, lty=ifelse(i>8, 4, 1))
}
bolding <- rep(1, length(OtherResults))
bolding[which.max(aurocs)] <- 2
legend('bottomright', title="AUROCs", title.font=2, bty='n', text.font=bolding,
       legend=paste0(newalgnames[names(OtherResults)], " (", sprintf("%.03f", aurocs), ")"),
       col=seq_len(length(OtherResults)), lty=1, cex=LEGEND_CEX)


## Partial ROCs
aurocs <- vapply(OtherResults, \(x) calc_partial_roc(x$Result$TPR, x$Result$FPR, CUTOFF), numeric(1L))
## use previous ordering, don't need to reorder here
plot(c(0,1), c(0,1), type='l', lty=2, xlim=c(0,CUTOFF), ylim=c(0,CUTOFF_Y), xaxs='i', yaxs='i',
     xlab="False positive rate", ylab='True positive rate', main="")
for(i in seq_along(OtherResults)){
  lines(x=c(0,OtherResults[[i]]$Result$FPR,1), y=c(0,OtherResults[[i]]$Result$TPR,1), col=i, lty=ifelse(i>8, 4, 1))
}

bolding <- rep(1, length(OtherResults))
bolding[which.max(aurocs)] <- 2
legend('topleft', title="Partial AUROCs", title.font=2, bty='n', text.font=bolding,
       legend=paste0(newalgnames[names(OtherResults)], " (", sprintf("%.03f", aurocs), ")"),
       col=seq_len(length(OtherResults)), lty=1, cex=LEGEND_CEX)

dev.off(dev.list()['pdf'])
