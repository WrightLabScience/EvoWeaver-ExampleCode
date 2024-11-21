load(file.path(datadir, "SupplementalData", "CORUM", "CORUMNuclearPredictions.RData"))
load(file.path(datadir, "Modules", "ModuleStatistics.RData"))
source(file.path(datadir,"HelperScripts","ColorPalettes.R"))
LEGEND_CEX <- 1

pdf(file=file.path(figdir, "SXX_CORUMrocs.pdf"), width=4.3*2, height=4.3*2)
par(mar=c(3,3,1,0.5)+0.1, mgp=c(1.5,0.5,0))
layout(matrix(1:4, nrow=2, byrow=TRUE))
## Plot main component ROCs
aurocs <- vapply(roc_list, \(x) x$AUROC, numeric(1L))
#o <- order(aurocs, decreasing=TRUE)
#roc_list <- roc_list[o]
#aurocs <- aurocs[o]
roc_list <- roc_list[c(3,4,1:2,5:7,10:12,8:9)]
cols <- rep(1:4, times=c(4,3,3,2))
subcols <- c(1:4, 1:3, 1:3, 1:2)
all_cols <- character(length(roc_list))
for(i in seq_along(roc_list)){
  all_cols[i] <- EW_shades[[cols[i]]][subcols[i]]
}
plot(c(0,1), c(0,1), type='l', lty=2, xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i',
     xlab="False positive rate", ylab='True positive rate', main="CORUM Nuclear Transport Complexes")
for(i in seq_along(roc_list)){
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
# legend('bottomright', title="Component Algorithms", title.font=2,
#        legend=paste0(new_cnames[names(roc_list)], " (", round(aurocs, 3), ")"), col=seq_len(length(roc_list)), lty=rep(c(1,4), times=c(8,length(roc_list)-8)), cex=0.5)


plot(c(0,1), c(0,1), type='l', lty=2, xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i',
     xlab="False positive rate", ylab='True positive rate', main="CORUM Ensemble Models")
old_names <- c("Logit", "RandomForest", "NeuralNetwork")
new_names <- c("Logistic Regression", "Random Forest", "Neural Network")
names(new_names) <- old_names

## Plot ensembles trained on CORUM
o <- order(vapply(ens_rocs, \(x) x$AUROC, numeric(1L)), decreasing=TRUE)
enscols <- c('brown', 'pink', 'violet')
for(i in seq_len(3)){
  lines(x=ens_rocs[[i]]$FPR, y=ens_rocs[[i]]$TPR, col=enscols[i])
}

## Plot ensembles trained on KEGG
o <- order(vapply(TransferRes, \(x) x$AUROC, numeric(1L)), decreasing=TRUE)
for(j in seq_len(3)){
  v <- TransferRes[[j]]
  lines(v$FPR, v$TPR, col=enscols[j], lty=5)
}

## Train on CORUM, test on KEGG
## RawScores is the KEGG data, RawResults is the CORUM data
set.seed(819L)
cln <- colnames(RawScores)
cln[cln=="SequenceInfo"] <- "SequenceInfo"
colnames(RawScores) <- cln
RawScores[is.na(RawScores)] <- 0
train_data <- cbind(as.data.frame(RawResults), isTP=FinalDataset$isTP)
rf <- randomForest(as.data.frame(RawResults), y=as.factor(FinalDataset$isTP), maxnodes=25, )
RandForestPreds <- predict(rf, RawScores, type='prob')[,'TRUE']

lr <- glm(isTP ~ ., data=train_data, family='binomial')
LogRegPreds <- predict(lr, RawScores)

nn <- neuralnet(isTP ~ ., hidden=12, threshold=0.5, stepmax=1e6,
                data=train_data, linear.output=FALSE)
NeuralNetPreds <- predict(nn, RawScores)[,1]

corum_to_kegg <- list(Logit=vcheckans(LogRegPreds, RawScores$isTP),
                      RandomForest=vcheckans(RandForestPreds, RawScores$isTP),
                      NeuralNetwork=vcheckans(NeuralNetPreds, RawScores$isTP))
names(corum_to_kegg) <- new_names[names(corum_to_kegg)]
names(EnsembleBlockStatistics) <- new_names[names(EnsembleBlockStatistics)]
EnsembleBlockStatistics <- EnsembleBlockStatistics[names(corum_to_kegg)]

plot(c(0,1), c(0,1), type='l', lty=2, xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i',
     xlab="False positive rate", ylab='True positive rate', main="KEGG Module Ensemble Models")

for(i in seq_along(EnsembleBlockStatistics)){
  lines(x=EnsembleBlockStatistics[[i]]$FPR,
        y=EnsembleBlockStatistics[[i]]$TPR,
        col=enscols[i], lty=1)
  lines(x=corum_to_kegg[[i]]$FPR,
        y=corum_to_kegg[[i]]$TPR,
        col=enscols[i], lty=5)
}

plot(NULL, xlim=c(0,1), ylim=c(0,1),
     axes=FALSE, frame.plot=FALSE,
     xlab='', ylab='')
COLS <- c(-0.1, 0.55, 0.8)
ROWS <- c(-0.05,0.65)
LINE_WIDTH <- 1.5
auroc_adjoffset <- 0.25
legend('topleft', title="\nComponent Algorithms", title.font=2,
       legend=new_cnames[names(roc_list)],
       col=all_cols, lwd=LINE_WIDTH,
       lty=ltys[subcols], cex=LEGEND_CEX,
       inset=c(COLS[1],ROWS[1]), xpd=NA, bty='n', title.adj=0)
legend('topleft', title="AUROCs\n(CORUM)", title.font=2, title.adj=0.5,
       legend=sprintf("%.03f", aurocs[names(roc_list)]), adj=auroc_adjoffset,
       cex=LEGEND_CEX, xpd=NA, bty='n',
       inset=c(COLS[2],ROWS[1]))
legend('topleft', title="\nEnsemble Models", title.font=2, title.adj=0,
       legend=c(new_names, paste0(new_names, " (Transfer)")),
       lty=rep(c(1,5), each=3), col=rep(enscols, 2),
       lwd=LINE_WIDTH,
       cex=LEGEND_CEX, inset=c(COLS[1],ROWS[2]), xpd=NA, bty='n')

aurocs_corum <- c(vapply(ens_rocs, \(x) x$AUROC, 0),
                  vapply(TransferRes, \(x) x$AUROC, 0))

aurocs_kegg <- c(vapply(EnsembleBlockStatistics, \(x) x$AUROC, 0),
                 vapply(corum_to_kegg, \(x) x$AUROC, 0))
legend('topleft', title="AUROCs\n(CORUM)", title.font=2, title.adj=0.5,
       legend=sprintf("%.03f", aurocs_corum), adj=auroc_adjoffset,
       cex=LEGEND_CEX, inset=c(COLS[2],ROWS[2]), xpd=NA, bty='n')
legend('topleft', title="AUROCs\n(KEGG)", title.font=2, title.adj=0.5,
       legend=sprintf("%.03f", aurocs_kegg), adj=auroc_adjoffset,
       cex=LEGEND_CEX, inset=c(COLS[3],ROWS[2]), xpd=NA, bty='n')


## Lastly, we'll inset the low FPR region on each plot
lrbt <- c(0.29, 0.47, 0.59, 0.77)
lrbt_coords <- list(lrbt,
                    lrbt + c(0.5,0.5,0,0),
                    lrbt + c(0,0,-0.5,-0.5))

## top left inset
CUTOFF <- 0.02
for(PLOT_NUM in seq_along(lrbt_coords)){
  par(fig=lrbt_coords[[PLOT_NUM]], new=TRUE, mar=c(0,0,0,0))
  plot(c(0,1), c(0,1), type='l', xaxs='i', yaxs='i',
       col='black', lty=2, lwd=1, xlab='', ylab='',
       ylim=c(0, 1), xlim=c(0, CUTOFF),
       main='', axes=FALSE, frame=TRUE,
       oma=c(0,0,0,0))
  # x axis
  axis(1, seq(0,0.02,0.005),
       labels = c("0.00", '', '0.01', '', '0.02'),
       cex.axis=1, tck=-0.05, mgp=c(3,0.5,0), cex.lab=1)
  axis(2, seq(0,1.0,0.2),
       labels= c("0.0", '', '0.4', '', '0.8', ''),
       cex.axis=1, tck=-0.05, mgp=c(3,0.5,0), cex.lab=1)
  if(PLOT_NUM == 1){
    # top left
    for(i in seq_along(roc_list)){
      yv <- roc_list[[i]]$TPR
      xv <- roc_list[[i]]$FPR
      lines(x=c(0,xv,1), y=c(0,yv,1), col=i, lty=ifelse(i>8, 4, 1))
    }
  } else if (PLOT_NUM == 2){
    for(j in seq_len(3)){
      v <- TransferRes[[j]]
      lines(x=ens_rocs[[j]]$FPR, y=ens_rocs[[j]]$TPR, col=enscols[j])
      lines(v$FPR, v$TPR, col=enscols[j], lty=5)
      #cat(names(TransferRes)[j], ':', v$AUROC, '\n')
    }
  } else if (PLOT_NUM == 3){
    for(i in seq_along(EnsembleBlockStatistics)){
      lines(x=EnsembleBlockStatistics[[i]]$FPR,
            y=EnsembleBlockStatistics[[i]]$TPR,
            col=enscols[i], lty=1)
      lines(x=corum_to_kegg[[i]]$FPR,
            y=corum_to_kegg[[i]]$TPR,
            col=enscols[i], lty=5)
    }
  }
}

dev.off(dev.list()['pdf'])
