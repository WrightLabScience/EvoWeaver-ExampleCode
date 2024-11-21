## Transferability is only a problem with ensemble methods
## so how much is performance impacted by just ignoring the ensemble?

basedir <- getwd()
source(file.path(basedir, "Data", "HelperScripts", "PredictionCheck.R"))
source(file.path(basedir, "Data", "HelperScripts", "ColorPalettes.R"))
scalefac <- 4.3
pdf(file.path(basedir, "OutputFigures", "SupplFigures", "SXX_EnsembleComponentComparison.pdf"),
    width=scalefac*2, height=scalefac*2)
par(mar=c(3,2.75,1,1.25)+0.1, mgp=c(1.5,0.5,0))
layout(matrix(1:4, byrow=TRUE, nrow=2))

cvec <- c(1,5,9,12)
all_cols <- c(unlist(EW_shades)[c(15,cvec+0,16,cvec+1, cvec+2)])
pch_types <- c(19,10,7)
cutoffs <- c(0.10, 0.01)
legend_cex <- 0.8
text_cex <- 0.75
med_choice <- 7L
ctr <- 1L
allx <- c()
ally <- c()
bestother <- c()
medianother <- c()
plotpchs <- c()
plotcols <- c()
key <- c("Sum of Algorithm Scores", "Best Ensemble",
         "Best Component Algorithm", "Median Component Algorithm",
         "Sum of Algorithm Z-Scores",
         "Ensemble Trained on CORUM")

partialROC <- function(tpr, fpr, cutoff){
  require(DescTools)
  suppressWarnings(AUC(fpr[fpr <= cutoff], tpr[fpr<=cutoff])) / cutoff
}
build_res <- function(input, cutoffs){
  v <- input$AUROC
  for(i in seq_along(cutoffs)){
    v <- c(v, partialROC(input$TPR, input$FPR, cutoffs[i]))
  }
  v
}

add_change_line <- function(x, y){
  xpos1 <- 0.5
  xdelt <- 0.05
  pctimp <- mean((y - x) / x)
  rawimp <- mean((y - x))
  xp2 <- xpos1-0.5*rawimp
  yp2 <- xpos1+0.5*rawimp

  s <- ifelse(sign(rawimp)==1L, "+", "")
  abline(a=rawimp, b=1, col='grey50', lty=4)
  if(abs(rawimp) > 0.05)
    arrows(x0=xpos1, x1=xp2,
           y0=xpos1, y1=yp2,
           length=0.05, col='grey50')

  yadj <- ifelse(sign(rawimp)==1L, 0.01, -0.025)
  text(x=xp2-yadj, y=yp2+yadj,
       labels=paste0(s, sprintf("%.02f", rawimp)),
       #" (", s, sprintf("%.01f", pctimp*100), "%)"),
       cex=text_cex, col='grey50', srt=45, adj=c(0.5,0))
}

make_plot_for_pair <- function(x, y, xlab, ylab, main=''){
  plot(c(0,1), c(0,1), type='l', lty=2,
       xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i',
       main=main,
       xlab=xlab,
       ylab=ylab)
  points(x, y, pch=plotpchs, col=plotcols)
  #plotlegend()
  add_change_line(x, y)
}

make_curve_for_pair <- function(vals, xlab="False positive rate",
                                ylab="True positive rate", main='',
                                ensembleIsKEGG=FALSE){
  if(ensembleIsKEGG) key[6] <- "Ensemble Trained on KEGG"
  rearr <- c(2,5,1,3,4)
  if(length(vals) == 6){
    rearr <- c(2,5,1,3,6,4)
  }
  vals <- vals[rearr]
  key <- key[rearr]
  plot(c(0,1), c(0,1), type='l', lty=2,
       xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i',
       main=main,
       xlab=xlab,
       ylab=ylab)
  col_rearr <- c(1:4,6,5)
  for(i in seq_along(vals)){
    lines(x=vals[[i]]$FPR, y=vals[[i]]$TPR, col=all_cols[col_rearr[i]])
  }
  leg <- paste0(key, " (", vapply(vals, \(x) sprintf("%.03f", x$AUROC), ''), ')')
  if(length(vals) == 5){
    leg <- c(leg, '')
    all_cols[6] <- NA
  }
  legend('bottomright',
         legend=leg, col=all_cols[col_rearr], lty=1, bty='n', cex=legend_cex, lwd=2)

  # v <- par(c("fig", "new", "mar"))
  # # left, right, bottom, top
  # lrbt <- c(0.30, 0.47, 0.63, 0.77)
  # par(fig=lrbt, new=TRUE, mar=c(0,0,0,0))
  # plot(c(0,1), c(0,1), type='l', xaxs='i', yaxs='i',
  #      col='black', lty=2, lwd=1, xlab='', ylab='',
  #      ylim=c(0, 1), xlim=c(0, 0.02),
  #      main='', axes=FALSE, frame=TRUE,
  #      oma=c(0,0,0,0))
  # # x axis
  # axis(1, seq(0,0.02,0.005),
  #      labels = c("0.00", '', '0.01', '', '0.02'),
  #      cex.axis=1, tck=-0.05, mgp=c(3,0.5,0), cex.lab=1)
  # axis(2, seq(0,1.0,0.2),
  #      labels= c("0.0", '', '0.4', '', '0.8', ''),
  #      cex.axis=1, tck=-0.05, mgp=c(3,0.5,0), cex.lab=1)
  # for(i in seq_along(vals)){
  #   lines(x=vals[[i]]$FPR, y=vals[[i]]$TPR, col=all_cols[i])
  # }
  # par(fig=v$fig, mar=v$mar, new=FALSE)
}

## CORUM Ensemble models
load(file.path(basedir, "Data", "SupplementalData", "CORUM", "CORUMNuclearPredictions.RData"))
library(randomForest)
library(neuralnet)
set.seed(819L)
cln <- colnames(RawResults)
RawResults[is.na(RawResults)] <- 0
train_data <- cbind(as.data.frame(RawResults), isTP=FinalDataset$isTP)

rf <- randomForest(as.data.frame(RawResults), y=as.factor(FinalDataset$isTP), maxnodes=25, )
lr <- glm(isTP ~ ., data=train_data, family='binomial')
nn <- neuralnet(isTP ~ ., hidden=12, threshold=0.5, stepmax=1e6,
                data=train_data, linear.output=FALSE)

## KEGG Complexes
load(file.path(basedir, "Data", "SupplementalData", "Complexes", "ComplexStatistics.RData"))
RawScores[is.na(RawScores)] <- 0
newres <- vcheckans(rowSums(RawScores[,-13L]), RawScores$isTP)
rankres <- vcheckans(rowSums(apply(RawScores[,-13], 2, \(x) (x-mean(x)) / sd(x))), RawScores$isTP)
oldres <- EnsembleComplexStatistics[which.max(vapply(EnsembleComplexStatistics, \(x) x$AUROC, numeric(1L)))][[1]]
aurocs <- vapply(ComplexStatistics, \(x) x$AUROC, numeric(1L))
otheralg <- ComplexStatistics[which.max(aurocs)][[1]]
medalg <- ComplexStatistics[order(aurocs)][[med_choice]]
RandForestPreds <- predict(rf, RawScores, type='prob')[,'TRUE']
LogRegPreds <- predict(lr, RawScores)
NeuralNetPreds <- predict(nn, RawScores)[,1]
corum_to_kegg <- list(Logit=vcheckans(LogRegPreds, RawScores$isTP),
                      RandomForest=vcheckans(RandForestPreds, RawScores$isTP),
                      NeuralNetwork=vcheckans(NeuralNetPreds, RawScores$isTP))
t_oldres <- corum_to_kegg[which.max(vapply(corum_to_kegg, \(x) x$AUROC, numeric(1L)))][[1]]

make_curve_for_pair(list(newres, oldres, otheralg, medalg, rankres, t_oldres), main="KEGG Complexes")

## KEGG Modules
load(file.path(basedir, "Data", "Modules", "ModuleStatistics.RData"))
RawScores[is.na(RawScores)] <- 0
newres <- vcheckans(rowSums(RawScores[,-13L]), RawScores$isTP)
rankres <- vcheckans(rowSums(apply(RawScores[,-13], 2, \(x) (x-mean(x)) / sd(x))), RawScores$isTP)
oldres <- EnsembleBlockStatistics[which.max(vapply(EnsembleBlockStatistics, \(x) x$AUROC, numeric(1L)))][[1]]
aurocs <- vapply(BlockStatistics, \(x) x$AUROC, numeric(1L))
otheralg <- BlockStatistics[which.max(aurocs)][[1]]
medalg <- BlockStatistics[order(aurocs)][[med_choice]]

RandForestPreds <- predict(rf, RawScores, type='prob')[,'TRUE']
LogRegPreds <- predict(lr, RawScores)
NeuralNetPreds <- predict(nn, RawScores)[,1]
corum_to_kegg <- list(Logit=vcheckans(LogRegPreds, RawScores$isTP),
                      RandomForest=vcheckans(RandForestPreds, RawScores$isTP),
                      NeuralNetwork=vcheckans(NeuralNetPreds, RawScores$isTP))

t_oldres <- corum_to_kegg[which.max(vapply(corum_to_kegg, \(x) x$AUROC, numeric(1L)))][[1]]

make_curve_for_pair(list(newres, oldres, otheralg, medalg, rankres, t_oldres), main="KEGG Modules")

## CORUM data
load(file.path(basedir, "Data", "SupplementalData", "CORUM", "CORUMNuclearPredictions.RData"))
newres <- vcheckans(rowSums(RawResults), FinalDataset$isTP)
rankres <- vcheckans(rowSums(apply(RawResults[,-13], 2, \(x) (x-mean(x)) / sd(x))), FinalDataset$isTP)
oldres <- ens_rocs[which.max(vapply(ens_rocs, \(x) x$AUROC, numeric(1L)))][[1]]
aurocs <- vapply(roc_list, \(x) x$AUROC, numeric(1L))
otheralg <- roc_list[which.max(aurocs)][[1]]
medalg <- roc_list[order(aurocs)][[med_choice]]
t_oldres <- TransferRes[which.max(vapply(TransferRes, \(x) x$AUROC, numeric(1L)))][[1]]
make_curve_for_pair(list(newres, oldres, otheralg, medalg, rankres, t_oldres),
                    main="CORUM Complexes", ensembleIsKEGG=TRUE)

## STRING
load(file.path(basedir, "Data", "Multiclass", "StringComparisonData.RData"))
EWScores[is.na(EWScores)] <- 0
newres <- vcheckans(rowSums(EWScores[,1:12]), EWScores$Category<=3)
rankres <- vcheckans(rowSums(apply(EWScores[,1:12], 2, \(x) (x-mean(x)) / sd(x))), EWScores$Category<=3)
oldres <- vcheckans(rowSums(FullSubscoresString[,14:16]), EWScores$Category<=3)
aurocs <- apply(EWScores[,1:12], 2, \(x)vcheckans(x, EWScores$Category<=3)$AUROC)
otheralg <- vcheckans(EWScores[,which.max(aurocs)], EWScores$Category<=3)
medalg <- vcheckans(EWScores[,order(aurocs)][,med_choice], EWScores$Category<=3)

RandForestPreds <- predict(rf, EWScores[,1:12], type='prob')[,'TRUE']
LogRegPreds <- predict(lr, EWScores[,1:12])
NeuralNetPreds <- predict(nn, EWScores[,1:12])[,1]
corum_to_kegg <- list(Logit=vcheckans(LogRegPreds, EWScores$Category<=3),
                      RandomForest=vcheckans(RandForestPreds, EWScores$Category<=3),
                      NeuralNetwork=vcheckans(NeuralNetPreds, EWScores$Category<=3))
t_oldres <- corum_to_kegg[which.max(vapply(corum_to_kegg, \(x) x$AUROC, numeric(1L)))][[1]]

make_curve_for_pair(list(newres, oldres, otheralg, medalg, rankres, t_oldres), main="STRING Comparison")

# make_plot_for_pair(bestother, allx, "Best Component Method", "Sum of Evidence Streams")
# make_plot_for_pair(bestother, ally, "Best Component Method", "Best Ensemble Method")
# make_plot_for_pair(medianother, allx, "Median Component Method", "Sum of Evidence Streams")
# make_plot_for_pair(medianother, ally,  "Median Component Method", "Best Ensemble Method")
# make_plot_for_pair(ally, allx, "Best Ensemble Method", "Sum of Evidence Streams")

# plot.new()
# second_legends <- paste0("pAUROC (FPR < ", floor(cutoffs*100), "%)")
# key <- c("Sum of Evidence", "Best Ensemble",
#          "Best Component Algorithm", "Median Component Algorithm",
#          "Rank-Sum of Evidence",
#          "Ensemble Transfer Learning")
# legend('center', bty='n', legend=key,
#        col=all_cols[seq_along(key)], lty=1, lwd=2)
# legend('center', bty='n',
#        legend=c(key, "", "AUROC",
#                 second_legends,
#                 "", "Equal Performance", "Mean Difference in Performance"),
#        pt.bg=all_cols[seq_along(key)],
#        cex=legend_cex,
#        pch=c(rep(22,length(key)),NA, pch_types,NA,NA,NA),
#        lty=c(rep(NA, length(key)+length(pch_types)+2L), 2L, 4L),
#        pt.cex=1.5)


dev.off(dev.list()['pdf'])
