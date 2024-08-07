library(DescTools)

check_ans <- function(p, a, cutoff, vectorIn=F, balance=F){
  if (!vectorIn){
    triscores <- p[upper.tri(p)]
    trueMat <- a[upper.tri(a)]

    nonmissing <- !is.na(triscores)
    triscores <- triscores[nonmissing]
    trueMat <- trueMat[nonmissing]
  } else {
    triscores <- p
    trueMat <- a
  }

  if (balance){
    negs <- which(!trueMat)
    pos <- which(trueMat)
    smaller <- min(length(negs), length(pos))

    pos <- c(sample(negs, smaller), sample(pos, smaller))
    triscores <- triscores[pos]
    trueMat <- trueMat[pos]
  }

  cutoff <- quantile(triscores, cutoff)

  checkScores <- triscores > cutoff

  TP <- sum(checkScores & trueMat)
  TN <- sum(!checkScores & !trueMat)
  FP <- sum(checkScores & !trueMat)
  FN <- sum(!checkScores & trueMat)

  sens <- TP / (TP + FN)
  spec <- TN / (FP + TN)
  prec <- TP / (TP + FP)

  if(is.na(sens)) sens <- 1
  if(is.na(spec)) spec <- 1
  if(TP + FP == 0) prec <- NA
  res <- c(sens, spec, prec, TP, TN, FP, FN)
  names(res) <- c('Recall', 'Specificity', 'Precision', 'TP', 'TN', 'FP', 'FN')

  #return(list(Scores=coev_scores, Pos=cogs, Neg=negatives))
  return(res)
}

gen_control <- function(testnums){
  load('/Users/aidan/Nextcloud/RStudioSync/streptomyces_data/PATHWAYS.RData')
  testset <- PathwayToCog[testnums]
  allcogs <- unique(unlist(testset))
  control <- matrix(F, nrow=length(allcogs), ncol=length(allcogs))
  for (j in seq_along(testset)){
    p <- testset[[j]]
    idxs <- sapply(p, function(x) which(x==allcogs))
    allcombos <- combn(idxs, 2)
    for (i in 1:ncol(allcombos)){
      i1 <- allcombos[1,i]
      i2 <- allcombos[2,i]
      control[i1, i2] <- T
      control[i2, i1] <- T
    }
  }
  rownames(control) <- colnames(control) <- allcogs
  return(control)
}

gen_roc <- function(p, a, nsub=100, checkfxn = check_ans, vectorIn=F, balance=F){
  ticks <- nsub + 2
  gradations <- seq(0, 1, length.out=ticks)
  TPR <- rep(0, ticks)
  FPR <- rep(0, ticks)
  Prec <- rep(NA, ticks)
  for ( i in 2:(ticks-1)){
    s <- checkfxn(p, a, gradations[i], vectorIn=vectorIn, balance=balance)
    #Recall, Specificity, Precision, TP, TN, FP, FN
    TPR[i] <- s[1]
    FPR[i] <- 1-s[2]
    Prec[i] <- s[3]
  }
  TPR[1] <- 1
  FPR[1] <- 1
  TPR <- rev(TPR)
  FPR <- rev(FPR)
  area <- round(AUC(FPR, TPR), 3)
  puc <- round(AUC(TPR[!is.na(Prec)], Prec[!is.na(Prec)]), 3)


  return(list(TPR=TPR, FPR=FPR, Prec=Prec, auc=area, puc=puc))
}

vcheckans <- function(preds, control){
  op <- preds
  control <- control[!is.na(preds)]
  preds <- preds[!is.na(preds)]
  o <- order(preds, decreasing=T)
  control <- control[o]
  preds <- preds[o]
  control <- control == 1
  t <- rep(F, length(preds))
  TP <- TN <- FP <- FN <- rep(0, length(preds))

  for ( i in 1:length(preds) ){
    cval <- preds[i]
    t[preds >= cval] <- T
    TP[i] <- sum(t & control)
    TN[i] <- sum(!t & !control)
    FP[i] <- sum(t & !control)
    FN[i] <- sum(!t & control)
  }

  sens <- TP / (TP + FN)
  spec <- TN / (FP + TN)
  prec <- TP / (TP + FP)

  #res <- list(sens, spec, prec, TP, TN, FP, FN)
  #names(res) <- c('Recall', 'Specificity', 'Precision', 'TP', 'TN', 'FP', 'FN')

  TPR <- c(0, sens)
  FPR <- c(0, 1-spec)
  prec <- c(1,prec)
  area <- suppressWarnings(round(AUC(FPR, TPR), 3))
  puc <- suppressWarnings(round(AUC(TPR, prec), 3))

  res <- list(TPR, FPR, prec, area, puc, op)
  names(res) <- c('TPR', 'FPR', "Precision", 'AUROC', 'AUPRC', 'Predictions')
  return(res)
}

find_elbow <- function(statistics_entry){
  x <- statistics_entry$FPR
  y <- statistics_entry$TPR
  firstpoint <- c(x[1], y[1])
  lastpoint <- c(x[length(x)], y[length(y)])
  d <- function(pt){
    v1 <- firstpoint - lastpoint
    v2 <- pt - firstpoint
    m <- cbind(v1, v2)
    abs(det(m)) / sqrt(sum(v1 * v1))
  }

  dists <- vapply(seq_along(x), \(i) d(c(x[i], y[i])), numeric(1L))
  pos <- which.max(dists)
  p <- c(x[pos], y[pos], statistics_entry$Predictions[pos])
  names(p) <- c("FPR", "TPR", "Cutoff")
  return(p)
}
