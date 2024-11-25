load(file.path(datadir, "StringComparisonData.RData"))
plot1_name <- file.path(figdir, "MainFigures", "4_FigStringEW.pdf")
plot2_name <- file.path(figdir, "SupplFigures", "StringComparisonExtended.pdf")

## ROC Curves
make_rocs <- function(remake_ensemble=FALSE){
  set.seed(768L)
  resultsList <- vector('list', 2)

  cutoff_positive_case <- 3L

  posTRUE <- which(FullSubscoresString$ActualCat <= cutoff_positive_case)
  posFALSE <- which(FullSubscoresString$ActualCat > cutoff_positive_case)

  ss <- sample(posFALSE, length(posTRUE), replace=FALSE)
  posFALSE <- ss

  EWScores <- EWScores[c(posTRUE,posFALSE),]
  FullSubscoresString <- FullSubscoresString[c(posTRUE, posFALSE),]

  actual <- FullSubscoresString$ActualCat <= cutoff_positive_case
  for(i in seq_len(8)){
    #cat('\t', colnames(FullSubscoresString)[i+2], '\n')
    scores <- suppressWarnings(vcheckans(FullSubscoresString[,i+2L], actual))
    resultsList[[1]][[colnames(FullSubscoresString)[i+2L]]] <- scores
  }
  for(i in seq_len(12)){
    #cat('\t', colnames(EWScores)[i], '\n')
    scores <- suppressWarnings(vcheckans(EWScores[,i], actual))
    resultsList[[2]][[colnames(EWScores)[i]]] <- scores
  }
  scoresE <- suppressWarnings(vcheckans(rowSums(FullSubscoresString[,14:16]), actual))
  cat("EvoWeaver Random Forest: ", scoresE$AUROC, '\n')
  resultsList[[2]]$RandomForest <- scoresE

  names(resultsList) <- c("STRING", 'EvoWeaver')
  resultsList$EWScores <- EWScores
  resultsList$StringScores <- FullSubscoresString

  return(resultsList)
}

compute_string_combined_score <- function(data, prior=0.041){
  ## Scores are given as s*1000, need to remove this multiplier
  data <- data / 1000
  ## anything less than the prior gets set to zero
  data[data < prior] <- prior
  data <- (data - prior) / (1-prior)
  outscores <- vapply(seq_len(nrow(data)), \(i) Reduce(`*`, 1-data[i,]), numeric(1L))
  outscores <- (1-outscores) * (1-prior) + prior

  # String says they round, but they use Python's int()
  # which instead truncates towards zero
  return(trunc(outscores*1000))
}

calc_running_STRINGscore <- function(data, actual, usemin=TRUE, prior=0.041){
  possibles <- seq_len(ncol(data))
  totalROCscores <- numeric(ncol(data))
  ordering <- integer(ncol(data))
  cat("STRING Cumulative Ordering:\n")
  for(i in seq_len(ncol(data))){
    if(i != 1){
      curcols <- ordering[seq_len(i-1L)]
    } else {
      curcols <- c()
    }
    nextvals <- vapply(possibles, \(x){
      os <- compute_string_combined_score(data[,c(curcols,x), drop=FALSE])
      vcheckans(os, actual)$AUROC
    }, numeric(1L))
    if(usemin)
      nextchoice <- which.min(nextvals)
    else
      nextchoice <- which.max(nextvals)
    cat('\t', colnames(data)[possibles[nextchoice]], '(', nextvals[nextchoice], ')\n')
    ordering[i] <- possibles[nextchoice]
    possibles <- possibles[-nextchoice]
    totalROCscores[i] <- nextvals[nextchoice]
  }
  names(totalROCscores) <- colnames(data)[ordering]
  return(totalROCscores)
}

plot_heatmap <- function(rocdata){
  # Heatmap plotting
  par(xpd=NA)
  ewd <- rocdata$EvoWeaver
  strd <- rocdata$STRING
  cormat <- matrix(NA, nrow=length(strd)-1, ncol=length(ewd)-1)
  pTP <- rocdata$StringScores$ActualCat <= 3
  rocdata$EWScores[is.na(rocdata$EWScores)] <- 0L
  for(i in seq_len(nrow(cormat))){
    for(j in seq_len(ncol(cormat))){
      cor1 <- cor(rocdata$EWScores[pTP,j],
                  rocdata$StringScores[pTP,i+2L],
                  method = 'spearman')
      cor2 <- cor(rocdata$EWScores[!pTP,j],
                  rocdata$StringScores[!pTP,i+2L],
                  method = 'spearman')
      cormat[i,j] <- mean(ifelse(is.na(cor1), 0, cor1),
                          ifelse(is.na(cor2), 0, cor2))
    }
  }

  #rownames(cormat) <- names(strd)[seq_len(nrow(cormat))]
  colnames(cormat) <- c("P/A Overlap", "P/A Jaccard", "G/L Distance",
                        "G/L MI", "RP ContextTree", "RP MirrorTree",
                        "Tree Distance", "Gene Distance", "Moran's I",
                        "Orientation MI", "Gene Vector", "Sequence Info")
  cormat <- cormat[,c(2,3,1,4,5:7,8,10,9,12,11)]
  rownames(cormat) <- c("Gene Neighborhood", "Cooccurrence", "Gene Fusion",
                        "Experimental Data", "Coexpression Data",
                        "Databases", "Text Mining")
  cormat <- cormat[c(3,4,2,1,5,6,7),]
  #colnames(cormat) <- names(ewd)[seq_len(ncol(cormat))]
  htmcols <- c("#1c6376", '#008585', "#4E998A","#7CAC94","#A4BEA5","#C3D4B1",'white',
               "#eee2a9","#eecd93","#edb588","#ea9c89","#DB6577","#a42f54")

  cols <- colorRampPalette(htmcols)(200)
  plot(NULL, xlab='', ylab='', xlim=c(0,1), ylim=c(0,1),
       type='n', frame.plot=FALSE, xaxt='n',yaxt='n')
  cormat <- t(cormat)
  cormat[seq_len(nrow(cormat)),] <- cormat[rev(seq_len(nrow(cormat))),]
  rownames(cormat) <- rev(rownames(cormat))
  xstart <- 0.1
  xend <- 1
  ystart <- 0.05
  boxwidth <- (xend-xstart) / ncol(cormat)
  boxheight <- boxwidth*0.7
  for(i in seq_len(nrow(cormat))){
    text(x=xstart-0.15*boxwidth,
         y=ystart+(i-0.5)*boxheight,
         labels=rownames(cormat)[i],
         srt=0, cex=0.75, adj=c(1,0.5))
    for(j in seq_len(ncol(cormat))){
      s <- cormat[i,j]
      if(s < 0 && s > -0.05)
        s <- 0
      flab <- sprintf("%.02f",s)
      if(nchar(flab) == 4)
        flab <- paste0(' ', flab)
      rect(xleft=xstart+(j-1)*boxwidth,
           xright=xstart+j*boxwidth,
           ytop=ystart+(i-1)*boxheight,
           ybottom=ystart+i*boxheight,
           col=cols[(round(s, 2)+1)*100],
           border='grey60')
      text(x=xstart+(j-0.5)*boxwidth,
           y=ystart+(i-0.5)*boxheight,
           labels=flab,
           cex=0.75)

      if(i==1){
        text(x=xstart+(j-0.5)*boxwidth,
             y=ystart-0.3*boxheight,
             labels=colnames(cormat)[j],
             srt=30, cex=0.75, adj=c(1,0.5))
      }
    }
  }

  # draw a box around the Cooccurrence:PP methods
  # PP is columns 1:4
  xpos_cooc <- which(colnames(cormat) == "Cooccurrence")
  ys_cooc <- which(rownames(cormat) == "P/A Jaccard") + 1L
  ye_cooc <- which(rownames(cormat) == "G/L MI")
  rect(xleft=xstart+(xpos_cooc-1)*boxwidth, xright=xstart+xpos_cooc*boxwidth,
       ytop=ystart+(ys_cooc-1)*boxheight, ybottom=ystart+(ye_cooc-1)*boxheight,
       border='#D81B60', col=NULL, lwd=2, lty=1)

  ## Gene distance predictors
  xpos_go <- which(colnames(cormat) == "Gene Neighborhood")
  ys_go <- which(rownames(cormat) == "Gene Distance") + 1L
  ye_go <- which(rownames(cormat) == "Moran's I")
  rect(xleft=xstart+(xpos_go-1)*boxwidth, xright=xstart+xpos_go*boxwidth,
       ytop=ystart+(ys_go-1)*boxheight, ybottom=ystart+(ye_go-1)*boxheight,
       border='#45A649', col=NULL, lwd=2, lty=1)

  #title(main=substitute(bold("Correlation of Scores ("~rho~")")), line=2, adj=0.5)
  title(xlab="STRING Evidence Streams", mgp=c(2.25,1,0), cex.lab=0.75, col.lab="#2B6DA8")
  title(ylab="EvoWeaver Component Predictors", mgp=c(2.85,1,0),
        cex.lab=0.75, col.lab="#45A649", adj=0.85)

  legxs <- 1.065
  legys <- 0.34
  legye <- legys+0.5
  legxe <- legxs+0.05
  legh <- (legye - legys)/length(cols)
  for(i in seq_along(cols)){
    rect(xleft=legxs, xright=legxe,
         ybottom=(i-1)*legh+legys, ytop=i*legh+legys,
         col=cols[i], border = NA)
  }
  lines(c(0, 0.015)+legxe, rep(legys,2), col='black', lwd=1.5)
  lines(c(0, 0.015)+legxe, rep(legye,2), col='black', lwd=1.5)
  lines(c(0, 0.015)+legxe, rep((legys+legye)/2,2), col='black', lwd=1.5)
  text(x=rep(0.05,3)+legxe, y=c(legys, legye, (legys+legye)/2)-6*legh,
       labels=c('-1', '+1', ' 0'), cex=0.75, adj=c(0.5,0))
  text(x=legxs-0.025,y=legys+(legye-legys)/2,
       labels=substitute("Spearman's"~rho),
       cex=0.75, srt=90, adj=c(0.5,0.5))
  rect(xleft=legxs, xright=legxe,ybottom=legys,ytop=legye,
       border='black', lwd=1.5)
}

plot_newFig <- function(rocdata, sr, stringdenovo, transferperf){
  ptsize <- 12
  label_cex <- 0.7

  tc0 <- 'grey40'
  tc1 <- EW_fonts[3]
  #tc2 <- EW_shades$red[3]
  tc2 <- EW_fonts[4]
  #tc3 <- EW_fonts[2]
  tc3 <- EW_shades$yellow[2]
  c1 <- '#2B6DA8'
  c3 <- 'white'
  c4 <- "#824484"
  #tc2 <- '#FFC107'
  #c2 <- '#E0A608'
  pdf(file=plot1_name,
      width=4.3*2, height=4.3, pointsize = ptsize, onefile=TRUE)
  layout(matrix(1:2, nrow=1))
  conversion <- c("Cooccurrence", "Experimental Data", "Gene Fusion",
                  "Gene Neighborhood", "Coexpression Data", "Text Mining",
                  "Databases")
  names(conversion) <- c("cooccurence", "experimental", "fusion", "neighborhood",
                         "coexpression", "textmining", "database")
  names(sr) <- conversion[names(sr)]
  names(sr) <- paste0(c('  ', rep("+ ", length(sr)-1)), names(sr))
  sr <- c(0.50, sr, stringdenovo$AUROC, rocdata$EvoWeaver$RandomForest$AUROC, transferperf$AUROC)
  names(sr)[c(1,length(sr)-2, length(sr)-1, length(sr))] <- c("  Random Guessing",
                                                "   STRING (de novo only)",
                                                "   EvoWeaver Random Forest",
                                                "   EvoWeaver Transfer")
  width <- 1
  space <- 0.2
  spacevec <- c(0.1,0.6,rep(0.2,length(sr)-5),0.6,0.2,0.2)
  plot(NULL, xlim=c(0,width*length(sr)+sum(spacevec)), ylim=c(0,1), axes=FALSE, xlab='', ylab='')
  abline(h=0.5,col='black',lty=1)
  barplot(sr, names.arg='', yaxt='n', ylim=c(0,1),
          #main="Performance Predicting KEGG Pathways",
          space=spacevec,
          col=c(tc0, rep(tc1, 7), c4, tc2, tc3), add=TRUE)
  title(ylab="Area Under the ROC Curve (AUROC)",
        mgp=c(1.5,1,0), cex.lab=0.75)
  axis(side=2, at=seq(0,1,0.1), mgp=c(0,0.5,0), cex.axis=0.75)

  offsets <- seq(from=0,by=space, length.out=length(sr))
  offsets <- offsets + 0.1
  offsets[2:length(offsets)] <- offsets[2:length(offsets)] + 0.4
  offsets[length(offsets)-2] <- offsets[length(offsets)-2]+0.4
  offsets[length(offsets)-1] <- offsets[length(offsets)-1]+0.4
  offsets[length(offsets)] <- offsets[length(offsets)]+0.4
  labs <- names(sr)

  text(x=seq_along(sr) + offsets-0.35*width, y=0.025,
       srt=90, labels=labs, cex=label_cex, adj=c(0,0), col=c3)
  text(x=seq_along(sr)+offsets-1*width, y=sr+0.02,
       labels=sprintf("%.02f",sr), cex=label_cex, adj=c(0,0))
  #lines(x=rep(length(sr)-1+offsets[length(sr)]+0.5*space,2),
  #      y=c(-0.15,1.1), col='#D81B60', xpd=NA, lwd=1.5)

  ## Arrow for x-axis
  yarrpos <- -0.06
  xarrst <- offsets[2]+1.1
  xarrend <- length(sr)-3+offsets[length(sr)-3]-0.1
  lineoff <- 0.0135
  arrows(y0=yarrpos, x0=xarrst, x1=xarrend,
         length=0.1, angle=20, col='black', lwd=3, xpd=NA)
  text(y=yarrpos-0.075, x=5.5, labels="Additional STRING Predictors",
       cex=label_cex, xpd=NA)
  lines(y=c(yarrpos+lineoff, yarrpos-lineoff), x=rep(xarrst-0.1,2),col='black', lwd=3, xpd=NA)
  lines(y=c(yarrpos+lineoff, yarrpos-lineoff), x=rep(xarrend+0.1,2),col='black', lwd=3, xpd=NA)

  ## Arrows inside text/database bars
  for(k in c(7,8)){
    arrows(y0=sr[k-1], y1=sr[k]-0.01, x0=k+offsets[k]-0.5*width,
           length=0.05, col=c3, lwd=2, angle=35)
    lines(x=c(0.345*width, 0.655*width)+k-1+offsets[k], y=rep(sr[k-1],2), col=c3, lwd=2)
    text(x=k+offsets[k]-0.5*width, y=sr[k-1]-0.075,
         labels=sprintf("+%.2f",sr[k]-sr[k-1]),
         cex=label_cex, srt=90, col=c3)
  }

  ## Labeling last two bars
  yoff <- 0.045
  xoff <- 3.95
  lines(x=c(5.25+offsets[8]+0.5, 8+offsets[8]-1.15*width),
        y=c(0.9225, sr[8]+0.0375), lwd=1.5, col=c1)
  text(x=5+offsets[8]-xoff, y=0.9+yoff, labels="STRING's Total Score",
       col=c1, cex=label_cex, adj=c(0,0.5))
  # lines(x=c(7+offsets[7], 9+offsets[9]-1.15*width),
  #       y=c(0.95, sr[9]+0.025), lwd=2, col=c2)
  # text(x=7+offsets[7]-xoff, y=0.95+yoff, labels='EvoWeaver Random Forest',
  #      col=c2, cex=0.75, adj=c(0,0.5))
  text(x=c(-2.5,13.5), y=c(1.15,1.15), labels=c("a", "b"),
       xpd=NA, font=2, cex=0.75, adj=c(0,0))

  plot_heatmap(rocdata)

  dev.off(dev.list()['pdf'])
}

plot_STRINGROCs <- function(rocdata, stringrocs, stringdenovo, actual){
  pdf(plot2_name, onefile=TRUE, width=4.3*2, height=4.3*2)
  # par(mar=c(0, 0, 0, 0) + 1,
  # oma=c(0, 0, 0, 0) + 2,
  # mgp=c(1,1,0))
  par(mgp=c(1.5,0.5,0),
      mar=c(1,0.5,1,0)+2.1,
      xpd=NA)
  ordering <- names(stringrocs)
  layout(matrix(1:4, byrow=TRUE, nrow=2))

  strdata <- rocdata$STRING
  strdata <- strdata[order(vapply(strdata, \(x) x$AUROC, numeric(1L)), decreasing=TRUE)]
  cols <- seq_along(strdata)
  names(cols) <- names(strdata)
  ## STRING Individual ROCs
  plot(c(0,1), c(0,1), type='l', lty=3,
       xaxs='i', yaxs='i',
       main="STRING Individual Predictors",
       xlab="False positive rate", ylab="True positive rate")
  for(i in seq_along(strdata)){
    lines(x=strdata[[i]]$FPR, y=strdata[[i]]$TPR, col=i)
  }
  remap <- c("Gene Neighborhood", "Cooccurrence", "Gene Fusion",
             "Experimental Data", "Coexpression Data", "Databases",
             "Text Mining", "Total Score")
  names(remap) <- c("neighborhood", "cooccurence", "fusion",
                    "experimental", "coexpression", "database",
                    "textmining", "combined_score")
  legnames <- remap[names(strdata)]
  legnames <- paste0(legnames, ' (', vapply(strdata, \(x) formatC(x$AUROC, digits=3,format='f'), character(1L)), ')')
  legnames <- c(legnames, "Random Guessing (0.500)")
  legend("bottomright", col=c(cols, 1), bty='n',
         legend=legnames, lty=c(rep(1,length(strdata)), 3), inset=c(0,0.025), cex=0.75)

  rawscores <- rocdata$StringScores
  STRING_cd <- vector('list', length(ordering))
  for(i in seq_along(ordering)){
    tmp <- compute_string_combined_score(rawscores[,ordering[seq_len(i)],drop=FALSE])
    STRING_cd[[i]] <- vcheckans(tmp, actual)
  }
  STRING_cd <- rev(STRING_cd)

  plot(c(0,1), c(0,1), type='l', lty=3,
       xaxs='i', yaxs='i',
       main="STRING Cumulative Predictions",
       xlab="False positive rate", ylab="True positive rate")
  for(i in seq_along(STRING_cd)){
    lines(x=STRING_cd[[i]]$FPR, y=STRING_cd[[i]]$TPR, col=i+1)
  }

  legnames <- remap[rev(ordering)]
  legnames <- paste0("+ ", legnames,  " (",
                     vapply(STRING_cd, \(x) formatC(x$AUROC, 3, format='f'), character(1L)), ')')
  legnames <- c('',legnames, "   Random Guessing (0.500)")
  legend("bottomright", inset=c(0,0.025),
         legend=legnames, col=c(0,cols[rev(ordering)], 1),
         lty=c(0,rep(1, length(ordering)), 3), cex=0.75, bty='n')

  arrows(x0=0.525, x1=0.525, y0=0.065, y1=0.33, length=0.05, lwd=2)
  text(x=c(0,0.025)+0.49, y=0.13, c("Additional","Predictors"),
       srt=90, cex=0.75,
       adj=c(0,0))

  ## plot EvoWeaver predictors
  colList <- list(red=c('#E82B70','#C80B50','#A80030','#780000'),
                  blue=c('#1E88E5','#0E68C5','#0048A5','#0028C5'),
                  yellow=c('#FFC107','#DFA100','#BF8100'),
                  other=c('#5CDC63','#3CBC83','#7CFC83'),
                  black=c('black', 'gray60', 'gray80'))
  ewdata <- rocdata$EvoWeaver
  oldnames <- c("PAOverlap", "PAJaccard", "GLDistance", "GLMI",
                "RPContextTree", "RPMirrorTree", "TreeDistance",
                "GeneDistance", "MoransI", "OrientationMI",
                "GeneVector", "SequenceInfo",
                "RandomForest")
  newnames <- c("P/A Overlap", "P/A Jaccard", "G/L Distance", "G/L MI",
                "RP ContextTree", "RP MirrorTree", "Tree Distance",
                "Gene Distance", "Moran's I", "Orientation MI",
                "Gene Vector", "Sequence Info",
                "Random Forest")
  names(newnames) <- oldnames
  names(ewdata) <- newnames[names(ewdata)]
  ewdata <- ewdata[c(13,3,4,2,1,6,5,7,8,9,10,11,12)]
  rocs <- vapply(ewdata, \(x) x$AUROC, numeric(1L))
  orig_rocs <- rocs
  #o <- order(rocs, decreasing=TRUE)
  #ewdata <- ewdata[o]
  #rocs <- sort(rocs, decreasing=TRUE)
  cols <- c(colList$black[1], colList$red, colList$blue[1:3], colList$yellow[1:3], colList$other[1:2])
  ltys <- c(1,1,5,3,4,1,5,3,1,5,3,1,5)
  plot(c(0,1), c(0,1), type='l', lty=3,
       xaxs='i', yaxs='i',
       main="EvoWeaver Predictors",
       xlab="False positive rate", ylab="True positive rate")
  textauroc <- sprintf("%.03f", orig_rocs)
  legend("bottomright", inset=c(0,0.025),
         legend=paste0(names(ewdata), ' (', textauroc, ')'),
         col=cols,
         lty=ltys, cex=0.75, bty='n')
  for(i in seq_along(ewdata)){
    lines(x=ewdata[[i]]$FPR, y=ewdata[[i]]$TPR, col=cols[i], lty=ltys[i])
  }

  par(xpd=FALSE)
  plot(c(0,1), c(0,1), type='l', lty=3,
       xaxs='i', yaxs='i',
       main="EvoWeaver vs. STRING (low FPR)",
       xlab="False positive rate", ylab="True positive rate",
       xlim=c(0,0.02))
  lines(x=ewdata[["Random Forest"]]$FPR, y=ewdata[["Random Forest"]]$TPR, col='black')
  lines(x=strdata$combined_score$FPR, y=strdata$combined_score$TPR, col='red')
  legend("topright", inset=c(0,0.025),
         legend=c("EvoWeaver Random Forest", "STRING Combined Score",
                  "Random Guessing"),
         lty=c(1,1,3), col=c("black", "red", "black"), cex=0.75, bty='n')
  dev.off(dev.list()['pdf'])
}

StringVsEW <- make_rocs(FALSE)
actual <- StringVsEW$StringScores$ActualCat <= 3
StringRes <- calc_running_STRINGscore(StringVsEW$StringScores[,3:9], actual)
StringDenovoRes <- compute_string_combined_score(StringVsEW$StringScores[,c("neighborhood", "cooccurence", "fusion")])
StringDenovoRes <- vcheckans(StringDenovoRes, actual)

## Reproducing predictions from RSRobustnessPrepPlot.R
## CORUM Ensemble models
load(file.path(SourceDir, "Data", "SupplementalData", "CORUM", "CORUMNuclearPredictions.RData"))
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
EWScores[is.na(EWScores)] <- 0
set.seed(768L)
resultsList <- vector('list', 2)
cutoff_positive_case <- 3L
posTRUE <- which(FullSubscoresString$ActualCat <= cutoff_positive_case)
posFALSE <- which(FullSubscoresString$ActualCat > cutoff_positive_case)
ss <- sample(posFALSE, length(posTRUE), replace=FALSE)
posFALSE <- ss
EWScores <- EWScores[c(posTRUE,posFALSE),]
FullSubscoresString <- FullSubscoresString[c(posTRUE, posFALSE),]

actual <- FullSubscoresString$ActualCat <= cutoff_positive_case
RandForestPreds <- predict(rf, EWScores[,1:12], type='prob')[,'TRUE']
LogRegPreds <- predict(lr, EWScores[,1:12])
NeuralNetPreds <- predict(nn, EWScores[,1:12])[,1]
corum_to_kegg <- list(Logit=vcheckans(LogRegPreds, actual),
                      RandomForest=vcheckans(RandForestPreds, actual),
                      NeuralNetwork=vcheckans(NeuralNetPreds, actual))
TransferLearning <- corum_to_kegg[which.max(vapply(corum_to_kegg, \(x) x$AUROC, numeric(1L)))][[1]]


plot_newFig(StringVsEW, StringRes, StringDenovoRes, TransferLearning)
plot_STRINGROCs(StringVsEW, StringRes, StringDenovoRes, actual)
