## ROC Curves
make_rocs <- function(){
  require(neuralnet)
  require(randomForest)
  resultsList <- vector('list', 2)
  source("Fig3_and_Stats/PredictionCheck.R")
  load('StringPredictions50.RData')
  set.seed(123L)
  cutoff_positive_case <- 3L
  nfold <- 5L

  posTRUE <- which(FullSubscoresString$ActualCat <= cutoff_positive_case)
  posFALSE <- which(FullSubscoresString$ActualCat > cutoff_positive_case)
  all_props <- EWScores$MinTreeBranch

  allbreaks <- seq(0, 8600, by=100)
  num_bins <- length(allbreaks) - 1L
  count1s <- all_props[posTRUE]
  tempbins <- .bincode(count1s, allbreaks)
  tab1 <- tabulate(tempbins, nbins=num_bins)

  subprops <- all_props[posFALSE]
  num_newbins <- .bincode(subprops, allbreaks)
  sampprobs <- (tab1 / tabulate(num_newbins, nbins=num_bins))[num_newbins]
  ss <- sample(posFALSE, length(posTRUE), prob = sampprobs, replace=FALSE)
  posFALSE <- ss

  EWScores <- EWScores[c(posTRUE,posFALSE),]
  FullSubscoresString <- FullSubscoresString[c(posTRUE, posFALSE),]

  actual <- FullSubscoresString$ActualCat <= cutoff_positive_case
  for(i in seq_len(8)){
    cat('\t', colnames(FullSubscoresString)[i+2], '\n')
    scores <- suppressWarnings(vcheckans(FullSubscoresString[,i+2L], actual))
    resultsList[[1]][[colnames(FullSubscoresString)[i+2L]]] <- scores
  }
  for(i in seq_len(12)){
    cat('\t', colnames(EWScores)[i], '\n')
    scores <- suppressWarnings(vcheckans(EWScores[,i], actual))
    resultsList[[2]][[colnames(EWScores)[i]]] <- scores
  }

  ## RF
  scoresE <- suppressWarnings(vcheckans(rowSums(FullSubscoresString[,14:16]), actual))

  resultsList[[2]]$RandomForest <- scoresE

  names(resultsList) <- c("STRING", 'EvoWeaver')
  resultsList$EWScores <- EWScores
  resultsList$StringScores <- FullSubscoresString

  return(resultsList)
}

compute_string_combined_score <- function(data, prior=0.041){
  data <- data / 1000
  data[data > prior] <- (data[data>prior] - prior) / (1-prior)
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
  cat("Ordering:\n")
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
    cat('\t', colnames(data)[possibles[nextchoice]], '\n')
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
  pTP <- rocdata$StringScores$ActualCat <= 2
  for(i in seq_len(nrow(cormat))){
    for(j in seq_len(ncol(cormat))){
      cormat[i,j] <- mean(cor(rocdata$EWScores[pTP,j],
                              rocdata$StringScores[pTP,i+2L],
                              method = 'spearman'),
                          cor(rocdata$EWScores[!pTP,j],
                              rocdata$StringScores[!pTP,i+2L],
                              method = 'spearman'))
    }
  }

  ## Renaming and reordering for consistency with other plots
  colnames(cormat) <- c("G/L Correlation", "P/A Jaccard", "G/L Distance",
                        "P/A MI", "RP ContextTree", "RP MirrorTree",
                        "Tree Distance", "Gene Distance", "Moran's I",
                        "Transcription MI", "Gene Vector", "Sequence Info")
  cormat <- cormat[,c(1,3,2,4,5:7,8,10,9,12,11)]
  rownames(cormat) <- c("Neighborhood", "Cooccurrence", "Gene Fusion",
                        "Experimental", "Coexpression", "Databases", "Textmining")
  cormat <- cormat[c(3,2,4,1,5,7,6),]

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
             srt=45, cex=0.75, adj=c(1,0.5))
      }
    }
  }
  title(main=substitute(bold("Correlation of Scores ("~rho~")")), line=2, adj=0.5)
  title(xlab="STRING Scores", mgp=c(2.5,1,0))
  title(ylab="EvoWeaver Scores", mgp=c(3,1,0))

  legxs <- 1.025
  legys <- 0.05
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
  rect(xleft=legxs, xright=legxe,ybottom=legys,ytop=legye,
       border='black', lwd=1.5)
}

plot_Fig4 <- function(rocdata, sr){
  of <- "4_FigStringEW.pdf"
  ptsize <- 12
  odir <- './'
  c1 <- '#2B6DA8'
  #c2 <- '#E0A608'
  c3 <- 'white'
  pdf(file=file.path(odir, of),
      width=4.3*2, height=4.3, pointsize = ptsize, onefile=TRUE)
  layout(matrix(1:2, nrow=1))
  sr <- c(0.50, sr, rocdata$EvoWeaver$RandomForest$AUROC)
  names(sr) <- c("  Random Guessing", "  Gene Fusion", "+ Cooccurrence",
                        "+ Experimental", "+ Neighborhood", "+ Coexpression",
                        "+ Textmining", "+ Databases",
                 "   EvoWeaver Random Forest")
  width <- 1
  space <- 0.2
  barplot(sr, names.arg='', yaxt='n', ylim=c(0,1),
          main="Performance Predicting KEGG Pathways",
          space=c(0.1,0.6,rep(0.2,length(sr)-3),0.6),
          col=c('grey', rep('#1E88E5', 7), '#FFC107'))
  title(ylab="AUROC", mgp=c(1.5,1,0))
  axis(side=2, at=seq(0,1,0.1), mgp=c(0,0.5,0), cex.axis=0.75)

  offsets <- seq(from=0,by=space, length.out=length(sr))
  offsets <- offsets + 0.1
  offsets[2:length(offsets)] <- offsets[2:length(offsets)] + 0.4
  offsets[length(offsets)] <- offsets[length(offsets)]+0.4
  labs <- names(sr)

  text(x=seq_along(sr) + offsets-0.35*width, y=0.025,
       srt=90, labels=labs, cex=0.75, adj=c(0,0))
  text(x=seq_along(sr)+offsets-1*width, y=sr+0.01,
       labels=sprintf("%.02f",sr), cex=0.75, adj=c(0,0))

  ## Arrow for x-axis
  arrows(y0=-0.03, x0=offsets[2]+1, x1=length(sr)-1+offsets[length(sr)-1]-0.1,
         length=0.1, angle=20, col='black', lwd=3, xpd=NA)
  text(y=-0.075, x=5.5, labels="More STRING Predictors", xpd=NA)

  ## Arrows inside text/database bars
  for(k in c(7,8)){
    arrows(y0=sr[k-1]+0.01, y1=sr[k]-0.01, x0=k+offsets[k]-0.4*width, length=0.025, col=c3)
    lines(x=c(0.1*width, 0.9*width)+k-1+offsets[k], y=rep(sr[k-1],2), col=c3)
    text(x=k+offsets[k]-0.7*width, y=(sr[k]-sr[k-1])/2+sr[k-1],
         labels=sprintf("%.2f",sr[k]-sr[k-1]),
         cex=0.5, srt=90, col=c3)
  }

  ## Labeling last two bars
  yoff <- 0.03
  xoff <- 5.35
  lines(x=c(4+offsets[4]+0.5, 8+offsets[8]-1.15*width),
        y=c(0.8225, sr[8]+0.025), lwd=2, col=c1)
  text(x=5+offsets[5]-xoff, y=0.8+yoff, labels='STRING Total Score',
       col=c1, cex=0.75, adj=c(0,0.5))

  plot_heatmap(rocdata)

  dev.off(dev.list()['pdf'])
}

StringVsEW <- make_rocs()
actual <- StringVsEW$StringScores$ActualCat <= 3
StringRes <- calc_running_STRINGscore(StringVsEW$StringScores[,3:9], actual)
plot_Fig4(StringVsEW, StringRes)
