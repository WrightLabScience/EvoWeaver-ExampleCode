basepath <- './'

## ROC Curves
make_rocs <- function(){
  require(neuralnet)
  require(randomForest)
  resultsList <- vector('list', 2)
  source(file.path(basepath, "Fig3_and_Stats", "PredictionCheck.R"))
  load(file.path(basepath, 'Fig4', 'StringPredictions50.RData'))
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
  cat("\t Ensemble")
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
  colnames(cormat) <- c("G/L Correlation", "P/A Jaccard", "G/L Distance",
                        "P/A MI", "RP ContextTree", "RP MirrorTree",
                        "Tree Distance", "Gene Distance", "Moran's I",
                        "Orientation MI", "Gene Vector", "Sequence Info")
  cormat <- cormat[,c(1,3,2,4,5:7,8,10,9,12,11)]
  rownames(cormat) <- c("Gene Neighborhood", "Cooccurrence", "Gene Fusion",
                        "Experimental Data", "Coexpression Data",
                        "Databases", "Text Mining")
  cormat <- cormat[c(3,2,4,1,5,7,6),]
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
  #title(main=substitute(bold("Correlation of Scores ("~rho~")")), line=2, adj=0.5)
  title(xlab="STRING Evidence Streams", mgp=c(2.25,1,0), cex.lab=0.75)
  title(ylab="EvoWeaver Component Predictors", mgp=c(3,1,0), cex.lab=0.75)
  
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

plot_newFig <- function(rocdata, sr){
  of <- "4_FigStringEW.pdf"
  ptsize <- 12
  odir <- file.path(basepath, "Fig4")
  EW_shades <- list(red=c('#E82B70','#C80B50','#A80030','#780000'),
                    blue=c('#1E88E5','#0E68C5','#0048A5','#0028C5'),
                    yellow=c('#FFC107','#DFA100','#BF8100'),
                    green=c('#5CDC63','#3CBC83','#7CFC83'),
                    black=c('black', 'gray60', 'gray80'))
  EW_fonts <- c('#D81B60', '#E0A608', '#2B6DA8', '#45A649', 'black')
  tc0 <- 'grey40'
  tc1 <- EW_fonts[3]
  tc2 <- EW_shades$yellow[3]
  c1 <- '#2B6DA8'
  c3 <- 'white'
  #tc2 <- '#FFC107'
  #c2 <- '#E0A608'
  pdf(file=file.path(odir, of),
      width=4.3*2, height=4.3, pointsize = ptsize, onefile=TRUE)
  layout(matrix(1:2, nrow=1))
  sr <- c(0.50, sr, rocdata$EvoWeaver$RandomForest$AUROC)
  names(sr) <- c("  Random Guessing", "  Gene Fusion", "+ Cooccurrence",
                 "+ Experimental Data", "+ Gene Neighborhood",
                 "+ Coexpression Data",
                 "+ Text Mining", "+ Databases",
                 "   EvoWeaver Random Forest")
  width <- 1
  space <- 0.2
  barplot(sr, names.arg='', yaxt='n', ylim=c(0,1),
          #main="Performance Predicting KEGG Pathways",
          space=c(0.1,0.6,rep(0.2,length(sr)-3),0.6),
          col=c(tc0, rep(tc1, 7), tc2))
  title(ylab="Area Under ROC Curve (AUROC)",
        mgp=c(1.5,1,0), cex.lab=0.75)
  axis(side=2, at=seq(0,1,0.1), mgp=c(0,0.5,0), cex.axis=0.75)
  #abline(h=0.5, lty=2, col='black')
  offsets <- seq(from=0,by=space, length.out=length(sr))
  offsets <- offsets + 0.1
  offsets[2:length(offsets)] <- offsets[2:length(offsets)] + 0.4
  offsets[length(offsets)] <- offsets[length(offsets)]+0.4
  labs <- names(sr)
  
  text(x=seq_along(sr) + offsets-0.35*width, y=0.025,
       srt=90, labels=labs, cex=0.75, adj=c(0,0), col=c3)
  text(x=seq_along(sr)+offsets-1*width, y=sr+0.02,
       labels=sprintf("%.02f",sr), cex=0.75, adj=c(0,0))
  #lines(x=rep(length(sr)-1+offsets[length(sr)]+0.5*space,2),
  #      y=c(-0.15,1.1), col='#D81B60', xpd=NA, lwd=1.5)
  
  ## Arrow for x-axis
  yarrpos <- -0.03
  xarrst <- offsets[2]+1.1
  xarrend <- length(sr)-1+offsets[length(sr)-1]-0.1
  lineoff <- 0.0135
  arrows(y0=yarrpos, x0=xarrst, x1=xarrend,
         length=0.1, angle=20, col='black', lwd=3, xpd=NA)
  text(y=-0.075, x=5.5, labels="Additional STRING Predictors",
       cex=0.75, xpd=NA)
  lines(y=c(yarrpos+lineoff, yarrpos-lineoff), x=rep(xarrst-0.1,2),col='black', lwd=3, xpd=NA)
  lines(y=c(yarrpos+lineoff, yarrpos-lineoff), x=rep(xarrend+0.1,2),col='black', lwd=3, xpd=NA)
  
  ## Arrows inside text/database bars
  for(k in c(7,8)){
    arrows(y0=sr[k-1], y1=sr[k]-0.01, x0=k+offsets[k]-0.5*width,
           length=0.05, col=c3, lwd=2, angle=35)
    lines(x=c(0.345*width, 0.655*width)+k-1+offsets[k], y=rep(sr[k-1],2), col=c3, lwd=2)
    text(x=k+offsets[k]-0.5*width, y=sr[k-1]-(sr[k]-sr[k-1])*0.8,
         labels=sprintf("+%.2f",sr[k]-sr[k-1]),
         cex=0.75, srt=90, col=c3)
  }
  
  ## Labeling last two bars
  yoff <- 0.03
  xoff <- 4.2
  lines(x=c(5.25+offsets[8]+0.5, 8+offsets[8]-1.15*width),
        y=c(0.8225, sr[8]+0.035), lwd=1.5, col=c1)
  text(x=5+offsets[8]-xoff, y=0.8+yoff, labels="STRING's Total Score",
       col=c1, cex=0.75, adj=c(0,0.5))
  # lines(x=c(7+offsets[7], 9+offsets[9]-1.15*width),
  #       y=c(0.95, sr[9]+0.025), lwd=2, col=c2)
  # text(x=7+offsets[7]-xoff, y=0.95+yoff, labels='EvoWeaver Random Forest',
  #      col=c2, cex=0.75, adj=c(0,0.5))
  plot_heatmap(rocdata)
  
  dev.off(dev.list()['pdf'])
}

StringVsEW <- make_rocs()
actual <- StringVsEW$StringScores$ActualCat <= 3
StringRes <- calc_running_STRINGscore(StringVsEW$StringScores[,3:9], actual)
plot_newFig(StringVsEW, StringRes)
