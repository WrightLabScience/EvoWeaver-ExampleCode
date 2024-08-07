plot_heatmap <- function(vals, textval=vals,
                         cols=colorRampPalette(c('blue', '#F0F0F0', 'orange'))(100),
                         xposrange=c(0,1),
                         yposrange=c(0,1),
                         label_cells=TRUE,
                         border='gray40',
                         textcol='black',
                         valrange=NULL, add=FALSE,
                         upper_only=FALSE, inside_legend=FALSE,
                         xadjustment=0, yadjustment=0, ...){
  TEXT_CUTOFF <- 6L
  cur_mar <- par('mar')
  par(mar=c(0,0,0,0))
  stopifnot('incompatible dimensions'=length(dim(vals)) == length(dim(textval)))
  stopifnot('incompatible dimensions'=all(dim(vals) == dim(textval)))
  if(all(is.na(vals)))
    stop("'vals' must have at least one nonmissing entry")

  if(!add){
    plot.new()
  }
  if(is.null(valrange))
    valrange <- c(min(vals), max(vals))

  if(upper_only){
    vals[lower.tri(vals, diag=TRUE)] <- NA
    textval[lower.tri(textval, diag=TRUE)] <- NA
  }

  # we can just draw as if it's centered and then rescale at the end

  # First getting positions of the rectangles
  hmwidth <- ncol(vals)
  hmheight <- nrow(vals)
  boxwidth <- 1/hmwidth
  allvals <- seq(0,1,by=boxwidth)+xadjustment
  rightvals <- allvals[-1]
  leftvals <- allvals[-length(allvals)]

  boxheight <- 1/hmheight
  allvals <- seq(0,1,by=boxheight)+yadjustment
  bottomvals <- allvals[-1]
  topvals <- allvals[-length(allvals)]

  # Next getting positions of text
  rowannotx <- rep(rightvals[length(rightvals)]+hmwidth*0.0025, hmheight)
  rowannoty <- topvals + (boxheight / 2)
  colannoty <- rep(max(bottomvals)+hmheight*0.0025, hmwidth)
  colannotx <- leftvals + (boxwidth / 2)
  if(is.null(rownames(vals)))
    rowannot <- rep('', hmheight)
  else
    rowannot <- rownames(vals)

  if(is.null(colnames(vals)))
    colannot <- rep('', hmwidth)
  else
    colannot <- colnames(vals)

  if(upper_only){
    rowannot[length(rowannot)] <- ''
    colannot[1] <- ''
  }

  # We'll draw it column by column
  leftvals <- rep(leftvals, each=hmheight)
  rightvals <- rep(rightvals, each=hmheight)
  bottomvals <- rep(rev(bottomvals), length.out=hmheight*hmwidth)
  topvals <- rep(rev(topvals), length.out=hmheight*hmwidth)

  all_color_ranges <- seq(valrange[1], valrange[2], length.out=length(cols)-2)
  all_color_ranges <- c(valrange[1]-1, all_color_ranges, valrange[2]+1)

  bordercols <- rep(border, length(vals))
  posmissing <- which(is.na(vals))
  vals[posmissing] <- mean(vals, na.rm=TRUE)
  all_colors <- vapply(c(vals),
                       \(x) cols[which.min(x >= all_color_ranges)],
                       character(1L))
  vals[posmissing] <- NA
  all_colors[posmissing] <- '#00000000'
  bordercols[posmissing] <- '#00000000'

  # Draw cell annotations
  textposx <- textposy <- 0
  if(label_cells){
    posna <- which(is.na(c(textval)))
    if(!is.character(textval))
      textval <- as.character(textval)
    textval <- c(textval)
    textval[posna] <- ''
    longpos <- nchar(textval) > TEXT_CUTOFF
    textval[longpos] <- substring(textval[longpos], 1, TEXT_CUTOFF)
    textposx <- leftvals + (boxwidth / 2)
    textposy <- bottomvals - (boxheight / 2)
  }

  # Legend plotting
  if(inside_legend){
    m <- 4
    l_boxheight <- (1-m*boxheight) / length(cols)
    l_yp <- seq(0,1-m*boxheight, l_boxheight) + 0.25*m*boxheight
    l_left <- min(leftvals) + boxwidth/4 + boxwidth
    l_right <- min(leftvals) + 3*boxwidth/4 + boxwidth
    l_bottom <- l_yp[-length(l_yp)]
    l_top <- l_yp[-1]
    l_midpoint <- l_top[length(l_top) %/% 2 + 1] - l_boxheight
    l_textxpos <- l_right[1] + 0.05
    l_textypos <- c(l_bottom[2], l_midpoint, l_bottom[length(l_bottom)-1]) + l_boxheight/2
  } else {
    l_boxwidth <- 1 / length(cols)
    l_xp <- seq(0,1,l_boxwidth)
    l_left <- l_xp[-length(l_xp)]
    l_right <- l_xp[-1]
    l_top <- rep(min(topvals)-boxheight/4, length(l_left))
    l_bottom <- l_top-0.03
    l_midpoint <- l_left[length(l_left) %/% 2 + 1] - l_boxwidth
    l_textxpos <- c(l_left[2], l_midpoint, l_left[length(l_left)-1]) + l_boxwidth/2
    l_textypos <- l_bottom[1]-0.05
  }
  # First pass adjustments,
  # this brings everything within a (0,1) range
  sm <- 0.8
  yoffset <- 0.15*sm
  xoffset <- 0.0*sm

  # Second pass adjustments, relocate the whole thing for the user
  uo_x <- xposrange[1]
  uo_y <- yposrange[1]
  us_x <- xposrange[2] - xposrange[1]
  us_y <- yposrange[2] - yposrange[1]

  ## Offsets and adjustments
  leftvals <- (leftvals*sm+xoffset)*us_x+uo_x
  rightvals <- (rightvals*sm+xoffset)*us_x+uo_x
  rowannotx <- (rowannotx*sm+xoffset)*us_x+uo_x
  colannotx <- (colannotx*sm+xoffset)*us_x+uo_x
  textposx <- (textposx*sm+xoffset)*us_x+uo_x
  l_left <- (l_left*sm+xoffset)*us_x+uo_x
  l_right <- (l_right*sm+xoffset)*us_x+uo_x
  l_midpoint <- (l_midpoint*sm+xoffset)*us_x+uo_x
  l_textxpos <- (l_textxpos*sm+xoffset)*us_x+uo_x

  bottomvals <- (bottomvals*sm+yoffset)*us_y+uo_y
  topvals <- (topvals*sm+yoffset)*us_y+uo_y
  rowannoty <- (rowannoty*sm+yoffset)*us_y+uo_y
  colannoty <- (colannoty*sm+yoffset)*us_y+uo_y
  textposy <- (textposy*sm+yoffset)*us_y+uo_y
  l_bottom <- (l_bottom*sm+yoffset)*us_y+uo_y
  l_top <- (l_top*sm+yoffset)*us_y+uo_y
  l_textypos <- (l_textypos*sm+yoffset)*us_y+uo_y



  # Draw heatmap
  rect(leftvals, topvals,
       rightvals, bottomvals,
       col=all_colors, border = bordercols)

  # Draw row and col names
  text(x=rowannotx, y=rowannoty,
       labels=rev(rowannot),
       srt=-30, adj=c(0,0.5), xpd=NA, ...)
  text(x=colannotx, y=colannoty,
       labels=colannot,
       srt=45, adj=c(0,0), xpd=NA, ...)

  # Draw cell labels
  if(label_cells)
    text(x=textposx, y=textposy, labels=textval, col=textcol, cex=0.6,...)

  # Draw legend
  bh <- abs(l_bottom[1] - l_top[1]) / 2
  bw <- abs(l_left[1] - l_right[1]) / 2
  rect(l_left, l_bottom-bh*2, l_right, l_top-bh*2, col=cols, border='#00000000')
  text(x=l_textxpos,
       y=l_textypos,
       srt=0,
       #labels=round(c(valrange[1], mean(valrange), valrange[2]), 2), ...)
       labels=c('-1', '0', '+1'), ...)
  if(inside_legend){
    subset_ticks <- c(seq(1, length(l_bottom), by=12.5), length(l_bottom)-2)
    rect(rep(l_textxpos[1]-4.2*bw/2, length(subset_ticks)), l_bottom[subset_ticks]+bh/2,
         rep(l_textxpos[1]-3*bw/2, length(subset_ticks)), l_bottom[subset_ticks]+bh/2,
         col='black', border='black')
    text(x=l_textxpos-5.5*bw, y=l_textypos[2], srt=90,
         labels=expression(plain(Spearman)~plain(Correlation)~'('*rho*')'),
                      adj=c(0.5,0.5), font=5)
  }
  rect(l_left, min(l_bottom)-bh, l_right, max(l_top)-bh, border='black', lwd=1.5)

  par(mar=cur_mar)
}
plot_fig_2x2 <- function(MainAlgos, EnsembleAlgos, RawScores, filename,
                         isComplex=FALSE,
                         xadjustment=-0.025, yadjustment=-0.015, ...){
  curpar <- par(no.readonly = TRUE)
  aa <- c(MainAlgos,EnsembleAlgos)
  aacols <- vapply(aa, \(x) x$Color, numeric(1L))
  aa <- aa[order(aacols)]
  aacols <- aacols[order(aacols)]
  aalty <- unlist(lapply(rle(aacols)$lengths, seq_len))
  #aalty <- (c(1,2,4,0))[aalty]
  names(aalty) <- names(aacols)
  aalwd <- aalty
  aalwd[] <- 1L
  aalwd[names(aalwd) %in% names(EnsembleAlgos)] <- 2L
  .plot_curve <- function(isROC, isSmall, ylabel, ...){
    cols <- c('#D81B60', '#1E88E5', '#FFC107', '#5CDC63', 'black')
    colList <- list(red=c('#E82B70','#C80B50','#A80030','#780000'),
                    blue=c('#1E88E5','#0E68C5','#0048A5','#0028C5'),
                    yellow=c('#FFC107','#DFA100','#BF8100'),
                    other=c('#5CDC63','#3CBC83','#7CFC83'),
                    black=c('black', 'gray60', 'gray80'))
    xlab <- ifelse(isROC, "False positive rate", "Recall")
    ylab <- ifelse(isROC, "True positive rate", "Precision")
    xmaxv <- ifelse(isSmall, 0.05, 1)
    #### Plotting ROCs
    nowpar <- par(no.readonly = TRUE)
    if(isSmall){
      v <- par('plt')
      # left, right, bottom, top
      lrbt <- c(0.30, 0.47, 0.63, 0.77)
      par(fig=lrbt, new=TRUE, mar=c(0,0,0,0))
      plot(c(0,1), c(0,1), type='l', xaxs='i', yaxs='i',
           col='black', lty=2, lwd=1, xlab='', ylab='',
           ylim=c(0, 1), xlim=c(0, 0.02),
           main='', axes=FALSE, frame=TRUE,
           oma=c(0,0,0,0))
      # x axis
      axis(1, seq(0,0.02,0.005),
           labels = c("0.00", '', '0.01', '', '0.02'),
           cex.axis=1, tck=-0.05, mgp=c(3,0.5,0), cex.lab=1)
      axis(2, seq(0,1.0,0.2),
           labels= c("0.0", '', '0.4', '', '0.8', ''),
           cex.axis=1, tck=-0.05, mgp=c(3,0.5,0), cex.lab=1)
    } else {
      plot(c(0,1), c(0+((!isROC)*0.5),1-((!isROC)*0.5)), type='l', xaxs='i', yaxs='i',
           col='black', lty=2, lwd=1, ylab='', xlab='',
           ylim=c(ifelse(isROC, 0, 0.4), 1), xlim=c(0, xmaxv),
           main=ifelse(isROC, 'Receiver Operating Characteristic (ROC)', "Precision-Recall Curve (PRC)"),
           mar=c(0, 0, 0, 0) + 0.1,
           oma=c(0, 0, 0, 0) + 0.1,
           mgp=c(0,0.5,0),
           ...)
      title(xlab=xlab, mgp=c(2,1,0))
      title(ylab=ifelse(ylabel, ylab, ''), mgp=c(2,1,0))
    }
    for (i in seq_along(aa)){
      pn <- names(aa)[i]
      p <- aa[[i]]
      #lw <- aalwd[pn]
      lw <- ifelse(aacols[pn]!=5, 1, 1.5)
      #lty <- aalty[pn]
      lty <- (c(1,5,3,4))[aalty[pn]]
      col <- colList[[aacols[pn]]][aalty[pn]]
      if(isROC){
        yv <- p$TPR[p$FPR <= xmaxv]
        xv <- p$FPR[p$FPR <= xmaxv]
        lines(p$FPR, p$TPR, lwd=lw, lty=lty, col=col)
      }
      else{
        yv <- p$Precision[p$TPR <= xmaxv]
        xv <- p$TPR[p$TPR <= xmaxv]
        lines(p$TPR, p$Precision, lwd=lw, lty=lty, col=col)
      }
    }
    #par(nowpar)
    if(!isSmall & isROC){
      legnames <- c(names(aa), "Random Guessing")
      aa[["Random Guessing"]] <- list(AUROC=0.5,AUPRC=0.5)
      aacols[["Random Guessing"]] <- 5
      aalty[["Random Guessing"]] <- 1
      maxwidth <- max(strwidth(legnames)) + strwidth('.')
      rocadd <- vapply(names(aa), \(x){
        return(sprintf("(AUROC: %0.3f, AUPRC: %0.3f)",
                       aa[[x]]$AUROC, aa[[x]]$AUPRC))
      }, character(1L))
      aurocs <- vapply(names(aa), \(x) sprintf("%0.3f", aa[[x]]$AUROC), character(1L))
      auprcs <- vapply(names(aa), \(x) sprintf("%0.3f", aa[[x]]$AUPRC), character(1L))
      # legnames <- vapply(legnames, \(x){
      #   while(strwidth(x) < maxwidth){
      #     x <- paste0(x, ' ')
      #   }
      #   x
      # }, character(1L))
      # legnames <- paste(legnames, rocadd)
      legnames <- c("Phylogenetic Profiling", legnames[1:4], '',
                    "Phylogenetic Structure", legnames[5:7], '',
                    "Gene Organization", legnames[8:10], '',
                    "Sequence Level", legnames[11:12], '',
                    "Ensemble Methods", legnames[13:15], '',
                    "Other", legnames[16])
      pm <- rep(FALSE, length(legnames))
      pm[c(1,6,7,11,12,16,17,20,21,25,26)] <- TRUE

      rocBold <- prcBold <- rep(1, length(pm))
      aurocStr <- character(length(pm))
      aurocStr[!pm] <- aurocs
      aurocStr[1] <- "AUROC"
      suppressWarnings(rocBold[c(1,which.max(as.numeric(aurocStr)))] <- 2)
      rocBold[c(2,8,13,18)] <- prcBold[c(2,8,13,18)] <- 2
      auprcStr <- character(length(pm))
      auprcStr[!pm] <- auprcs
      auprcStr[1] <- "AUPRC"
      suppressWarnings(prcBold[c(1,which.max(as.numeric(auprcStr)))] <- 2)

      #pm[24:26] <- TRUE
      allcolors <- rep("#00000000", length(legnames))
      allcolors[!pm] <- vapply(names(aa),
                               \(x) colList[[aacols[x]]][aalty[x]], character(1L))
      allcolors[pm] <- "#00000000"
      #$allcolors[26] <- 'black'
      pchn <- rep(1, length(allcolors))
      pchn[pm] <- 0
      pchn[!pm] <- aalty
      alllty <- rep(1L, length(pchn))
      alllty[!pm] <- vapply(aalty, \(x) (c(1L,5L,3L,4L))[x], integer(1L))
      alllty[27] <- 2L
      flb <- rep(1L, length(pchn))
      flb[pm] <- 2L
      #flb[26] <- 1L
      txtcex <- rep(1L, length(flb))
      txtcex[legnames==''] <- 0.33
      textcols <- rep('black', length(alllty))
      textcols[c(1,7,12,17)] <- c('#D81B60', '#2B6DA8', '#E0A608', '#45A649')
      h <- -0.050
      legend('topright', inset=c(-0.700,h),
             legend=legnames,
             col=allcolors,
             lty=alllty,
             lwd=c(rep(1.5,26), 1),
             text.font=flb,
             bty='n',
             y.intersp=txtcex,
             text.col=textcols,
             xpd=NA,
             seg.len=2.5)
      legend('topright', inset=c(-1.04,h),
             legend=aurocStr,
             text.font=rocBold,
             adj=c(0.5,0.5),
             bty='n',
             y.intersp=txtcex,
             text.col='black',
             xpd=NA)
      legend('topright', inset=c(-1.315,h),
             legend=auprcStr,
             adj=c(0.5,0.5),
             text.font=prcBold,
             bty='n',
             y.intersp=txtcex,
             text.col='black',
             xpd=NA)
    }
  }
  #### Plotting Heatmap
  heatmapalgos <- names(MainAlgos)
  ll <- length(heatmapalgos)
  corrvals <- matrix(nrow=ll, ncol=ll)
  colnames(corrvals) <- rownames(corrvals) <- heatmapalgos
  corrvalsNeg <- corrvals
  for(algo1 in heatmapalgos){
    for(algo2 in heatmapalgos){
      corrvals[algo1,algo2] <- corrvals[algo2, algo1] <- cor(RawScores[RawScores$isTP,algo1],
                                                             RawScores[RawScores$isTP,algo2],
                                                             use='pairwise',
                                                             method='spearman')
      corrvalsNeg[algo1,algo2] <- corrvalsNeg[algo2, algo1] <- cor(RawScores[!RawScores$isTP,algo1],
                                                                   RawScores[!RawScores$isTP,algo2],
                                                                   use='pairwise',
                                                                   method='spearman')
    }
  }

  avgheat <- (corrvalsNeg + corrvals) / 2
  displaynums <- matrix(sprintf("%.2f", avgheat), nrow=nrow(avgheat))
  displaynums[displaynums=='-0.00'] <- "0.00"
  displaynums[is.na(avgheat)] <- ''

  pdf(filename, width=4.3*2, height=4.3*2, onefile=TRUE)

  layout(matrix(c(3,0,1,2), nrow=2, byrow=TRUE), widths=c(1,1), heights=c(1,1))
  par(mar=c(5,4,1,0.5)+0.1)
  .plot_curve(FALSE,FALSE,TRUE, ...)

  # bottom half of geyser, top half of temps, white in the middle
  geysercols <- c("#1c6376", '#008585', "#4E998A","#7CAC94","#A4BEA5","#C3D4B1",'white',
                  "#eee2a9","#eecd93","#edb588","#ea9c89","#DB6577","#a42f54")

  plot_heatmap(avgheat, displaynums, upper_only = TRUE, valrange=c(-1,1),
               yposrange = c(0,0.86), textcol = 'black',
               col=colorRampPalette(geysercols)(100),
               inside_legend=TRUE,
               xadjustment=xadjustment,
               yadjustment=yadjustment, ...)
  # Adding boxes around same algorithm heatmaps
  lbx <- 0.0675+xadjustment+0.005
  tby <- 0.7915+yadjustment+0.005
  blw <- 0.0666
  blh <- 0.0574
  .draw_rect <- function(l,r,t,b,col,lty=1){
    rect(lbx+l*blw,tby-b*blh,lbx+r*blw,tby-t*blh, border=col, lwd=2, lty=lty)
  }
  # PA lines
  .draw_rect(0,3,0,0, '#D81B60')
  .draw_rect(3,3,3,0, '#D81B60')
  .draw_rect(2,3,3,3, '#D81B60')
  .draw_rect(2,2,2,3, '#D81B60')
  .draw_rect(1,2,2,2, '#D81B60')
  .draw_rect(1,1,1,2, '#D81B60')
  .draw_rect(0,1,1,1, '#D81B60')
  .draw_rect(0,0,0,1, '#D81B60')

  # DM lines
  .draw_rect(4,6,4,4, '#2B6DA8')
  .draw_rect(6,6,4,6, '#2B6DA8')
  .draw_rect(5,6,6,6, '#2B6DA8')
  .draw_rect(5,5,5,6, '#2B6DA8')
  .draw_rect(4,5,5,5, '#2B6DA8')
  .draw_rect(4,4,4,5, '#2B6DA8')

  # GeneDistance lines
  .draw_rect(7,9,7,7, '#DFA100')
  .draw_rect(9,9,7,9, '#DFA100')
  .draw_rect(8,9,9,9, '#DFA100')
  .draw_rect(8,8,8,9, '#DFA100')
  .draw_rect(7,8,8,8, '#DFA100')
  .draw_rect(7,7,7,8, '#DFA100')

  .draw_rect(10,11,10,11, '#45A649')

  .plot_curve(TRUE,FALSE,TRUE, ...)
  .plot_curve(TRUE,TRUE,FALSE, ...)
  par(curpar)
  invisible(dev.off())
}
