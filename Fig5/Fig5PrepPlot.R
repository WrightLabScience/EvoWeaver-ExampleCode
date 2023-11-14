library(phytools)
library(SynExtend)
library(circlize)
library(dendextend)
source('PrepPlotGainLoss.R')
load('ModulePredsAllPairs.RData')
load('B3GNT5_ST6GAL1_plottingdata.RData', v=T)

outdir <- './'
outname <- '5_FigDiscoveryGraphs.pdf'
outviolinname <- '5_FigDiscoveryViolin.pdf'

ptsize <- 12
legendcex <- 0.75
axsize <- 0.9
titlecex <- 1
legendboxcol <- 'black'
legendbgcol <- 'white'
w <- 7.08
## All the titles
mjcol1 <- -25.5
mjcol2 <- 23.25
mjrow1 <- 69
mjrow2 <- -13.5
cvec <- c('#45A649','#D81B60', '#1E88E5','#FFC107')
pdf(file.path(outdir, outname), width=w, height=w, onefile=TRUE, pointsize=ptsize)
layout(matrix(c(rep(0,3),rep(0,5),
                rep(1,3),0,0,rep(2,3),
                rep(1,3),rep(0,5),
                rep(0,8),
                rep(3,4),rep(4,4)), nrow=5, byrow=TRUE),
       heights = c(0.05,0.95,0.05,0.025, 1),
       widths = c(1,1,1,0.05,0.05,1,1,1))

.plotColocAndLabels <- function(){
  # colocalization
  y <- read.tree(ofile)
  scale <- 1/5
  m <- 80
  ncd <- colocdiffs
  colocdiffs[colocdiffs > m] <- m+10
  #colocdiffs <- colocdiffs / 10
  #colocdiffs <- log(colocdiffs+1)
  #colocdiffs[colocdiffs>8] <- 8
  plotTree.wBars(y, colocdiffs,
                 tip.labels=FALSE,
                 col = c(cvec[3],cvec[3],cvec[4])[colocdirs+1],
                 scale=scale, fsize=0.8)
  par(xpd=NA)
  numv <- 9L+1
  lastval <- numv
  numv <- (seq_len(numv)-1) * 10
  bc <- c(24.17, -1)
  for(i in seq_along(numv)+1){
    if(i < lastval+1){
      lines(x=bc[1]+c(numv[i-1]*scale, numv[i-1]*scale),
            y=c(bc[2], length(colocdiffs)+0.5), lty=3, col='gray30')
    }
    # rect(xleft=bc[1]+numv[i-1]*scale,
    #      xright=bc[1]+numv[i]*scale,
    #      ytop=bc[2]+0.5,
    #      ybottom=bc[2],
    #      col='gray', border = 'black')
    if(i < (lastval+1)){
      rect(xleft=bc[1]+numv[i-1]*scale,
           xright=bc[1]+numv[i-1]*scale,
           ytop=bc[2]+0,
           ybottom=bc[2]-0.3,
           col='gray', border = 'black')
      text(x=bc[1]+numv[i-1]*scale, y=bc[2]-0.75,
           labels=as.character(numv[i-1]),
           adj=c(0,0.5), cex=axsize, srt=-45)
    }
  }
  lines(x=c(bc[1], bc[1]+((lastval-1)*scale*10)),y=c(rep(bc[2]+0.25,2)),
            lty=1, col='black', lwd=1.5)
  text(x=bc[1]+(m/2)*scale, y=bc[2]-4.5,
       labels=substitute(paste(bold('Number of Genes Apart'))),
       cex=axsize)

  # hardcoded line breaks
  pos_break <- names(colocdiffs)[which(colocdiffs>m)]
  vals <- ncd[pos_break]
  pos_break <- vapply(pos_break, \(l) which(l == labels(y)), integer(1L))
  #pos_break <- length(colocdiffs) - pos_break
  pos_break <- c(pos_break, -0.75)
  l <- 0.5
  xgap <- 2
  xw <- 2
  xp <- m+2
  for(i in seq_along(pos_break)){
    lines(x=c(bc[1]+(xp+1)*scale, bc[1]+(xp+xw+1)*scale),
          y=c(pos_break[i]-l,pos_break[i]+l),
          col='white', lwd=1)
    lines(x=c(bc[1]+xp*scale, bc[1]+(xp+xw)*scale),
          y=c(pos_break[i]-l,pos_break[i]+l),
          col='black', lwd=1)
    lines(x=c(bc[1]+(xp+xgap)*scale, bc[1]+(xp+xgap+xw)*scale),
          y=c(pos_break[i]-l,pos_break[i]+l),
          col='black', lwd=1)
    if(pos_break[i] != 0)
      text(x=bc[1]+(m+3)*scale, y=pos_break[i]-0.75,
           labels=as.character(vals[i]),
           cex=axsize, adj=c(0,1), xpd=NA)
  }

  # direction legend
  dirpos <- -1
  bh <- bc[2]+0.5
  arrowlen <- 5.5
  yoffset <- -4
  arrowhgt <- 2.5
  ac1 <- 0.5
  ac2 <- 1+ac1+2*arrowlen
  rect(xleft=dirpos, xright=dirpos+ac2+2*arrowlen+ac1,
       ybottom=bh+yoffset-arrowlen, ytop=bh+0.75,
       lty=1, border=legendboxcol, col=legendbgcol)
  gene1 <- 'B3GNT5'
  gene2 <- 'ST6GAL1'
  ttal <- 0.75*arrowlen
  hahgt <- arrowhgt / 2
  text(x=dirpos+ac1+1.1*arrowlen, y=bh-0.5, labels='Transcriptional Direction',
       cex=legendcex, adj=c(0,0.5))
  text(x=dirpos+ac1+arrowlen, y=bh+yoffset-1.525*arrowhgt,
       labels='Same', cex=legendcex, adj=c(0.5,0.5), col='black')
  text(x=dirpos+ac2+arrowlen, y=bh+yoffset-1.525*arrowhgt,
       labels='Opposite', cex=legendcex, adj=c(0.5,0.5), col='black')

  polygon(x=c(0,ttal,ttal,arrowlen,ttal,ttal,0)+dirpos+ac1,
          y=c(hahgt,hahgt,arrowhgt,0,-arrowhgt,-hahgt,-hahgt)+bh+yoffset,
          col=cvec[3])
  polygon(x=c(0,ttal,ttal,arrowlen,ttal,ttal,0)+dirpos+arrowlen+ac1,
          y=c(hahgt,hahgt,arrowhgt,0,-arrowhgt,-hahgt,-hahgt)+bh+yoffset,
          col=cvec[3])
  polygon(x=c(0,ttal,ttal,arrowlen,ttal,ttal,0)+dirpos+ac2,
          y=c(hahgt,hahgt,arrowhgt,0,-arrowhgt,-hahgt,-hahgt)+bh+yoffset,
          col=cvec[4])
  polygon(x=arrowlen-c(0,ttal,ttal,arrowlen,ttal,ttal,0)+dirpos+arrowlen+ac2,
          y=c(hahgt,hahgt,arrowhgt,0,-arrowhgt,-hahgt,-hahgt)+bh+yoffset,
          col=cvec[4])
  xoff <- 2.5

  text(x=dirpos+ac1+xoff,y=bh+yoffset, substitute(paste(bold(gene1))), col='white', cex=legendcex, adj=c(0.5,0.5))
  text(x=dirpos+ac1+xoff+arrowlen+0.1,y=bh+yoffset, substitute(paste(bold(gene2))), col='white', cex=legendcex, adj=c(0.5,0.5))
  text(x=dirpos+ac2+xoff,y=bh+yoffset, substitute(paste(bold(gene1))), col='black', cex=legendcex, adj=c(0.5,0.5))
  text(x=dirpos+ac2+xoff+arrowlen+0.4,y=bh+yoffset, substitute(paste(bold(gene2))), col='black', cex=legendcex, adj=c(0.5,0.5))
  text(x=mjcol1, y=mjrow1,
        labels=substitute(paste(bold('Phylogenetic Profiling'))),
       cex=titlecex, col='#D81B60', xpd=NA)
  #par(xpd=NA)
  text(x=mjcol2, y=mjrow1,
       labels=substitute(paste(bold('Gene Organization'))),
       cex=titlecex, col='#E0A608', xpd=NA)
  #par(xpd=NA)
  text(x=mjcol2, y=mjrow2,
       labels=substitute(paste(bold('Sequence Level'))),
       cex=titlecex, col='#45A649', xpd=NA)
  #par(xpd=NA)
  text(x=mjcol1, y=mjrow2,
       labels=substitute(paste(bold('Phylogenetic Structure'))),
       cex=titlecex, col='#317EC2', xpd=NA)

  text(x=-0.6, y=25, srt=90, cex=axsize,
       labels=substitute(paste(bold('Subset of Species Tree'))))
  par(mar=c(5.1,4.1,4.1,2.1), xpd=FALSE)
}

.plotNVDT <- function(){
  pchs <- c(19,15,17)
  x <- sqrt(abs(nvdtCorr60vecs[1,]))*sign(nvdtCorr60vecs[1,])
  y <- sqrt(abs(nvdtCorr60vecs[2,]))*sign(nvdtCorr60vecs[2,])
  plot(x=x, y=y, pch=rep(pchs, each=20), col=rep(cvec[4:2], each=20),
       xaxt='n', yaxt='n', xlim=c(-1.5,2), ylim=c(-1.5,2),
       xlab='', ylab='', mar=c(0,0,0,0))
  axis(side=2, at=seq(-2,2,by=0.5),
       tck=-0.01, mgp=c(3,0.5,0), las=3)
  axis(side=1, at=seq(-2,2,by=0.5),
       tck=-0.01, mgp=c(3,0.25,0))
  mtext("B3GNT5 Sequence Features",
        side=2, cex=axsize*par('cex'), line=1.6)
  mtext("SDT6GAL1 Sequence Features",
        side=1, cex=axsize*par('cex'), line=1.25)
  # NVDT LoBF
  lobf <- lm(y~x)
  abline(a=lobf$coefficients[1], b=lobf$coefficients[2],
         col=cvec[1], lty=2, lwd=2, untf=FALSE)
  text(x=1.7, y=0.65, parse(text=paste0("R^2 == ", round(summary(lobf)$r.squared, 3L))),
       col=cvec[1], cex=axsize)

  legend(x=-1.5,y=2,
         legend=c('Residue Counts','Residue Positions',
                  'Residue Variance', 'Best Fit'),
         pch=c(pchs, NA),
         lty=c(0,0,0,2),
         col=cvec[4:1],
         pt.cex=legendcex,
         lwd=c(0,1),
         cex=legendcex, box.col=legendboxcol, box.lty=1,
         bg=legendbgcol,
         seg.len = 1.5, ncol=1, x.intersp=0.75,
         y.intersp=1)
}

.plotMT <- function(){
  # MirrorTree
  labs <- names(colocdiffs)
  treeA <- subset(M1$Tree, labs)
  treeG <- subset(M2$Tree, labs)
  tanglegram(treeA, treeG, edge.lwd=1, lwd=1,
             common_subtrees_color_branches=TRUE,
             common_subtrees_color_lines_default_single_leaf_color='grey',
            sort=TRUE, just_one = FALSE, lab.cex=0.8,
            margin_inner=0, margin_outer=2,
            color_lines='grey')
  par(xpd=NA)
  y <- rep(18, 2)
  x <- c(-1.33, 0.62)
  text(x=x[1], y=y[1],
       labels = substitute(paste(bold("M1 Gene Tree"))),
       srt=90, cex=0.8)
  par(xpd=NA)
  text(x=x[2], y=y[2],
       labels = substitute(paste(bold("M2 Gene Tree"))),
       srt=-90, cex=0.8)
  par(xpd=FALSE)
}

.plotRPMT2ax <- function(){
  RPCVecs <- RPCVecs[abs(RPCVecs[,1]) < 1e4 & abs(RPCVecs[,2]) < 1e4,]
  x <- RPCVecs[,1]
  y <- RPCVecs[,2]
  c1 <- as.matrix(Cophenetic(M1$Tree))
  c2 <- as.matrix(Cophenetic(M2$Tree))
  subs <- intersect(rownames(c1), rownames(c2))

  c1 <- c1[subs,subs]
  c2 <- c2[subs,subs]
  c1 <- c1[upper.tri(c1)]
  c2 <- c2[upper.tri(c2)]
  par(mar=c(5,4,4,4)+0.3)
  rpcol <- '#d29d00'
  mtcol <- '#156dba'
  rpaxcol <- '#a87e00'
  mtaxcol <- '#115795'
  mtcol <- mtaxcol
  rpcol <- rpaxcol

  # MT data and line
  plot(x=c1, y=c2, col='#317EC208', pch=20, cex=0.5,
       axes=FALSE, xlab='', ylab='',
       ylim=c(0,1), xlim=c(0,1))
  axis(side=4, at=pretty(c(-0.1, range(c2), 1.05)),
       col=mtaxcol, col.axis=mtaxcol, mgp=c(3,0.5,0), tck=-0.01, las=3)
  axis(side=3, at=pretty(c(-0.1, range(c1), 1.05)),
       col=mtaxcol, col.axis=mtaxcol, mgp=c(3,0.25,0), tck=-0.01)
  mtext("B3GNT5 Pairwise Distances",
        col=mtaxcol, side=4, cex=axsize*par('cex'), line=1.6)
  mtext("SDT6GAL1 Pairwise Distances",
        col=mtaxcol, side=3, cex=axsize*par('cex'), line=1.1)

  lobf <- lm(c2~c1)
  abline(a=lobf$coefficients[1], b=lobf$coefficients[2],
         col=mtcol, lty=2, lwd=2)
  text(x=0.50, y=0.965, parse(text=paste0("R^2 == ", round(summary(lobf)$r.squared, 3L))),
       col=mtcol, cex=axsize)

  par(lheight=0.7)
  legend(x=0.7,y=0.29,
         legend=c('Projected\nDistance', 'Patristic\nDistance',
                  'RPCT Fit', 'MirrorTree Fit'),
         pch=c(rep(19,2), rep(NA,2)),
         lty=c(rep(0,2), rep(2,2)),
         col=c(cvec[4:3], mtcol, rpcol),
         pt.cex=legendcex,
         lwd=c(rep(0,2), rep(1,2)),
         cex=legendcex, box.col=legendboxcol, box.lty=1,
         bg=legendbgcol,
         seg.len = 1.5, ncol=1, x.intersp=0.75,
         y.intersp=c(1,1.6,0.5,1,0))
  par(lheight=1)

  # RP data and line
  par(new=TRUE)
  plot(x=x,y=y,col=rpcol, pch=1, cex=0.95,
       axes=FALSE, xlab='', ylab='', bty='n',
       xlim=c(-400,800), ylim=c(-400,800))
  points(x=x,y=y,col=cvec[4],pch=20,cex=0.9)
  axis(side=2, at=c(-450,pretty(range(y)), 850),
       col=rpaxcol, col.axis=rpaxcol, tck=-0.01, mgp=c(3,0.5,0), las=3)
  axis(side=1, at=pretty(range(x*1.2)),
       col=rpaxcol, col.axis=rpaxcol, tck=-0.01, mgp=c(3,0.25,0))
  mtext("B3GNT5 Projected Distances",
        col=rpaxcol, side=2, cex=axsize*par('cex'), line=1.6)
  mtext("SDT6GAL1 Projected Distances",
        col=rpaxcol, side=1, cex=axsize*par('cex'), line=1.25)
  lobf <- lm(y~x)
  abline(a=lobf$coefficients[1], b=lobf$coefficients[2],
         col=rpcol, lty=2, lwd=2)
  text(x=650, y=420, parse(text=paste0("R^2 == ", round(summary(lobf)$r.squared, 3L))),
       col=rpcol, cex=axsize)

  points(483,-130, col=rpcol, pch=1, cex=0.9)
  points(483,-230, col='#317EC2', pch=19, cex=0.9)

  par(mar=c(5,4,4,2)+0.1, new=FALSE)

}

.plotGL <- function(){
  # gain/loss
  ptcex <- 0
  tree1 <- PrepPlotGainLoss(subspec, attnum=1L,
                            pch=21, ptcex=ptcex,
                            lwd=c(0.5,0.1)*3,
                            cols=c(cvec[2:1], 'black'))
  tree2 <- PrepPlotGainLoss(subspec, attnum=2L,
                            pch=20, ptcex=ptcex, lty=c(2,0,1),
                            lwd=c(0.5,0.5)*3,
                            col=c(cvec[4:3], 'black'))

  r <- 0.8
  dth <- 0.97
  cm <- 0.5
  par(xpd=NA)
  circos.par(canvas.xlim=c(-r,r), canvas.ylim=c(-r,r),
             circle.margin=cm, cell.padding=c(0,0))
  circlize_dendrogram(tree1, labels=FALSE,
                      dend_track_height = dth)
  legend(x=0.49,y=-0.635,
         legend=c(paste0('+', M1$GeneName),
                  paste0('-', M1$GeneName),
                  '',
                  paste0('+', M2$GeneName),
                  paste0('-', M2$GeneName)),
         fill=c(cvec[1:2], NA, cvec[3:4]),
         border=c(rep('black',2), '#00000000', rep('black',2)),
         cex=legendcex, box.col=legendboxcol, box.lty=1,
         bg=legendbgcol,
         ncol=1, x.intersp=0.5,
         y.intersp=c(0.5,0.75,0.2,0.75,0.5))
  lines(c(0.5,0.85), rep(-0.778,2), lty=2, lwd=1, col='black')
  rect(xleft=0.46, xright=0.51,
       ytop=-0.63, ybottom=-0.68, col='white', border='white')
  rect(xleft=0.51, xright=0.57,
       ytop=-0.61, ybottom=-0.64, col='white', border='white')
  rect(xleft=0.8213,xright=0.84,
       ybottom=0, ytop=0.05, col='white', border = NA)
  par(new=TRUE)
  circos.par(canvas.xlim=c(-r,r), canvas.ylim=c(-r,r),
             circle.margin=cm, cell.padding=c(0,0))
  circlize_dendrogram(tree2, labels=FALSE,
                      dend_track_height = dth)

  v1 <- rapply(subspec, \(x) attr(x, 'FitchState')[1])
  v2 <- rapply(subspec, \(x) attr(x, 'FitchState')[2])
  names(v1) <- names(v2) <- labels(subspec)
  # outer labels
  mh <- attr(tree1, 'height')
  rw <- 4
  par(new=TRUE)
  r <- 0.9
  circos.par(canvas.xlim=c(-r,r), canvas.ylim=c(-r,r),
             circle.margin=cm, cell.padding=c(0,0))
  circos.initialize(1, xlim=c(0,length(v1)))
  circos.track(ylim=c(0,mh+4), track.height=0.3, bg.border=NA, panel.fun=\(x,y){
    l <- length(v1)
    la <- labels(tree1)
    nh <- mh-8
    circos.rect((1:l)-1, rep(nh,l),
                1:l, rep(nh+rw,l),
                col=c(cvec[2:1])[v1[la]+1],
                border=FALSE)
    circos.rect((1:l)-1, rep(nh+rw,l),
                1:l, rep(nh+rw*2,l),
                col=c(cvec[4:3])[v2[la]+1],
                border=FALSE)

  })
  circos.clear()
  text(x=-0.05, y=-0.075, srt=0, cex=legendcex,
       labels=substitute('Species Tree'))
  lines(x=c(0.905,0.973), y=c(0.015,0.016), col='black', lwd=1)
  lines(x=c(0.905,0.973), y=c(0.00,-0.001), col='black')
  arrows(x0=rep(0.62,2), y0=rep(0.2,2),
         x1=c(0.11,0.45), y1=c(-0.04,-0.0525),
         col='black', length=0.03, lwd=1.5)
  arrows(x0=0, y0=0.695,
         x1=-0.15, y1=0.33,
         col='black', length=0.03, lwd=1.5)
  par(lheight=0.7)
  text(x=0.7, y=0.23, "Simultaneous\ngains", col='black', cex=legendcex)
  text(x=0.06, y=0.72, "Simultaneous\nloss", col='black', cex=legendcex)
  par(xpd=FALSE, new=FALSE, mar=c(5.1,4.1,4.1,2.1), lheight=1)
}

.plotViolins <- function(){
  require(sm)
  labs <- c("G/L Correlation", "G/L Distance", "P/A Jaccard", "P/A MI",
            "RP ContextTree", "RP MirrorTree", "Tree Distance",
            "Gene Distance", "Transcription MI", "Moran's I",
            "Sequence Info", "Gene NV")
  layout(matrix(c(0,1), nrow=1),
         widths=c(0.2,1), heights=0.9)
  AllPairs <- AllPairs[AllPairs$NoKOOverlapCategory!=0 & !AllPairs$HasComplex,]
  d <- AllPairs[,c(3:6,8,10:16)]
  d[,4:5] <- d[,4:5]*AllPairs$PAPV
  ResultsLine <- AllPairs[which(AllPairs$Mod1=='M133_1_1' & AllPairs$Mod2=='M339_4_3'),]
  r <- ResultsLine[,c(3:6,8,10:16)]
  r[,4:5] <- r[,4:5]*ResultsLine$PAPV
  for(i in seq_len(ncol(d))){
    r[,i] <- r[,i] - min(d[,i], na.rm=TRUE)
    d[,i] <- d[,i] - min(d[,i], na.rm=TRUE)
    r[,i] <- r[,i] / max(d[,i], na.rm=TRUE)
    d[,i] <- d[,i] / max(d[,i], na.rm=TRUE)
  }
  rearrv <- c(7,6,4,5,2,1,3,9,11,10,12,8)
  d <- d[,rearrv]
  r <- r[,rearrv]
  interspacing <- 1
  plot(NULL, xlim=c(0,1), ylim=c(0, interspacing*(length(rearrv))+1),
       axes=FALSE, xlab='', ylab='', mar=c(0,6,0,0))
  axis(2, at=interspacing*seq_along(rearrv), labels = rev(labs),
       cex.axis=legendcex*par('cex'), las=2, xpd=NA, tck=-0.01, mgp=c(3,0.25,0))

  totalarea <- 0.25 / 4
  rw <- 0.04 / 3
  allcols <- rep(cvec[c(1,4,3,2)], times=c(2,3,3,4))
  maxcol <- ncol(d) + 1L
  for(i in seq_along(rearrv)){
    v <- d[,maxcol-i]
    v <- v[!is.na(v)]
    den <- sm.density(v, 0.025, display='none', verbose=0)
    xv <- den$eval.points
    yv <- den$estimate
    crit_points <- quantile(v, c(0.025,0.5,0.975), na.rm=TRUE)
    curarea <- DescTools::AUC(xv,yv)
    yv <- yv*(totalarea/curarea)
    polygon(x=c(xv,rev(xv)), y=c(i*interspacing+yv,rep(i, length(yv))),#rev(i*interspacing-yv)),
            col=allcols[i], border='black', lwd=0.5)
    points(x=crit_points[2], y=i*interspacing-0.05, pch=20, col='grey', cex=1.25)
    lines(x=c(-0.05,1), y=c(i,i), col='grey30', xpd=NA)
    rect(xleft=crit_points[1], xright=crit_points[3],
         ybottom=i*interspacing-rw-0.05,ytop=i*interspacing+rw-0.05, col='black')
    lines(x=rep(r[1,maxcol-i],2),
           y=c(i*interspacing+0.5, y1=i*interspacing-0.5), col='grey40',
          lwd=3, lty="11")
  }
  xpos <- 0.25
  ypos <- 0.4
  legend(x=xpos, y=ypos, legend=c("Inner 95% Range",
                             "Mean Score",
                             "B3GNT5/ST6GAL1 Score"),
         lty=c('solid',NA,'11'),
         lwd=c(2,0,2),
         col=c('black', 'black', 'grey40'),
         xpd=NA, cex=legendcex,
         box.col=legendboxcol, box.lty=1,
         bg=legendbgcol,
         y.intersp=0.8)
  points(xpos+0.09, ypos-0.94, pch=19, col='grey', xpd=NA)
  layout(1)
}

.plotGL()
.plotColocAndLabels()
.plotRPMT2ax()
.plotNVDT()
dev.off()
layout(1)

pdf(file.path(outdir, outviolinname), width=w*2/3, height=w*0.8, onefile=TRUE, pointsize=ptsize)
.plotViolins()
dev.off()
