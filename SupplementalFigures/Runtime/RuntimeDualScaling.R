load(file.path(datadir, 'EvoWeaver_RuntimeBenchmark.RData'))
outfile <- 'SXX_RuntimePlot.pdf'

benchmark_timings <- benchmark_timings[1:400,]

npairs <- benchmark_timings[,5]
npairs <- npairs / max(npairs)

results <- matrix(nrow=4, ncol=ncol(benchmark_timings)-5L)
set.seed(179L)
for(i in seq(6,ncol(benchmark_timings))){
  res <- NULL
  ntaxa <- benchmark_timings[,1]
  ntaxa <- ntaxa / max(ntaxa)
  timeval <- benchmark_timings[,i] / max(benchmark_timings[,i])
  res <- nlsLM(timeval ~ x*(ntaxa**a) * (npairs**b) + z,
             start=list(a=1,b=1,x=5,z=10),
             lower=unlist(list(a=0,b=0,x=0.001,z=-Inf)),
             control=nls.lm.control(maxiter=1000))
  v <- summary(res)$coefficients[c("a","b","x","z"),1]
  results[,i-5L] <- v
}

colnames(results) <- colnames(benchmark_timings)[-(1:5)]
rownames(results) <- c("TaxaExp", "PairsExp", "TP_Mult", "Intercept")
results

## example runtime
NUM_TAXA <- max(benchmark_timings[,1])
NUM_PAIRS <- max(benchmark_timings[,5])
exData <- benchmark_timings[benchmark_timings[,1]==NUM_TAXA &
                              benchmark_timings[,5]==NUM_PAIRS,]
example_point <- apply(exData[,-(1:5)], 2, mean)
exData <- benchmark_timings[benchmark_timings[,1]==NUM_TAXA &
                                      benchmark_timings[,5]==1,]
example_point2 <- apply(exData[,-(1:5)], 2, mean)

name_mapping <- c("P/A Jaccard", "G/L MI", "P/A Overlap", "G/L Distance",
                  "RP MirrorTree", "RP ContextTree", "Tree Distance",
                  "Gene Distance", "Moran's I", "Orientation MI",
                  "Gene Vector (10)", "Sequence Info (10)",
                  "Gene Vector (100)", "Sequence Info (100)",
                  "Gene Vector (500)", "Sequence Info (500)")
colnames(results) <- name_mapping
names(example_point) <- names(example_point2) <- name_mapping

## Reorder to the same ordering as in the paper
reorder_vec <- c(1,4,3,2,
                 6,5,7,
                 8,10,9,
                 12,14,16,
                 11,13,15)
results <- results[,reorder_vec]
example_point <- example_point[reorder_vec]
example_point2 <- example_point2[reorder_vec]

pdf(file.path(figdir, outfile), width=4.3*2,height=4.3*2.25, onefile=TRUE)
par(mar=c(6,3,2,1)+0.1, mgp=c(2,0.6,0), oma=c(2,0,0,0))
titles <- c("Time Complexity by Number of Taxa",
            "Time Complexity by Number of Pairs")
side_mgp <- c(2,0.6,0)
bottom_mgp <- c(2,0.85,0)
cols <- c('#45A649','#D81B60')
textcex <- 1
numcex <- 0.9
title_line <- -0.75
layout(matrix(1:4, nrow=2, byrow=TRUE))
offsets <- cumsum(rep(1.1,ncol(results)))
for(i in 1:2){
  barplot(results[i,], ylim=c(0,2.5), xlim=c(0,max(offsets)+1),
          space=0.1, width=1, yaxt='n', xaxt='n', xpd=TRUE, cex.axis=textcex,
          main='', ylab=ifelse(i==1, "Scaling Exponent", ''))
  title(main=titles[i], line=title_line)
  axis(side=1, at=offsets-0.5,tick=FALSE, labels=colnames(results),
       mgp=bottom_mgp,
       las=3, gap.axis=-1, line=-0.75, cex.axis=textcex, xpd=NA)
  axis(side=2, at=seq(0,2,0.5), mgp=side_mgp,
       las=1, gap.axis=-1, cex.axis=textcex, xpd=NA)
  for(j in seq_along(cols)){
    text(x=max(offsets)+0.6, y=j+0.005,
         labels=bquote("O(x"^.(j)*")"), cex=textcex, adj=c(0,0.5),
         col=cols[j], xpd=NA)
    lines(x=c(-0.5,max(offsets)+0.5), y=rep(j,2), col=cols[j], lwd=1.5, lty=3)
  }
  rect(xleft=offsets-1,xright=offsets,
       ytop=results[i,]+0.1,ybottom=results[i,]+0.01, col='white', border=NA)
  to_print <- results[i,]
  to_print[to_print >= 100] <- round(to_print)
  text(x=offsets-0.5,
       y=results[i,]+0.025,
       font=2,
       srt=90,
       adj=0,
       labels=sprintf("%.2f",to_print), cex=numcex)
}

examples <- rbind(example_point, example_point2)
for(i in 1:2){
  example_point <- examples[i,]
  inter_space <- 0.1
  bar_width <- 1
  offsets <- cumsum(rep(bar_width + inter_space,length(example_point)))
  barplot(example_point, ylim=c(0.025,2500), xlim=c(0,max(offsets)+1),
          space=inter_space, width=bar_width, xaxt='n', yaxt='n', xpd=TRUE, cex.axis=textcex,
          log='y', main='',
          ylab=ifelse(i==1, "Runtime (sec.)", ''))
  title(main=paste0('Runtime for ', NUM_TAXA, " Taxa and ",
                    ifelse(i==1, paste0(NUM_PAIRS, " Pairs"), "1 Pair")),
        line=title_line)
  axis(side=1, at=offsets-(bar_width/2),tick=FALSE, labels=colnames(results),
       mgp=bottom_mgp,
       las=3, gap.axis=-1, line=-0.75, cex.axis=textcex, xpd=NA)
  axis(side=2, at=c(0.05,0.1,1,10,100,500), labels=c("0.05","0.1","1","10","100","500"),
       mgp=side_mgp,
       las=1, gap.axis=-1, cex.axis=textcex, xpd=NA)
  time_breakpoints <- c(1,60)
  for(j in seq_along(time_breakpoints)){
    text(x=max(offsets)+1.5, y=time_breakpoints[j]*1.025,
         labels=ifelse(j==1, '1 sec.', "1 min."),
         col=ifelse(j==1, cols[1], cols[2]), cex=textcex, xpd=NA)
    lines(x=c(-0.5,max(offsets)+0.5), y=rep(time_breakpoints[j],2),
          col=cols[j], lwd=1.5, lty=3)
  }

  labs <- vapply(example_point, \(x){
    if(x >= 100){
      sprintf("%.0f", x)
    } else if(x >= 10){
      sprintf("%.1f",x)
    } else {
      sprintf("%.2f",x)
    }
  }, character(1L))
  yp <- log(example_point,base=10)
  yp <- yp + 0.05
  rect(xleft=offsets-1,xright=offsets,
       ytop=10**(yp+0.33),ybottom=10**yp, col='white', border=NA)
  text(x=offsets-0.5,
       y=10**(yp),
       font=2,
       srt=90,
       adj=0,
       labels=labs, cex=numcex)
}
dev.off(dev.list()['pdf'])
layout(1)
