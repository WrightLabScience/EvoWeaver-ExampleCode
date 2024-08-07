outfile <- 'SXX_ExampleTreesPP.pdf'
TEXT_CEX <- 0.75
CEX_MULT <- 1.33
pdf(file=file.path(figdir, outfile), onefile=TRUE, width=4.3*2,height=4.3*2)
par(mai=c(0,0,0,0)+0.2, mgp=c(0,0,0), oma=c(2,3,2,0))
layout(matrix(1:30, nrow=6), widths=c(3, rep(1, 4L)))

set.seed(777L)
n_taxa <- 32
d <- dist(seq_len(n_taxa))
d[] <- runif(length(d))
tree <- TreeLine(myDistMatrix=d, method='NJ', verbose=FALSE)
ctr <- 1L
tree <- dendrapply(tree, \(x){
  if(is.leaf(x)){
    attr(x,'label') <- as.character(ctr)
    attr(x, 'height') <- max(0, attr(x,'height')-0.01)
    ctr <<- ctr + 1
  }
  x
}, how='post')

labs <- labels(tree)
ll <- length(labs)
# 5 scenarios we want to compare:
#   1. high conservation
#   2. Low Conservation
#   3. mutually exclusive
#   4. 50/50 conservation, balanced across tree
#   5. 50/50 conservation alternating
#   6. random pattern

vecs <- list(
  # high conservation
  a1=seq_len(ll),
  a2=seq_len(ll),

  # mutually exclusive
  b1=seq(1,ll,by=2),
  b2=seq(2,ll,by=2),

  # balanced 50/50 by clade
  c1=seq_len(attr(tree[[1]], 'members'))[-8],
  c2=seq_len(attr(tree[[1]], 'members'))[-1],

  # balanced 50/50 alternating (pair w b1)
  d1=seq(1,ll,by=2),

  # random
  e1=sample(ll, 11),
  e2=sample(ll, 11),

  f1=-3:0 + ll,
  f2=-1:0 + ll
)

vecs <- lapply(vecs, \(x) labs[x])

p1 <- c(1,10,3,5,7,8)
p2 <- c(2,11,4,6,3,9)

descrip <- c("High Conservation", "Low Conservation", "Mutual Exclusion",
             "Clade Balance", "Alternating Balance", "Random")
y <- .Call("initCDend", tree, PACKAGE='SynExtend')
for(i in seq_along(p1)){
  v1 <- vecs[[p1[i]]]
  v2 <- vecs[[p2[i]]]
  pa1 <- .Call("calcGainLoss", y, v1, FALSE, PACKAGE="SynExtend")
  pa2 <- .Call("calcGainLoss", y, v2, FALSE, PACKAGE="SynExtend")
  v3 <- 2L*pa1 + pa2 + 1L
  v4 <- .Call("cladeCollapsePA", y, v3)
  ctr <- 1L
  cols <- c("grey",'#2B6DA8', '#D81B60',"#824484")
  pdend <- dendrapply(tree, \(x){
    attr(x, "nodePar") <- list(col=cols[v3[ctr]], pch=20)
    if(!is.leaf(x)){
      curcol <- attr(x,"nodePar")$col
      rcol <- attr(x[[1]],"nodePar")$col
      lcol <- attr(x[[2]],"nodePar")$col

      if(rcol != curcol){
        newcol <- cols[bitwXor(match(rcol, cols)-1L, match(curcol,cols)-1L) + 1L]
        attr(x[[1]],"edgePar")$col <- ifelse(newcol=='grey', 'black', newcol)
      }
      if(lcol != curcol){
        newcol <- cols[bitwXor(match(lcol, cols)-1L, match(curcol,cols)-1L) + 1L]
        attr(x[[2]],"edgePar")$col <- ifelse(newcol=='grey', 'black', newcol)
      }
    }
    ctr <<- ctr + 1L
    x
  }, how='post.order')
  #plot(pdend, main=paste0("Scenario ", i, ": ", descrip[i]), yaxt='n', leaflab='none')
  plot(pdend, main='', cex.main=0.75, yaxt='n', leaflab='none',
       horiz=FALSE)
  pdend <- dendrapply(pdend, \(x){attr(x, 'edgePar') <- list(lty=0); x})
  par(new=TRUE)
  plot(pdend, main='', cex.main=0.75, yaxt='n', leaflab='none',
       horiz=FALSE)
  mtext(descrip[i], side=2, font=2, cex=TEXT_CEX)

}
rm(y)

legendcols <- c("grey",'#2B6DA8', '#D81B60',"#824484","black")
legend('bottomright', pch=c(rep(20,4), rep(NA, 4)),
        lty=c(rep(0,4), rep(1,4)), pt.cex=1.5, lwd=1.5,
        col=legendcols[c(2:3,1,4,2:3,5,4)], cex=TEXT_CEX*CEX_MULT, ncol=4,
        legend=c("Gene 1 Present", "Gene 2 Present",
                 "Both Absent", "Both Present",
                 "Change in Gene 1", "Change in Gene 2",
                 "Change in Neither Gene", "Change in Both Genes"),
        xpd=NA, inset=c(-1.25,-0.525), bty='n')
cols_legend <- c('#D81B60', '#E0A608', '#2B6DA8', '#45A649', 'grey40', 'black')
score_colors <- cols_legend[c(4,5,1)]
score_choices <- c("Positive", "Weak", "Negative")
signif_colors <- cols_legend[c(3,2)]
signif_choices <- c("High", "Low")
alg_choices <- c("P/A Jaccard", "P/A Overlap", "G/L MI", "G/L Distance")

alg_scores <- matrix(c(1,1,3,1,1,2,  # jaccard
                       1,1,3,1,1,2,  # overlap
                       2,1,3,1,2,2,  # MI
                       2,1,3,1,1,2), # gainloss
                     ncol=4)
alg_signif <- matrix(c(2,2,1,2,1,2,  # jaccard
                       2,1,1,1,1,2,  # overlap
                       2,2,1,2,1,2,  # MI
                       2,1,1,1,1,2), # gainloss
                     ncol=4)

pw <- EvoWeaver(vecs, MySpeciesTree=tree, NoWarn = TRUE)
res <- predict(pw, Method=c("PAJaccard", "PAOverlap", "GLMI", "GLDistance"),
               Subset=cbind(p1,p2), Verbose=FALSE)[c(1,6,2,4,3,5),-c(1:2)]

SIGNIF_POS <- c(0.5,0.5)
SCORE_POS <- c(0.5,0.7)
RESULT_POS <- c(0.5,0.3)
for(i in seq_len(4L)){
  for(j in seq_len(6L)){
    plot(NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE, frame.plot=FALSE,
         xlab='', ylab='')
    if(j == 1){
      title(main=alg_choices[i], cex.main=1.25)
    }
    # score text
    score_c <- alg_scores[j,i]
    text(x=SCORE_POS[1], y=SCORE_POS[2],
         labels=paste0(score_choices[score_c], " Score"),
         col=score_colors[score_c], xpd=NA, font=2, cex=TEXT_CEX*CEX_MULT)

    signif_c <- alg_signif[j,i]
    # signif text
    text(x=SIGNIF_POS[1], y=SIGNIF_POS[2],
         labels=paste0(signif_choices[signif_c], " Significance"),
         col=signif_colors[signif_c], xpd=NA, font=2, cex=TEXT_CEX*CEX_MULT)

    actual_score <- res[j,i]
    text(x=RESULT_POS[1], y=RESULT_POS[2],
         labels=paste0("Result: ", ifelse(actual_score<0, ' ',''), sprintf("%.02f", actual_score)),
         col='black', xpd=NA, font=2, cex=TEXT_CEX*CEX_MULT)
  }
}

dev.off(dev.list()['pdf'])
