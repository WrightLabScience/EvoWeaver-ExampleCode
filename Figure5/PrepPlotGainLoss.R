PrepPlotGainLoss <- function(dend, attnum=1L,
                         pch, ptcex, lty=c(1,1,1), lwd=c(1,1),
                         cols=c('#D81B60', '#1CFE63', 'black')){
  require(SynExtend)
  dendrapply(dend, \(x){
    ai <- attnum
    s <- attr(x, 'FitchState')
    l <- list()
    if(is.leaf(x)){
      l$col <- ifelse(s[ai]==0L, cols[1], cols[2])
      l$pch <- pch
      l$cex <- ptcex
      if(all(s == 1L))
        l$cex <- 0
      attr(x, 'nodePar') <- l
    } else {
      for(i in seq_along(x)){
        sc <- attr(x[[i]], 'FitchState')
        l <- list()
        if(s[ai] != sc[ai]){
          l$col <- ifelse(s[ai] == 1L, cols[1], cols[2])
          l$lwd <- lwd[1]
          if(!all(s != sc))
            l$lty <- lty[3]
          else
            l$lty <- lty[1]
        } else {
          l$col <- cols[3]
          l$lwd <- lwd[2]
          l$lty <- lty[2]
        }
        attr(x[[i]], 'edgePar') <- l
      }
    }
    x
  }, how='post')
}
