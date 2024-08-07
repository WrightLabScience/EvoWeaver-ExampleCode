BuildCaseStudyForInputNetworkNoWeight <- function(AllPairs, subpreds, allpredictions,
                                          Modules, Legend,
                                          ActualGraph, marginlabel, plotHeader,
                                          ## these extra args are used in supplemental figs
                                          VERT_AXIS=0.8, HORIZ_AXIS=1.25,
                                          V_SIZE=13, E_WIDTH=1.5, VLC=0.9,
                                          highlight_labels=NULL, wrong_lty=2,
                                          LDIST_MULT=2.45,
                                          TITLE_OFFSET=1.0, MARGIN_OFFSET=0,
                                          GRAYLINE_VALS=c(-2.25,20.5), VERTNAMESCENTER=FALSE){
  ## This function assumes that a dev is already open with the correct layout()
  ## it will neither open nor close a dev, we're just writing plots to the screen
  ## Inputs:
  ##  - AllPairs, subpreds, allpredictions: data from module classification
  ##  - Modules: Output of find_all_blocks, with KEGGToGene filled in
  ##  - Namemap: (named character) conversion from module codes to gene names
  ##  - Legend: (named character) colors to use for each module point
  ##  - ActualGraph: (int matrix) connections in KEGG
  USE_DIRECTED <- FALSE
  MAX_DEGREE <- 1L
  require(igraph)
  encoded_genes <- Modules$Encoded
  comm <- Modules$Blocks
  mapping <- Modules$EncodeToKEGG


  p <- (AllPairs$Mod1 %in% encoded_genes) & (AllPairs$Mod2 %in% encoded_genes)
  commScores <- subpreds[p,]
  ensemble_scores <- rowSums(allpredictions[p,1,drop=F])
  #ensemble_scores <- colSums(t(allpredictions[p,1:2,drop=F]) * c(1,0.5))
  #ensemble_scores <- colSums(t(allpredictions[p,]) * c(0.814, 0.137, 0.015, 0.022, 0.011))
  OrigModules <- AllPairs[p, c(1,2,19,20)]
  m <- matrix(NA_real_, nrow=length(comm), ncol=length(comm))
  mpos <- matrix(NA_integer_, nrow=length(comm), ncol=length(comm))
  for(i in seq_len(nrow(commScores))){
    p1 <- mapping[OrigModules$Mod1[i]]
    p2 <- mapping[OrigModules$Mod2[i]]
    p1pos <- which(p1 == comm)
    p2pos <- which(p2 == comm)
    mpos[p1pos,p2pos] <- mpos[p2pos,p1pos] <- i
  }

  rownames(mpos) <- colnames(mpos) <- Modules$KEGGToGene[comm]


  m <- matrix(0, nrow=nrow(mpos), ncol=ncol(mpos))
  rownames(m) <- colnames(m) <- rownames(mpos)

  PercentileRank <- vector('list', ncol(subpreds)-1)
  names(PercentileRank) <- colnames(subpreds)[seq_along(PercentileRank)]
  for(i in seq_len(ncol(subpreds)-1)){
    m[] <- 0
    diag(m) <- 1
    q <- commScores[,i]
    rng <- unlist(subpreds[,i])
    rng <- rng[!is.na(rng)]
    q <- vapply(q, \(x) sum(x >= rng[]) / length(rng), numeric(1L))
    m[] <- q[mpos]
    m[is.na(m)] <- 0
    diag(m) <- 1
    PercentileRank[[i]] <- m
    #heatmap(m, Rowv=NA, Colv=NA, main=colnames(subpreds)[i], scale="none",
    #        symm=TRUE)
  }
  m[] <- 0
  m[] <- ensemble_scores[mpos]
  m[is.na(m)] <- 0
  diag(m) <- 1
  PercentileRank[[ncol(subpreds)]] <- m

  rownames(mpos)
  nr <- nrow(mpos)
  layout_mat <- t(vapply(seq_len(nr), \(x) c(sin((x-1)/(nr)*2*pi), cos((x-1)/(nr)*2*pi)),
                         numeric(2L)))
  layout_mat[,1] <- layout_mat[,1] * HORIZ_AXIS
  layout_mat[,2] <- layout_mat[,2] * VERT_AXIS
  rownames(layout_mat) <- rownames(mpos)

  ActualGraph[] <- rownames(mpos)[c(ActualGraph)]
  actual_g <- graph_from_edgelist(ActualGraph, TRUE)

  all_deg <- igraph::degree(actual_g, mode='out')
  #all_deg[] <- round(mean(igraph::degree(actual_g, mode='out')))
  if(!is.null(MAX_DEGREE))
    all_deg[] <- MAX_DEGREE
  colnames(subpreds)
  pred_categories <- list(1:4, 5:7, 8:10, 11:12, 1:12, 13)
  pred_categories <- list(1:4, 5:7, 8:10, 11:12, 13)
  lpc <- length(pred_categories)+1L
  catnames <- c("Component Predictors", "",
                "", "", #"Consensus",
                "Ensemble Predictions", "KEGG Connections")
  cols <- c('#D81B60', '#E0A608', '#2B6DA8', '#45A649', "#4c4c4c", '#000000')
  total_edgelist <- data.frame(V1=NULL, V2=NULL, color=NULL, id=NULL)
  for(j in seq_along(pred_categories)){
    TestMat <- Reduce(`+`, PercentileRank[pred_categories[[j]]]) / length(pred_categories[[j]])
    diag(TestMat) <- 0
    for(i in seq_len(nrow(TestMat))){
      rnam <- rownames(TestMat)[i]
      cutoff_ind <- max(all_deg[rnam], 1)
      #    cutoff_ind <- 2
      v <- TestMat[i,]
      tmp_val <- ifelse(cutoff_ind==0, Inf, sort(v,TRUE)[cutoff_ind])
      v[v<tmp_val] <- 0
      v[v>=tmp_val] <- 1
      TestMat[i,] <- v
    }

    adjPos <- which(TestMat==TRUE, arr.ind=TRUE)
    adjPos[] <- rownames(TestMat)[c(adjPos)]
    if(!USE_DIRECTED){
      adjPos <- t(apply(adjPos, 1, sort))
      adjPos <- unique(adjPos)
    }
    adjPos <- as.data.frame(adjPos)
    colnames(adjPos) <- c("from", "to")
    rownames(adjPos) <- NULL
    adjPos$color <- cols[j]
    adjPos$id <- j
    total_edgelist <- rbind(total_edgelist, adjPos)
  }
  tmp <- as.data.frame(ActualGraph)
  tmp <- cbind(tmp, 'black', length(pred_categories)+1L)
  #total_edgelist$width <- 0.5
  colnames(tmp) <- colnames(total_edgelist)
  total_edgelist <- rbind(total_edgelist, tmp)
  V(actual_g)
  v_degs <- ((0:nrow(mpos[-1,]) / nrow(mpos)) * 2*pi) - pi
  v_degs <- v_degs + pi/2
  v_degs[1] <- -0.9*pi
  v_fills <- Legend

  names(v_degs) <- names(v_fills) <- rownames(mpos)
  for(i in seq_len(lpc)){
    tmpdf <- total_edgelist[total_edgelist$id==i,]
    if(i != lpc){
      cmplines <- ActualGraph
      if(!USE_DIRECTED){
        tmpdf[,1:2] <- t(apply(tmpdf[,1:2], 1, sort))
        tmpdf <- tmpdf[!duplicated(tmpdf[,1:2]),]
        cmplines <- t(apply(ActualGraph, 1, sort))
        cmplines <- unique(cmplines)
      }
      tmpdf$lty <- wrong_lty
      for(j in seq_len(nrow(tmpdf))){
        if(any((tmpdf[j,1] == cmplines[,1] & tmpdf[j,2] == cmplines[,2]) |
          (tmpdf[j,2] == cmplines[,1] & tmpdf[j,1] == cmplines[,2]))){
          tmpdf$lty[j] <- 1
        }
      }
    }
    g <- graph_from_data_frame(tmpdf,
                               directed=(i==lpc) || USE_DIRECTED,
                               vertices=names(V(actual_g)))
    par(xpd=NA)
    if(i < lpc-1){
      plot(g, vertex.label=NA,
           vertex.color=v_fills[names(V(actual_g))],
           vertex.size=V_SIZE*1.25, layout=t(t(layout_mat[names(V(actual_g)),]) * c(1.3,1.1)),
           edge.arrow.size=0.1, edge.width=E_WIDTH,
           xpd=NA, rescale=FALSE)
    } else {
      vvv <- V(actual_g)
      vl_col <- rep("grey30", length(vvv))
      names(vl_col) <- names(vvv)
      if(i!=lpc && (!VERTNAMESCENTER)){
        vl_col[] <- "#00000000"
      }
      if(!is.null(highlight_labels)){
        vl_col[highlight_labels] <- cols[1]
      }
      plot(g, vertex.label.cex=VLC,
           vertex.label.color=vl_col,
           vertex.color=v_fills[names(V(actual_g))],
           vertex.label.dist=2.8+LDIST_MULT*abs(layout_mat[,1]),
           vertex.label.degree=c(-pi/2,v_degs[names(V(actual_g))][-1])-(pi/16),
           vertex.size=V_SIZE, layout=t(t(layout_mat[names(V(actual_g)),]) * c(1,0.9)),
           edge.arrow.size=0.35, vertex.label.font=2, edge.width=E_WIDTH,
           xpd=NA, rescale=FALSE)
    }
    if(plotHeader && i %in% c(1, lpc-1,lpc)){
      mtext(catnames[i], cex=0.75, font=2, xpd=NA, adj=ifelse(i==1,-2+TITLE_OFFSET,0.5), padj=-1)
    }
    if(i==1){
      text(x=-2.75+MARGIN_OFFSET, y=-1, marginlabel,cex=1, xpd=NA, font=2, srt=90)
      if(!plotHeader)
        lines(x=GRAYLINE_VALS, y=rep(1.0625,2), xpd=NA, col='grey80', lty=1)
    }
  }
}
