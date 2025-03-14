plot_heatmaps <- FALSE

# precomputed statistics for each tree
load(file.path(datadir, "ModuleTreeStatistics.RData"))
load(file.path(datadir, "ModToKO.RData"))
load(file.path(datadir, 'ModulePredsAllPairs.RData'))

if(USE_KEGG){
  subpreds <- KEGGAllPairs
  outfile <- file.path(outdir, 'MulticlassModuleDataKEGG.RData')
} else {
  subpreds <- AllPairs
  outfile <- file.path(outdir, 'MulticlassModuleData.RData')
}
Pairings <- subpreds[,c("Mod1", "Mod2")]

subpreds$Category <- subpreds$NoKOOverlapCategory
subpreds <- subpreds[,c(10,9,8,6,4,3,5,12:14,11,15)]
subpreds[is.na(subpreds)] <- 0L
subpreds$Category <- AllPairs$NoKOOverlapCategory
Pairings <- Pairings[subpreds$Category>0,]
subpreds <- subpreds[subpreds$Category>0,]

ap1 <- tree_members[Pairings[,1]]
ap2 <- tree_members[Pairings[,2]]
all_props <- vapply(seq_along(ap1), \(x) min(ap1[x], ap2[x]), integer(1L))

n_fold <- 5L
testSets <- vector('list', n_fold)
groupings <- rep(0, nrow(subpreds))

cat("Building CV Folds:\n")
set.seed(724L)
# First assign each set to a group
allcats <- unique(subpreds$Category)
for(i in allcats){
  pos <- which(subpreds$Category == i)
  assignment <- rep(seq_len(n_fold), length.out=length(pos))
  assignment <- sample(assignment)
  groupings[pos] <- assignment
}

# Then create the train/test sets
for(i in seq_len(n_fold)){
  l <- list()
  pos <- which(groupings==i)
  l$test <- pos
  trainset_whole <- which(groupings!=i)
  trainsetsub <- subpreds[trainset_whole,]
  train_props <- all_props[trainset_whole]
  min_counts <- min(table(trainsetsub$Category))
  trainset <- rep(0, min_counts * length(allcats))

  # Finding the propensity of cat 1s
  allbreaks <- seq(0, 8600, by=100)
  num_bins <- length(allbreaks) - 1L
  count1s <- train_props[trainsetsub$Category==1L]
  tempbins <- .bincode(count1s, allbreaks)
  tab1 <- tabulate(tempbins, nbins=num_bins)

  for(j in seq_along(allcats)){
    subprops <- train_props[trainsetsub$Category==j]
    num_newbins <- .bincode(subprops, allbreaks)
    sampprobs <- (tab1 / tabulate(num_newbins, nbins=num_bins))[num_newbins]
    ss <- sample(which(trainsetsub$Category==j), min_counts, prob = sampprobs, replace=FALSE)
    trainset[((j-1)*min_counts+1):(j*min_counts)] <- trainset_whole[ss]
  }
  l$train <- trainset
  testSets[[i]] <- l
}
cat("Training Distribution:\n")
print(t(vapply(testSets, \(x) table(subpreds$Category[x$train]), allcats)))
cat("Testing Distribution:\n")
print(t(vapply(testSets, \(x) table(subpreds$Category[x$test]), allcats)))
cat("Totals:\n")
print(colSums(t(vapply(testSets, \(x) table(subpreds$Category[x$test]), allcats))))
print(table(subpreds$Category))

# Checking cdfs
cat("checking ecdfs of gene prevalence per fold:\n")
for(j in seq_along(testSets)){
  alltrain_ecdf <- testSets[[j]]$train
  alltrain_cat <- subpreds$Category[alltrain_ecdf]
  alltrain_minv <- vapply(alltrain_ecdf, \(x) min(tree_members[unlist(Pairings[x,])]), integer(1L))
  for(i in seq_along(allcats)){
    plot.ecdf(alltrain_minv[alltrain_cat==i], col=i, add=i!=1, pch='.', verticals=TRUE)
  }
}

cat("Testing on CV folds:\n")
set.seed(464L)
for(i in seq_along(testSets)){
  cat(i, '/', length(testSets), '\n')
  trainData <- subpreds[testSets[[i]]$train,]
  trainData$Category <- as.factor(trainData$Category)

  model <- randomForest(trainData[,1:12], trainData[,13], importance=TRUE, maxnodes = 100)
  #model <- randomForest(trainData[,14:18], trainData[,13], importance=TRUE)
  testData <- subpreds[testSets[[i]]$test,]
  actuals <- subpreds$Category[testSets[[i]]$test]

  preds <- predict(model, testData, type='prob')
  testSets[[i]]$model <- model
  testSets[[i]]$predictions <- preds
  testSets[[i]]$actuals <- actuals

  #print(model)
}

cat("Calculating feature importance:\n")
if(!USE_KEGG){
  fname <- file.path(figdir, "MainFigures", "3b_featimportance.pdf")
  pdf(file=fname, width=4, height=2.115, pointsize = 6)
  cvec <- c('#45A649','#D81B60', '#1E88E5','#FFC107')
  allconf <- vapply(testSets, \(x) importance(x$model, type=1L, scale=FALSE), numeric(nrow(testSets[[1]]$model$importance)))
  rownames(allconf) <- colnames(subpreds)[1:12]
  rownames(allconf) <- c("G/L Distance", "P/A Conservation",
                         "G/L MI", "P/A Jaccard",
                          "RP ContextTree", "RP MirrorTree", "Tree Distance",
                          "Gene Distance", "Moran's I", "Orientation MI",
                          "Gene Vector", "Sequence Info")
  allconf <- allconf[c(4,1,2,3,5:7,8,10,9,12,11),]
  layout(matrix(c(0,1), nrow=1), widths=c(0.15,1))
  w <- 0.4
  allconf <- allconf[rev(seq_len(nrow(allconf))),]
  barplot(rowMeans(allconf), beside=TRUE, names.arg=rownames(allconf), horiz=TRUE, las=2,
          cex.names=1, cex.axis=1, col=c(rev(rep(cvec[c(2,3,4,1)], times=c(4,3,3,2)))), space=w,
          xaxt='n', xlab='Mean Decrease in Accuracy', mgp=c(2,0.5,0), xlim=c(0,0.07))
  alabv <- seq(0,0.14,0.02)
  axis(side=1, at=alabv, labels=(paste0(alabv*100, '%')), mgp=c(3,0.5,0))
  yacc <- seq_len(nrow(allconf))*1 + seq(0, w*(nrow(allconf)-1), by=w) - 0.1
  for(i in seq_len(nrow(allconf))){
    #lines(y=rep(yacc[i], 2), x=c(min(allconf[i,]), max(allconf[i,])), lwd=2,xpd=NA)
    points(y=rep(yacc[i], ncol(allconf)), x=allconf[i,], pch=21, cex=0.75, xpd=NA)
  }
  dev.off()
  layout(matrix(1))
}

allpredictions <- matrix(0, nrow=nrow(subpreds), ncol=length(allcats))
for(i in seq_along(testSets)){
  allpredictions[testSets[[i]]$test,] <- testSets[[i]]$predictions
}

cat("Plotting confusion matrices by confidence:\n")
breaks <- c(0,0.5)
HeatmapsList <- vector('list', length(breaks))
if(!USE_KEGG){
  fnames <- c(file.path(figdir, "SupplFigures", "0ConfHeatmap.png"),
              file.path(figdir, "MainFigures", "3a_Heatmap.png"))
} else {
  fnames <- c(file.path(figdir, "SupplFigures", "0ConfKEGGHeatmap.png"),
              file.path(figdir, "SupplFigures", "50ConfKEGGHeatmap.png"))
}

for(j in seq_along(breaks)){
  cutoff <- breaks[j]
  print(cutoff)
  results_hm <- matrix(0, nrow=length(allcats), ncol=length(allcats))
  results_prob_hm <- matrix(0, nrow=length(allcats), ncol=length(allcats))
  key <- c("Adjacent", "Same Module", "Same Path", "Same GlobalPath", "None")
  colnames(results_hm) <- rownames(results_hm) <- key

  rh2 <- matrix(0, nrow=4,ncol=4)
  for(i in seq_len(nrow(allpredictions))){
    actual <- subpreds$Category[i]
    if(max(allpredictions[i,]) < cutoff) next
    choice <- which.max(allpredictions[i,])[1]
    results_hm[actual, choice] <- results_hm[actual,choice] + 1L
    results_prob_hm[actual,] <- results_prob_hm[actual,] + allpredictions[i,]
    rh2[min(actual,4), min(choice,4)] <- rh2[min(actual,4), min(choice,4)] + 1L
  }

  ## Plotting
  norm_results_hm <- results_hm / rowSums(results_hm)
  labs <- matrix(paste0(formatC(norm_results_hm*100, digits=1, format='f'), '%'),
                 nrow=nrow(norm_results_hm))
  htmcols <- c('#b5b5b5','#cbcbcb',
               hcl.colors(10, 'Terrain 2',rev=TRUE),
               '#02701b','#026318')
  htmcols <- c('#b5b5b5','#e0a683','#dba872','#cba764','#b7a555',
               '#97974e','#798642','#5e7b31','#3f6f22','#026318')
  output_plot <- pheatmap::pheatmap(#rh2 / rowSums(rh2),
    norm_results_hm,
    display_numbers = labs,
    number_color = 'black',
    cluster_cols=FALSE,
    cluster_rows=FALSE,
    labels_row = '',
    labels_col='',
    #main=paste0('Cutoff: ', cutoff, ' (n=',sum(results_hm), ')'),
    fontsize_number=26,
    number_format="%.3f",
    col=colorRampPalette(htmcols)(100),
    breaks=seq(0,1,0.01),
    drop_levels=FALSE,
    fontsize=20,
    legend=FALSE)
  ## EXTREMELY hacky workaround because the authors of `pheatmap`
  ## made some...interesting implementation decisions
  png(file=fnames[j], height=650, width=650)
  grid::grid.newpage()
  grid::grid.draw(output_plot$gtable)
  dev.off()
  HeatmapsList[[j]] <- results_hm
  print(rowSums(results_hm))
}
names(HeatmapsList) <- paste0(breaks*100, '%')

# This is the example for Fig. 5
Fig5Row <- which(Pairings$Mod1=='M133_1_1' &
                        Pairings$Mod2=='M339_4_3')
cat("Fig 5 prediction probabilities:\n")
print(allpredictions[Fig5Row,])

save(testSets, subpreds, Pairings,
     allconf, allpredictions, HeatmapsList, Fig5Row,
     file=outfile)
if(USE_KEGG){
  KEGGAllPairs$EnsemblePredictedClass <- vapply(seq_len(nrow(allpredictions)), \(x) which.max(allpredictions[x,]), integer(1L))
  KEGGAllPairs$EnsemblePredictedProb <- vapply(seq_len(nrow(allpredictions)), \(x) max(allpredictions[x,]), numeric(1L))
} else {
  AllPairs$EnsemblePredictedClass <- vapply(seq_len(nrow(allpredictions)), \(x) which.max(allpredictions[x,]), integer(1L))
  AllPairs$EnsemblePredictedProb <- vapply(seq_len(nrow(allpredictions)), \(x) max(allpredictions[x,]), numeric(1L))
}
cat("Saving updated AllPairs object...\n")
save(AllPairs, KEGGAllPairs, file=file.path(datadir, 'ModulePredsAllPairs.RData'))

if(!USE_KEGG){
  cat("Building network for Fig 3c:\n")
  ClusteringData <- AllPairs[!is.na(AllPairs$EnsemblePredictedClass),]
  nams <- unique(c(ClusteringData$Mod1, ClusteringData$Mod2))
  mprob <- m <- matrix(0, nrow=length(nams), ncol=length(nams))
  colnames(m) <- rownames(m) <- nams
  colnames(mprob) <- rownames(mprob) <- nams
  cutoff <- 0

  for(i in seq_len(nrow(ClusteringData))){
    m1 <- ClusteringData$Mod1[i]
    m2 <- ClusteringData$Mod2[i]
    r <- ClusteringData$EnsemblePredictedClass[i]
    p <- ClusteringData$EnsemblePredictedProb[i]
    if(p > 0.5)
      m[m1,m2] <- m[m2,m1] <- r
    mprob[m1,m2] <- mprob[m2,m1] <- p
  }
  newrn <- vapply(rownames(m), \(x){
    paste(ModMap[[x]], collapse=';')
  }, character(1L))
  rownames(m) <- colnames(m) <- newrn
  rownames(mprob) <- colnames(mprob) <- newrn

  testMat <- m
  for(i in seq_len(nrow(testMat))){
    r <- rep(0, ncol(testMat))
    r2 <- mprob[i,]
    r2[m[i,] != 1] <- 0
    r2[which(colnames(testMat) == rownames(testMat)[i])] <- 0
    r2p1 <- r2p2 <- NULL
    if(!all(r2 == 0)){
      # get top two direct connections
      r2p1 <- which.max(r2)
      r2[r2p1] <- 0
      r2p2 <- which.max(r2)
      r[c(r2p1,r2p2)] <- 1
    }
    # add back in "same module"s
    r[m[i,] == 2] <- 2
    testMat[i,] <- r
  }
  em <- testMat
  em[em==2] <- 0.5


  gmexp <- graph_from_adjacency_matrix(em, mode='max', weighted=FALSE)
  set.seed(125L)
  grps <- cluster_louvain(gmexp)
  finalgrps <- communities(grps)
  finalgrps <- lapply(finalgrps, sort)
  finalgrps <- lapply(finalgrps, unique)
  finalgrps <- finalgrps[order(lengths(finalgrps), decreasing=TRUE)]
  finalgrps <- finalgrps[lengths(finalgrps) > 1]
  Fig3Ex <- finalgrps[[sample(length(finalgrps), 1)]]
  Fig3Mat <- testMat[Fig3Ex, Fig3Ex]
  posOnes <- which(Fig3Mat==1, arr.ind=TRUE)
  rownames(posOnes) <- NULL
  posOnes <- cbind(rownames(Fig3Mat)[posOnes[,1]], colnames(Fig3Mat)[posOnes[,2]])
  for(i in seq_len(nrow(posOnes)))
    posOnes[i,] <- sort(posOnes[i,])
  posOnes <- unique(posOnes)
  posOnes[order(posOnes[,1], posOnes[,2]),]
  sort(unique(c(posOnes)))
  conversion <- c(
    "pigA", "pigJ", "pigH", "pigM", "pigF", #both modules 2-6, 1 is a complex
    "pigD", "pigE", "pigB", #M00837_2_1-3,
    "redJ", "redL", "redK"
    )
  names(conversion) <- sort(unique(c(posOnes)))
  posOnes[,1] <- conversion[posOnes[,1]]
  posOnes[,2] <- conversion[posOnes[,2]]
  gfig3ex <- graph_from_edgelist(posOnes, directed = FALSE)
  colmapping <- list(c("pigA", "pigJ", "pigH", "pigM", "pigF"),
                     c("pigD", "pigE", "pigB"),
                     c('redJ',"redL", "redK"))
  V(gfig3ex)$color <- vapply(names(V(gfig3ex)), \(x) {which(vapply(colmapping, \(y) x%in%y, logical(1L)))}, integer(1L))
  # plot them in a circle
  nr <- nrow(posOnes)
  vertex_pos <- cbind(rep(1:4, each=3), rep(1:4, 3))
  (plot(gfig3ex, vertex.label.cex=0.5, vertex.cex=4, edge.width=2))

  fullConnections <- m[Fig3Ex, Fig3Ex]
  rownames(fullConnections) <- colnames(fullConnections) <- conversion[rownames(fullConnections)]

  cat("Building list of misclasses:\n")
  to_write <- cbind(AllPairs, allpredictions)
  colnames(to_write)[(-4:0) + ncol(to_write)] <- paste0("EnsembleProb", as.character(1:5))
  to_write <- to_write[order(allpredictions[,1],
                             allpredictions[,2],
                             allpredictions[,3],
                             allpredictions[,4],
                             allpredictions[,5], decreasing = TRUE),]
  to_write <- to_write[to_write$EnsemblePredictedClass != to_write$NoKOOverlapCategory,]
  to_write$KO1Orig <- vapply(to_write$Mod1Orig, \(x) paste(unique(unlist(backmapping[x])), collapse=','), character(1L))
  to_write$KO2Orig <- vapply(to_write$Mod2Orig, \(x) paste(unique(unlist(backmapping[x])), collapse=','), character(1L))
  to_write$Mod1Orig <- vapply(to_write$Mod1Orig, paste, character(1L),collapse=',')
  to_write$Mod2Orig <- vapply(to_write$Mod2Orig, paste, character(1L),collapse=',')
  to_write <- to_write[,c(1:2,19:20,21,18,23,29:30,3:6,8:15,24:28)]
  write.csv(to_write, row.names = FALSE,
            file = file.path(outdir, "MisclassesByEnsembleProbability.csv"))
}
