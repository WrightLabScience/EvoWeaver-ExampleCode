library(randomForest)
outfile <- 'MulticlassModuleData.RData'

plot_heatmaps <- TRUE
outdir <- '/Users/aidan/Nextcloud/RStudioSync/KEGGStuff/ModuleAllPairsJob'

# precomputed statistics for each tree
load("ModuleTreeStatistics.RData")
load('ModulePredsAllPairs.RData', v=T)
subpreds <- AllPairs[!AllPairs$HasComplex,]
Pairings <- subpreds[,c("Mod1", "Mod2")]
subpreds$Jaccard <- subpreds$Jaccard*subpreds$PAPV
subpreds$MutualInformation <- subpreds$MutualInformation*subpreds$PAPV
subpreds$Hamming <- subpreds$Hamming*subpreds$PAPV
subpreds <- subpreds[,c(10,6,11,8,4,3,5,13:15,12,16)]
subpreds[is.na(subpreds)] <- 0L
subpreds$Category <- AllPairs$NoKOOverlapCategory[!AllPairs$HasComplex]
#subpreds$Category <- AllPairs$Category
Pairings <- Pairings[subpreds$Category>0,]
subpreds <- subpreds[subpreds$Category>0,]

ap1 <- tree_members[Pairings[,1]]
ap2 <- tree_members[Pairings[,2]]
all_props <- vapply(seq_along(ap1), \(x) min(ap1[x], ap2[x]), integer(1L))

n_fold <- 5L
testSets <- vector('list', n_fold)
groupings <- rep(0, nrow(subpreds))

set.seed(123L)
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
t(vapply(testSets, \(x) table(subpreds$Category[x$train]), allcats))
cat("Testing Distribution:\n")
t(vapply(testSets, \(x) table(subpreds$Category[x$test]), allcats))
cat("Totals:\n")
colSums(t(vapply(testSets, \(x) table(subpreds$Category[x$test]), allcats)))
table(subpreds$Category)

# Checking cdfs
#alltrain_ecdf <- unlist(lapply(testSets, \(x) x$train))
for(j in seq_along(testSets)){
  alltrain_ecdf <- testSets[[j]]$train
  alltrain_cat <- subpreds$Category[alltrain_ecdf]
  alltrain_minv <- vapply(alltrain_ecdf, \(x) min(tree_members[unlist(Pairings[x,])]), integer(1L))
  for(i in seq_len(5)){
    plot.ecdf(alltrain_minv[alltrain_cat==i], col=i, add=i!=1, pch='.', verticals=TRUE)
  }
}

for(i in seq_along(testSets)){
  cat(i, '/', length(testSets), '\n')
  trainData <- subpreds[testSets[[i]]$train,]
  trainData$Category <- as.factor(trainData$Category)

  #model <- randomForest(trainData[,1:12], trainData[,13], importance=TRUE)
  model <- randomForest(trainData[,14:18], trainData[,13], importance=TRUE)
  testData <- subpreds[testSets[[i]]$test,]
  actuals <- subpreds$Category[testSets[[i]]$test]

  preds <- predict(model, testData, type='prob')
  testSets[[i]]$model <- model
  testSets[[i]]$predictions <- preds
  testSets[[i]]$actuals <- actuals

  print(model)
}

fname <- 'FeatImportance.pdf'
pdf(file=fname, width=3.5, height=2, pointsize = 6)
cvec <- c('#45A649','#D81B60', '#1E88E5','#FFC107')
allconf <- vapply(testSets, \(x) importance(x$model, type=1L, scale=FALSE), numeric(nrow(testSets[[1]]$model$importance)))
rownames(allconf) <- colnames(subpreds)[1:12]
rownames(allconf) <- c("G/L Correlation", "P/A Jaccard",
                       "G/L Distance", "P/A MI",
                        "RP ContextTree", "RP MirrorTree", "Tree Distance",
                        "Gene Distance", "Moran's I", "Transcription MI",
                        "Gene Vector", "Sequence Info")
layout(matrix(c(0,1), nrow=1), widths=c(0.15,1))
w <- 0.4
allconf <- allconf[rev(seq_len(nrow(allconf))),]
allconf <- allconf[c(2,1,4,3,5,6,7,8,9,11,10,12),]
barplot(rowMeans(allconf), beside=TRUE, names.arg=rownames(allconf), horiz=TRUE, las=2,
        cex.names=1, cex.axis=1, col=c(rev(rep(cvec[c(2,3,4,1)], times=c(4,3,3,2)))), space=w,
        xaxt='n', xlab='Mean Decrease in Accuracy', mgp=c(2,0.5,0), xlim=c(0,0.15))
alabv <- seq(0,0.14,0.02)
axis(side=1, at=alabv, labels=(paste0(alabv*100, '%')), mgp=c(3,0.5,0))
yacc <- seq_len(nrow(allconf))*1 + seq(0, w*(nrow(allconf)-1), by=w) - 0.1
for(i in seq_len(nrow(allconf))){
  lines(y=rep(yacc[i], 2), x=c(min(allconf[i,]), max(allconf[i,])), lwd=2,xpd=NA)
}
#lines(x=rep(1/nrow(allconf),2), y=c(0,yacc[nrow(allconf)]+1+1.5*w), lty=4, xpd=NA)
dev.off()
layout(matrix(1))

allpredictions <- matrix(0, nrow=nrow(subpreds), ncol=length(allcats))
for(i in seq_along(testSets)){
  allpredictions[testSets[[i]]$test,] <- testSets[[i]]$predictions
}

breaks <- c(0,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,0.95)
breaks <- c(0,0.5)
HeatmapsList <- vector('list', length(breaks))
for(j in seq_along(breaks)){
  cutoff <- breaks[j]
  print(cutoff)
  results_hm <- matrix(0, nrow=length(allcats), ncol=length(allcats))
  results_prob_hm <- matrix(0, nrow=length(allcats), ncol=length(allcats))
  key <- c("Adjacent", "Same Module", "Same Path", "Same GlobalPath", "None")
  colnames(results_hm) <- rownames(results_hm) <- key

  rh2 <- matrix(0, nrow=4,ncol=4)
  for(i in seq_len(nrow(allpredictions))){
    #actual <- actualcats[i]
    actual <- subpreds$Category[i]
    if(max(allpredictions[i,]) < cutoff) next
    choice <- which.max(allpredictions[i,])[1]
    results_hm[actual, choice] <- results_hm[actual,choice] + 1L
    results_prob_hm[actual,] <- results_prob_hm[actual,] + allpredictions[i,]
    rh2[min(actual,4), min(choice,4)] <- rh2[min(actual,4), min(choice,4)] + 1L
  }

  if(plot_heatmaps){
    norm_results_hm <- results_hm / rowSums(results_hm)
    labs <- matrix(paste0(as.character(round(norm_results_hm*100,1)), '%'),
                   nrow=nrow(results_hm))
    htmcols <- c('#b5b5b5','#cbcbcb',
                 hcl.colors(10, 'Terrain 2',rev=TRUE),
                 '#02701b','#026318')
    htmcols <- c('#b5b5b5','#e0a683','#dba872','#cba764','#b7a555',
                 '#97974e','#798642','#5e7b31','#3f6f22','#026318')
    pheatmap::pheatmap(#rh2 / rowSums(rh2),
      norm_results_hm,
      #display_numbers = labs,
      number_color = 'black',
      cluster_cols=FALSE,
      cluster_rows=FALSE,
      labels_row = '',
      labels_col='',
      #main=paste0('Cutoff: ', cutoff, ' (n=',sum(results_hm), ')'),
      fontsize_number=16,
      number_format="%.3f",
      col=colorRampPalette(htmcols)(100),
      breaks=seq(0,1,0.01),
      drop_levels=FALSE,
      fontsize=20)
  } else {
    print(round(results_hm / rowSums(results_hm), 2))
  }
  HeatmapsList[[j]] <- results_hm
}
names(HeatmapsList) <- paste0(breaks*100, '%')

# This is the example for Fig. 4
Fig4Row <- which(Pairings$Mod1%in%c('M133_1_1', "M134_1_1") &
                        Pairings$Mod2==c('M339_4_3'))[1L]
allpredictions[Fig4Row,]

save(testSets, subpreds, Pairings,
     allconf, allpredictions, HeatmapsList, Fig4Row,
     file=outfile)

v <- rep(NA, nrow(AllPairs))
v[!AllPairs$HasComplex & AllPairs$NoKOOverlapCategory > 0] <-
  vapply(seq_len(nrow(allpredictions)), \(x) which.max(allpredictions[x,]), integer(1L))
AllPairs$EnsemblePredictedClass <- v

v <- rep(NA, nrow(AllPairs))
v[!AllPairs$HasComplex & AllPairs$NoKOOverlapCategory > 0] <-
  vapply(seq_len(nrow(allpredictions)), \(x) max(allpredictions[x,]), numeric(1L))
AllPairs$EnsemblePredictedProb <- v

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
  if(p > cutoff)
    m[m1,m2] <- m[m2,m1] <- r
}
newrn <- vapply(rownames(m), \(x){
  paste(ModMap[[x]], collapse=';')
}, character(1L))
rownames(m) <- colnames(m) <- newrn

em <- m
em[em > 2] <- 0L
em[em>0] <- 1L
library(igraph)
set.seed(106L)
gmexp <- graph_from_adjacency_matrix(em, mode='undirected', weighted=TRUE)
grps <- cluster_label_prop(gmexp)
finalgrps <- communities(grps)
finalgrps <- lapply(finalgrps, sort)
finalgrps <- lapply(finalgrps, unique)

# random group for Fig 3c
Fig3Ex <- finalgrps[[sample(seq_along(finalgrps), 1L)]]

Fig3Mat <- m[Fig3Ex, Fig3Ex]
posOnes <- which(Fig3Mat==1, arr.ind=TRUE)
for(i in seq_len(nrow(posOnes)))
  posOnes[i,] <- sort(posOnes[i,])
posOnes <- unique(posOnes)
posOnes[] <- rownames(Fig3Mat)[posOnes[]]
rownames(posOnes) <- NULL
posOnes

# Let's get a list of all the misclasses
tmp <- vapply(seq_len(nrow(allpredictions)), \(x) which.max(allpredictions[x,]), integer(1L))
locmispreds <- which(subpreds$Category > 3 & tmp%in%c(1L,2L) & allpredictions[,1]>=0.5)
mispreds <- cbind(subpreds[locmispreds,], Pairings[locmispreds,], rowSums(allpredictions[locmispreds,1:2]), allpredictions[locmispreds,1], allpredictions[locmispreds,2])
colnames(mispreds) <- c(colnames(subpreds), "Mod1", "Mod2", "CombinedProb", 'ProbDirect', 'ProbSameModule')
mispreds <- mispreds[order(mispreds$CombinedProb, decreasing=TRUE),]
head(mispreds, n=15)
mispreds <- cbind(mispreds, AllPairs[rownames(mispreds),c("Mod1Orig", "Mod2Orig")])
mispreds$Mod1Orig <- vapply(mispreds$Mod1Orig, paste, character(1L),collapse=',')
mispreds$Mod2Orig <- vapply(mispreds$Mod2Orig, paste, character(1L),collapse=',')
mcdirname <- 'CV_misclass.csv'
write.csv(mispreds[seq_len(15),c(18,19,14,15,13,16,17)], file=file.path(mcdirname, 'Misclass_prop6.csv'))
