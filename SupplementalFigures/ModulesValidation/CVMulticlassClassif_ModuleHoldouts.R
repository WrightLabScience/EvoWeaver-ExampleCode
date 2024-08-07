outfile <- 'MulticlassModuleData_ModuleHoldout.RData'
featimportance_filename <- 'featimportance_moduleholdouts.pdf'

#plot_heatmaps <- TRUE

# precomputed statistics for each tree
load(file.path(datadir, "Modules", "ModuleTreeStatistics.RData"))
load(file.path(datadir, "Modules", 'ModulePredsAllPairs.RData'))
subpreds <- AllPairs
Pairings <- subpreds[,c("Mod1", "Mod2")]
subpreds <- subpreds[,c(6,8:10, 3:5, 12:14, 11,15)]
subpreds[is.na(subpreds)] <- 0L
subpreds$Category <- AllPairs$NoKOOverlapCategory
Pairings <- Pairings[subpreds$Category>0,]
subpreds <- subpreds[subpreds$Category>0,]

ap1 <- tree_members[Pairings[,1]]
ap2 <- tree_members[Pairings[,2]]
all_props <- vapply(seq_along(ap1), \(x) min(ap1[x], ap2[x]), integer(1L))
allcats <- unique(subpreds$Category)

## Get a list of all modules -- we'll sample from these for our test group
all_modules <- unique(unlist(Pairings))
all_main_modules <- gsub("(M[0-9]{3})_.*", '\\1', all_modules)
all_main_modules <- match(all_main_modules, unique(all_main_modules))
names(all_main_modules) <- all_modules
num_main_modules <- length(unique(all_main_modules))

n_fold <- 10L
testSets <- vector('list', n_fold)

set.seed(123L)
## assign each module to a group
module_partitioning <- sample(n_fold, num_main_modules, replace=TRUE)
module_partitioning <- module_partitioning[all_main_modules]
names(module_partitioning) <- names(all_main_modules)
groupings <- cbind(module_partitioning[Pairings[[1]]], module_partitioning[Pairings[[2]]])
rownames(groupings) <- NULL

# Then create the train/test sets
for(i in seq_len(n_fold)){
  l <- list()
  mods <- names(module_partitioning[module_partitioning==i])
  ## find the test set -- all pairs involving the module
  pos <- which(groupings[,1]==i | groupings[,2] == i)
  l$test <- pos
  trainset_whole <- which(groupings[,1]!=i & groupings[,2] != i)
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

  model <- randomForest(trainData[,1:12], trainData[,13], importance=TRUE, maxnodes=100L)
  testData <- subpreds[testSets[[i]]$test,]
  actuals <- subpreds$Category[testSets[[i]]$test]

  # random forest
  preds <- predict(model, testData, type='prob')

  testSets[[i]]$model <- model
  testSets[[i]]$predictions <- preds
  testSets[[i]]$actuals <- actuals
}

allconf <- vapply(testSets, \(x) importance(x$model, type=1L, scale=FALSE), numeric(nrow(testSets[[1]]$model$importance)))
rownames(allconf) <- colnames(subpreds)[1:12]
rownames(allconf) <- c("P/A Jaccard", "G/L MI",
                       "P/A Conservation", "G/L Distance",
                        "RP ContextTree", "RP MirrorTree", "Tree Distance",
                        "Gene Distance", "Moran's I", "Transcription MI",
                        "Gene Vector", "Sequence Info")
allconf <- allconf[rev(seq_len(nrow(allconf))),]
allconf <- allconf[c(2,1,4,3,5,6,7,8,9,11,10,12),]
if(PLOT_FEAT_IMPORTANCE_MULTICLASS){
  fname <- file.path(figdir, featimportance_filename)
  pdf(file=fname, width=3.5, height=2, pointsize = 6)
  cvec <- c('#45A649','#D81B60', '#1E88E5','#FFC107')
  layout(matrix(c(0,1), nrow=1), widths=c(0.15,1))
  w <- 0.4
  barplot(rowMeans(allconf), beside=TRUE, names.arg=rownames(allconf), horiz=TRUE, las=2,
          cex.names=1, cex.axis=1, col=c(rev(rep(cvec[c(2,3,4,1)], times=c(4,3,3,2)))), space=w,
          xaxt='n', xlab='Mean Decrease in Accuracy', mgp=c(2,0.5,0), xlim=c(0,0.07))
  alabv <- seq(0,0.14,0.02)
  axis(side=1, at=alabv, labels=(paste0(alabv*100, '%')), mgp=c(3,0.5,0))
  yacc <- seq_len(nrow(allconf))*1 + seq(0, w*(nrow(allconf)-1), by=w) - 0.1
  for(i in seq_len(nrow(allconf))){
    lines(y=rep(yacc[i], 2), x=c(min(allconf[i,]), max(allconf[i,])), lwd=2,xpd=NA)
  }

  dev.off()
  layout(matrix(1))
}

allpredictions <- matrix(0, nrow=nrow(subpreds), ncol=length(allcats))
counts_per_row <- integer(nrow(allpredictions))
for(i in seq_along(testSets)){
  allpredictions[testSets[[i]]$test,] <- allpredictions[testSets[[i]]$test,] + testSets[[i]]$predictions
  counts_per_row[testSets[[i]]$test] <- counts_per_row[testSets[[i]]$test] + 1L
}
allpredictions <- allpredictions / counts_per_row

fnames <- c(file.path(figdir, "0ConfModuleHoldoutHeatmap.png"),
            file.path(figdir, "50ConfModuleHoldoutHeatmap.png"))
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

  norm_results_hm <- results_hm / rowSums(results_hm)
  labs <- matrix(paste0(formatC(norm_results_hm*100, digits=1, format='f'), '%'),
                 nrow=nrow(norm_results_hm))
  htmcols <- c('#b5b5b5','#cbcbcb',
               hcl.colors(10, 'Terrain 2',rev=TRUE),
               '#02701b','#026318')
  htmcols <- c('#b5b5b5','#e0a683','#dba872','#cba764','#b7a555',
               '#97974e','#798642','#5e7b31','#3f6f22','#026318')
  output_plot <- pheatmap::pheatmap(
    norm_results_hm,
    display_numbers = labs,
    number_color = 'black',
    cluster_cols=FALSE,
    cluster_rows=FALSE,
    labels_row = '',
    labels_col='',
    legend=FALSE,
    fontsize_number=26,
    number_format="%.03f",
    col=colorRampPalette(htmcols)(100),
    breaks=seq(0,1,0.01),
    drop_levels=FALSE,
    fontsize=20)
  ## EXTREMELY hacky workaround because the authors of `pheatmap`
  ## made some...interesting implementation decisions
  png(file=fnames[j], height=650, width=650)
  grid::grid.newpage()
  grid::grid.draw(output_plot$gtable)
  dev.off()
  HeatmapsList[[j]] <- results_hm
}
names(HeatmapsList) <- paste0(breaks*100, '%')

# This is the example for Fig. 4
Fig5Row <- which(Pairings$Mod1=='M133_1_1' &
                        Pairings$Mod2=='M339_4_3')[1L]

save(testSets, subpreds, Pairings,
     allconf, allpredictions, HeatmapsList, Fig5Row,
     file=file.path(outdir, outfile))
