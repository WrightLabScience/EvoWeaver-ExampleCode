## This script should not be run directly!
## this is called by GenerateDataAndFig_Fig2.R

load(file.path(datadir, "Modules", 'ModulePredsAllPairs.RData'))
load(file.path(datadir, "Modules", 'ModuleTreeStatistics.RData'))

if(IS_KEGG){
  subpreds <- KEGGAllPairs
} else {
  subpreds <- AllPairs
}
Pairings <- subpreds[,c("Mod1", "Mod2")]

AllCats <- subpreds$NoKOOverlapCategory
subpreds <- subpreds[,c(10,9,8,6,4,3,5,12:14,11,15)]
colors <- rep(1:4, times=c(4,3,3,2))
names(colors) <- colnames(subpreds)

RFModels <- LRModels <- NNModels <- vector('list', 5L)
set.seed(6698L)
#posfalse <- sample(which(AllCats<=4), sum(AllCats==1))
all_props <- vapply(seq_len(nrow(subpreds)), \(x) min(tree_members[Pairings$Mod1[x]], tree_members[Pairings$Mod2[x]]), integer(1L))
posfalse <- which(AllCats == 4)
postrue <- which(AllCats==1)

allbreaks <- seq(0, 8600, by=100)
num_bins <- length(allbreaks) - 1L
count1s <- all_props[postrue]
tempbins <- .bincode(count1s, allbreaks)
tab1 <- tabulate(tempbins, nbins=num_bins)

subprops <- all_props[posfalse]
num_newbins <- .bincode(subprops, allbreaks)
sampprobs <- (tab1 / tabulate(num_newbins, nbins=num_bins))[num_newbins]
ss <- sample(posfalse, length(postrue), prob = sampprobs, replace=FALSE)
posfalse <- ss

subpreds <- subpreds[c(postrue,posfalse),]
key <- rep(FALSE, nrow(subpreds))
key[seq_along(postrue)] <- TRUE

BlockStatistics <- list()
for(algo in colnames(subpreds)){
  print(algo)
  #nrow(subpreds[,algo])
  BlockStatistics[[algo]] <- vcheckans(subpreds[,algo], key)
  BlockStatistics[[algo]]$Color <- colors[algo]
  print(BlockStatistics[[algo]]$AUROC)
}

EnsembleBlockStatistics <- list()

testingResults <- subpreds
testingResults$isTP <- key
testingResults[is.na(testingResults)] <- 0
posv <- which(key)
negv <- which(!key)
possamp <- sample(1:5, length(posv), replace=TRUE)
negsamp <- sample(1:5, length(negv), replace=TRUE)
anscol <- which(colnames(testingResults) == 'isTP')

# Logit
print("Logistic Regression")
predictvals <- numeric(nrow(testingResults))
for(i in 1:5){
  rowsin <- c(posv[possamp!=i], negv[negsamp!=i])
  rowtest <- c(posv[possamp==i], negv[negsamp==i])
  x <- suppressWarnings(glm(isTP ~ ., data=testingResults[rowsin,], family='binomial'))
  LRModels[[i]] <- x
  pred <- predict(x, testingResults[rowtest,-anscol])
  predictvals[rowtest] <- pred
}
EnsembleBlockStatistics$Logit <- vcheckans(predictvals, key)
EnsembleBlockStatistics$Logit$Color <- 5L

# RandomForest
print("Random Forest")
library(randomForest)
predictvals <- numeric(nrow(testingResults))
testingResults$isTP <- as.factor(testingResults$isTP)
for(i in 1:5){
  rowsin <- c(posv[possamp!=i], negv[negsamp!=i])
  rowtest <- c(posv[possamp==i], negv[negsamp==i])
  x <- randomForest(isTP ~ ., data=testingResults[rowsin,], maxnodes=25) # consider maxnodes=25 to reduce overfitting
  RFModels[[i]] <- x
  pred <- predict(x, testingResults[rowtest,-anscol], type='prob')[,'TRUE']
  predictvals[rowtest] <- pred
}
EnsembleBlockStatistics$RandomForest <- vcheckans(predictvals, key)
EnsembleBlockStatistics$RandomForest$Color <- 5L

# Neural Network
print("Neural Network")
library(neuralnet)
for(i in 1:5){
  rowsin <- c(posv[possamp!=i], negv[negsamp!=i])
  rowtest <- c(posv[possamp==i], negv[negsamp==i])
  x <- neuralnet(isTP ~ ., hidden=12, threshold=0.5, stepmax=1e6,
                 data=testingResults[rowsin,], linear.output=FALSE)
  NNModels[[i]] <- x
  pred <- predict(x, testingResults[rowtest, -anscol])[,1]
  predictvals[rowtest] <- pred
}
EnsembleBlockStatistics$NeuralNetwork <- vcheckans(1-predictvals, key)
EnsembleBlockStatistics$NeuralNetwork$Color <- 5L
print(vapply(BlockStatistics, \(x) x$AUROC, numeric(1L)))
print(vapply(EnsembleBlockStatistics, \(x) x$AUROC, numeric(1L)))

RawScores <- subpreds
RawScores$isTP <- key

# This is used in Fig2_FigS1
save(BlockStatistics, EnsembleBlockStatistics, RawScores,
     file=file.path(outdir, ifelse(IS_KEGG, "ModuleStatisticsKEGG.RData", 'ModuleStatistics.RData')))
save(LRModels, NNModels, RFModels,
     file=file.path(outdir, ifelse(IS_KEGG, "EnsembleModelsKEGG.RData", 'EnsembleModels.RData')))

## Extra testing for other Ensemble methods
n <- names(EnsembleBlockStatistics)
n[n=="RandomForest"] <- "RandomForest25"
n[n=="NeuralNetwork"] <- "NeuralNetwork1"
names(EnsembleBlockStatistics) <- n
set.seed(955L)
for(maxnode in c(0, 2, 5, 10, 100)){
  cat(maxnode, 'nodes...\n')
  predictvals <- numeric(nrow(testingResults))
  for(i in 1:5){
    rowsin <- c(posv[possamp!=i], negv[negsamp!=i])
    rowtest <- c(posv[possamp==i], negv[negsamp==i])
    if(maxnode == 0)
      x <- randomForest(isTP ~ ., data=testingResults[rowsin,])
    else
      x <- randomForest(isTP ~ ., data=testingResults[rowsin,], maxnodes=maxnode)
    #  RFModels[[i]] <- x
    pred <- predict(x, testingResults[rowtest,-anscol], type='prob')[,'TRUE']
    predictvals[rowtest] <- pred
  }
  EnsembleBlockStatistics[[paste0("RandomForest", ifelse(maxnode==0, "Inf", maxnode))]] <- vcheckans(predictvals, key)
}
print(vapply(EnsembleBlockStatistics, \(x) x$AUROC, numeric(1L)))

library(e1071)
for(kernel in c("radial", "linear", "sigmoid", "polynomial")){
  predictvals <- numeric(nrow(testingResults))
  testingResults$isTP <- as.factor(testingResults$isTP)
  for(i in 1:5){
    rowsin <- c(posv[possamp!=i], negv[negsamp!=i])
    rowtest <- c(posv[possamp==i], negv[negsamp==i])
    x <- svm(isTP ~ ., data=testingResults[rowsin,], kernel=kernel, probability=TRUE)
    #  RFModels[[i]] <- x
    pred <- predict(x, testingResults[rowtest,-anscol], probability=TRUE)
    predictvals[rowtest] <- attr(pred, 'probabilities')[,'TRUE']
  }
  EnsembleBlockStatistics[[paste0("SVM", kernel)]] <- vcheckans(predictvals, key)
}

cat("neural network hidden layers:")
for(nlayer in seq(2,4)){
  cat(" ", nlayer)
  for(i in 1:5){
    rowsin <- c(posv[possamp!=i], negv[negsamp!=i])
    rowtest <- c(posv[possamp==i], negv[negsamp==i])
    x <- neuralnet(isTP ~ ., hidden=rep(12,nlayer), threshold=0.5, stepmax=1e6,
                   data=testingResults[rowsin,], linear.output=FALSE)
    pred <- predict(x, testingResults[rowtest, -anscol])[,1]
    predictvals[rowtest] <- pred
  }
  EnsembleBlockStatistics[[paste0("NeuralNetwork", nlayer)]] <- vcheckans(1-predictvals, key)
}
cat('\n')
cat("neural network, larger hidden layers:")
for(nlayer in seq(1,4)){
  cat(" ", nlayer)
  for(i in 1:5){
    rowsin <- c(posv[possamp!=i], negv[negsamp!=i])
    rowtest <- c(posv[possamp==i], negv[negsamp==i])
    x <- neuralnet(isTP ~ ., hidden=rep(24,nlayer), threshold=0.5, stepmax=1e6,
                   data=testingResults[rowsin,], linear.output=FALSE)
    pred <- predict(x, testingResults[rowtest, -anscol])[,1]
    predictvals[rowtest] <- pred
  }
  EnsembleBlockStatistics[[paste0("NeuralNetworkBig", nlayer)]] <- vcheckans(1-predictvals, key)
}
cat('\n')

print(vapply(EnsembleBlockStatistics, \(x) x$AUROC, numeric(1L)))

save(EnsembleBlockStatistics,
     file=file.path(outdir,
                    ifelse(IS_KEGG, "ExtendedModuleEnsembleMethodsKEGG.RData",
                           'ExtendedModuleEnsembleMethods.RData')))
