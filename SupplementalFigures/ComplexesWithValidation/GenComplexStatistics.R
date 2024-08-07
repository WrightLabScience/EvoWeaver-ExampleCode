load(file.path(datadir, 'ComplexPredsAllPairs.RData'))

if(USE_KEGG){
  pastev <- 'KEGG'
  subpreds <- KEGGAllPairs[KEGGAllPairs$Category!=2,]
} else {
  pastev <- ''
  subpreds <- AllPairs[AllPairs$Category!=2,]
}
all_props <- subpreds$MinMembers
Pairings <- subpreds[,c("KO1", "KO2")]

AllCats <- subpreds$Category
subpreds <- subpreds[,3:14]
colors <- rep(1:4, times=c(4,3,3,2))
names(colors) <- colnames(subpreds)

set.seed(6698L)
posfalse <- which(AllCats == 3)
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

LRModels <- RFModels <- NNModels <- vector('list', 5L)

cat("evaluating component algorithms...\n")
ComplexStatistics <- list()
for(algo in colnames(subpreds)){
  ComplexStatistics[[algo]] <- vcheckans(subpreds[,algo], key)
  ComplexStatistics[[algo]]$Color <- colors[algo]
}

EnsembleComplexStatistics <- list()

testingResults <- subpreds
testingResults$isTP <- key
testingResults[is.na(testingResults)] <- 0
posv <- which(key)
negv <- which(!key)
possamp <- sample(1:5, length(posv), replace=TRUE)
negsamp <- sample(1:5, length(negv), replace=TRUE)
anscol <- which(colnames(testingResults) == 'isTP')

cat("training ensemble models...\n")
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
EnsembleComplexStatistics$Logit <- vcheckans(predictvals, key)
EnsembleComplexStatistics$Logit$Color <- 5L

# RandomForest
print("Random Forest")
predictvals <- numeric(nrow(testingResults))
testingResults$isTP <- as.factor(testingResults$isTP)
for(i in 1:5){
  rowsin <- c(posv[possamp!=i], negv[negsamp!=i])
  rowtest <- c(posv[possamp==i], negv[negsamp==i])
  x <- randomForest(isTP ~ ., data=testingResults[rowsin,], maxnodes=25)
  RFModels[[i]] <- x
  pred <- predict(x, testingResults[rowtest,-anscol], type='prob')[,'TRUE']
  predictvals[rowtest] <- pred
}
EnsembleComplexStatistics$RandomForest <- vcheckans(predictvals, key)
EnsembleComplexStatistics$RandomForest$Color <- 5L

# Neural Network
print("Neural Network")
for(i in 1:5){
  rowsin <- c(posv[possamp!=i], negv[negsamp!=i])
  rowtest <- c(posv[possamp==i], negv[negsamp==i])
  x <- neuralnet(isTP ~ ., hidden=12, threshold=0.5, stepmax=1e6,
                 data=testingResults[rowsin,], linear.output=FALSE)
  NNModels[[i]] <- x
  pred <- predict(x, testingResults[rowtest, -anscol])[,1]
  predictvals[rowtest] <- pred
}
EnsembleComplexStatistics$NeuralNetwork <- vcheckans(1-predictvals, key)
EnsembleComplexStatistics$NeuralNetwork$Color <- 5L
print(vapply(EnsembleComplexStatistics, \(x) x$AUROC, numeric(1L)))
print(vapply(ComplexStatistics, \(x) x$AUROC, numeric(1L)))

RawScores <- subpreds
RawScores$isTP <- key

# This file is used in Fig2_FigS1
save(ComplexStatistics, EnsembleComplexStatistics, RawScores,
     file=file.path(outdir, paste0('ComplexStatistics', pastev, '.RData')))

save(LRModels, RFModels, NNModels,
     file=file.path(outdir, paste0("ComplexEnsembleModels", pastev, ".RData")))

if(!USE_KEGG){
  cat("evaluating extended ensemble models...\n")
  ## Extra testing for other Ensemble methods
  n <- names(EnsembleComplexStatistics)
  n[n=="RandomForest"] <- "RandomForest25"
  n[n=="NeuralNetwork"] <- "NeuralNetwork1"
  names(EnsembleComplexStatistics) <- n
  set.seed(955L)
  cat("random forest nodes: ")
  for(maxnode in c(0, 2, 5, 10, 100)){
    cat(maxnode, ' ')
    predictvals <- numeric(nrow(testingResults))
    for(i in 1:5){
      rowsin <- c(posv[possamp!=i], negv[negsamp!=i])
      rowtest <- c(posv[possamp==i], negv[negsamp==i])
      if(maxnode==0)
        x <- randomForest(isTP ~ ., data=testingResults[rowsin,])
      else
        x <- randomForest(isTP ~ ., data=testingResults[rowsin,], maxnodes=maxnode)
      #  RFModels[[i]] <- x
      pred <- predict(x, testingResults[rowtest,-anscol], type='prob')[,'TRUE']
      predictvals[rowtest] <- pred
    }
    EnsembleComplexStatistics[[paste0("RandomForest", ifelse(maxnode==0, "Inf", maxnode))]] <- vcheckans(predictvals, key)
  }
  cat('\n')
  print(vapply(EnsembleComplexStatistics, \(x) x$AUROC, numeric(1L)))

  cat("svm kernels: ")
  for(kernel in c("radial", "linear", "sigmoid", "polynomial")){
    cat(kernel, ' ')
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
    EnsembleComplexStatistics[[paste0("SVM", kernel)]] <- vcheckans(predictvals, key)
  }
  cat('\n')

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
    EnsembleComplexStatistics[[paste0("NeuralNetwork", nlayer)]] <- vcheckans(1-predictvals, key)
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
    EnsembleComplexStatistics[[paste0("NeuralNetworkBig", nlayer)]] <- vcheckans(1-predictvals, key)
  }
  cat('\n')
  print(vapply(EnsembleComplexStatistics, \(x) x$AUROC, numeric(1L)))

  save(EnsembleComplexStatistics,
       file=file.path(outdir, 'ExtendedComplexEnsembleMethods.RData'))
}
