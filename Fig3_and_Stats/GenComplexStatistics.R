load('ComplexPredsAllPairs.RData')

subpreds <- AllPairs[AllPairs$Category!=2,]
all_props <- subpreds$MinMembers
Pairings <- subpreds[,c("KO1", "KO2")]
subpreds$Jaccard <- subpreds$Jaccard*subpreds$PAPV
subpreds$MutualInformation <- subpreds$MutualInformation*subpreds$PAPV
subpreds$Hamming <- subpreds$Hamming*subpreds$PAPV
AllCats <- subpreds$Category
subpreds <- subpreds[,c(6,4,3,7,10,9,11,14,15,16,12,13)]
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

ComplexStatistics <- list()
source("PredictionCheck.R")
for(algo in colnames(subpreds)){
  print(algo)
  #nrow(subpreds[,algo])
  ComplexStatistics[[algo]] <- vcheckans(subpreds[,algo], key)
  ComplexStatistics[[algo]]$Color <- colors[algo]
  print(ComplexStatistics[[algo]]$AUROC)
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

# Logit
print("Logistic Regression")
predictvals <- numeric(nrow(testingResults))
for(i in 1:5){
  rowsin <- c(posv[possamp!=i], negv[negsamp!=i])
  rowtest <- c(posv[possamp==i], negv[negsamp==i])
  x <- glm(isTP ~ ., data=testingResults[rowsin,], family='binomial')
  pred <- predict(x, testingResults[rowtest,-anscol])
  predictvals[rowtest] <- pred
}
EnsembleComplexStatistics$Logit <- vcheckans(predictvals, key)
EnsembleComplexStatistics$Logit$Color <- 5L

# RandomForest
print("Random Forest")
library(randomForest)
predictvals <- numeric(nrow(testingResults))
testingResults$isTP <- as.factor(testingResults$isTP)
for(i in 1:5){
  rowsin <- c(posv[possamp!=i], negv[negsamp!=i])
  rowtest <- c(posv[possamp==i], negv[negsamp==i])
  x <- randomForest(isTP ~ ., data=testingResults[rowsin,])
  pred <- predict(x, testingResults[rowtest,-anscol], type='prob')[,'TRUE']
  predictvals[rowtest] <- pred
}
EnsembleComplexStatistics$RandomForest <- vcheckans(predictvals, key)
EnsembleComplexStatistics$RandomForest$Color <- 5L

# Neural Network
print("Neural Network")
library(neuralnet)
for(i in 1:5){
  rowsin <- c(posv[possamp!=i], negv[negsamp!=i])
  rowtest <- c(posv[possamp==i], negv[negsamp==i])
  x <- neuralnet(isTP ~ ., hidden=12, threshold=0.5, stepmax=1e6,
                 data=testingResults[rowsin,], linear.output=FALSE)
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
     file='ComplexStatistics.RData')
