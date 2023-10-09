load('ModulePredsAllPairs.RData')
subpreds <- AllPairs[!AllPairs$HasComplex,]
Pairings <- subpreds[,c("Mod1", "Mod2")]
subpreds$Jaccard <- subpreds$Jaccard*subpreds$PAPV
subpreds$MutualInformation <- subpreds$MutualInformation*subpreds$PAPV
subpreds$Hamming <- subpreds$Hamming*subpreds$PAPV
AllCats <- subpreds$NoKOOverlapCategory
subpreds <- subpreds[,c(10,6,11,8,4,3,5,13:15,12,16)]
colors <- rep(1:4, times=c(4,3,3,2))
names(colors) <- colnames(subpreds)
set.seed(6698L)
#posfalse <- sample(which(AllCats<=4), sum(AllCats==1))
all_props <- vapply(seq_len(nrow(subpreds)), \(x) min(tree_members[Pairings$Mod1[x]], tree_members[Pairings$Mod2[x]]), integer(1L))
posfalse <- which(AllCats <= 4)
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
source("PredictionCheck.R")
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
  x <- glm(isTP ~ ., data=testingResults[rowsin,], family='binomial')
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
  x <- randomForest(isTP ~ ., data=testingResults[rowsin,])
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
  pred <- predict(x, testingResults[rowtest, -anscol])[,1]
  predictvals[rowtest] <- pred
}
EnsembleBlockStatistics$NeuralNetwork <- vcheckans(1-predictvals, key)
EnsembleBlockStatistics$NeuralNetwork$Color <- 5L
print(vapply(EnsembleBlockStatistics, \(x) x$AUROC, numeric(1L)))

RawScores <- subpreds
RawScores$isTP <- key

# This is used in Fig2_FigS1
save(BlockStatistics, EnsembleBlockStatistics, RawScores,
     file='ModuleStatistics.RData')

