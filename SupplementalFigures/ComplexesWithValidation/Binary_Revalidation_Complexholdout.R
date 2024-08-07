## This script will redo binary complex classification using complex holdouts
load(file.path(datadir, 'ComplexPredsAllPairs.RData'))
load(file.path(datadir, "KEGGModuleComplexes.RData"))

## I'm going to hold out based on modules, which should have the same effect as holding out complexes
## each complex should always appear in the same module, so holding out groups that appear in
## the same module will be the same as holding out entire complexes (if not better)

NUM_CV_FOLD <- 10L
subpreds <- AllPairs[AllPairs$Category!=2,]
all_props <- subpreds$MinMembers
Pairings <- subpreds[,c("KO1", "KO2")]
AllCats <- subpreds$Category
subpreds <- subpreds[,3:14]
colors <- rep(1:4, times=c(4,3,3,2))
names(colors) <- colnames(subpreds)

set.seed(891L)
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
subPairings <- Pairings[c(postrue,posfalse),]
key <- rep(FALSE, nrow(subpreds))
key[seq_along(postrue)] <- TRUE

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

all_kos <- unique(unlist(subPairings))
subcompMapping <- ComplexMapping[ComplexMapping[,1] %in% all_kos,1:2]
all_modules <- unique(unlist(subcompMapping[,2]))
modulesamp <- sample(seq_len(NUM_CV_FOLD), length(all_modules), replace=TRUE)
names(modulesamp) <- all_modules
komapping <- modulesamp[subcompMapping[,2]]
names(komapping) <- subcompMapping[,1]

anscol <- which(colnames(testingResults) == 'isTP')

groupings <- cbind(komapping[subPairings[[1]]],
                    komapping[subPairings[[2]]])

# Names of all ensemble functions to test
function_names <- c("Logit",
                    "RandomForest2", "RandomForest5", "RandomForest10",
                    "RandomForest25", "RandomForest100", "RandomForestInf",
                    "NeuralNetwork1", "NeuralNetwork2",
                    "NeuralNetwork3", "NeuralNetwork4",
                    "NeuralNetworkBig1", "NeuralNetworkBig2",
                    "NeuralNetworkBig3", "NeuralNetworkBig4",
                    "SVMlinear", "SVMradial", "SVMpolynomial", "SVMsigmoid")

# train functions for each predictor
train_functions <- c(\(r) suppressWarnings(glm(isTP ~ ., data=testingResults[r,], family='binomial')),
                     \(r){
                       tmpdata <- testingResults
                       tmpdata$isTP <- as.factor(tmpdata$isTP)
                       randomForest(isTP ~ ., data=tmpdata[r,], maxnodes=2)
                     },
                     \(r){
                       tmpdata <- testingResults
                       tmpdata$isTP <- as.factor(tmpdata$isTP)
                       randomForest(isTP ~ ., data=tmpdata[r,], maxnodes=5)
                     },
                     \(r){
                       tmpdata <- testingResults
                       tmpdata$isTP <- as.factor(tmpdata$isTP)
                       randomForest(isTP ~ ., data=tmpdata[r,], maxnodes=10)
                     },
                     \(r){
                       tmpdata <- testingResults
                       tmpdata$isTP <- as.factor(tmpdata$isTP)
                       randomForest(isTP ~ ., data=tmpdata[r,], maxnodes=25)
                     },
                     \(r){
                       tmpdata <- testingResults
                       tmpdata$isTP <- as.factor(tmpdata$isTP)
                       randomForest(isTP ~ ., data=tmpdata[r,], maxnodes=100)
                     },
                     \(r){
                       tmpdata <- testingResults
                       tmpdata$isTP <- as.factor(tmpdata$isTP)
                       randomForest(isTP ~ ., data=tmpdata[r,])
                     },
                     \(r) neuralnet(isTP ~ ., hidden=rep(12,1), threshold=0.5, stepmax=1e6,
                                    data=testingResults[r,], linear.output=FALSE),
                     \(r) neuralnet(isTP ~ ., hidden=rep(12,2), threshold=0.5, stepmax=1e6,
                                    data=testingResults[r,], linear.output=FALSE),
                     \(r) neuralnet(isTP ~ ., hidden=rep(12,3), threshold=0.5, stepmax=1e6,
                                    data=testingResults[r,], linear.output=FALSE),
                     \(r) neuralnet(isTP ~ ., hidden=rep(12,4), threshold=0.5, stepmax=1e6,
                                    data=testingResults[r,], linear.output=FALSE),
                     \(r) neuralnet(isTP ~ ., hidden=rep(24,1), threshold=0.5, stepmax=1e6,
                                    data=testingResults[r,], linear.output=FALSE),
                     \(r) neuralnet(isTP ~ ., hidden=rep(24,2), threshold=0.5, stepmax=1e6,
                                    data=testingResults[r,], linear.output=FALSE),
                     \(r) neuralnet(isTP ~ ., hidden=rep(24,3), threshold=0.5, stepmax=1e6,
                                    data=testingResults[r,], linear.output=FALSE),
                     \(r) neuralnet(isTP ~ ., hidden=rep(24,4), threshold=0.5, stepmax=1e6,
                                    data=testingResults[r,], linear.output=FALSE),
                     \(r){
                       tmp <- testingResults
                       tmp$isTP <- as.factor(tmp$isTP)
                       svm(isTP ~ ., data=tmp[r,], probability=TRUE, kernel='linear')
                     },
                     \(r){
                       tmp <- testingResults
                       tmp$isTP <- as.factor(tmp$isTP)
                       svm(isTP ~ ., data=tmp[r,], probability=TRUE, kernel='radial')
                     },
                     \(r){
                       tmp <- testingResults
                       tmp$isTP <- as.factor(tmp$isTP)
                       svm(isTP ~ ., data=tmp[r,], probability=TRUE, kernel='polynomial')
                     },
                     \(r){
                       tmp <- testingResults
                       tmp$isTP <- as.factor(tmp$isTP)
                       svm(isTP ~ ., data=tmp[r,], probability=TRUE, kernel='sigmoid')
                     })

# predict() methods for each predictor
predict_functions <- c(\(r) predict(x, testingResults[r,-anscol]),
                       \(r) predict(x, testingResults[r,-anscol], type='prob')[,'TRUE'],
                       \(r) predict(x, testingResults[r,-anscol], type='prob')[,'TRUE'],
                       \(r) predict(x, testingResults[r,-anscol], type='prob')[,'TRUE'],
                       \(r) predict(x, testingResults[r,-anscol], type='prob')[,'TRUE'],
                       \(r) predict(x, testingResults[r,-anscol], type='prob')[,'TRUE'],
                       \(r) predict(x, testingResults[r,-anscol], type='prob')[,'TRUE'],
                       \(r) predict(x, testingResults[r, -anscol])[,1],
                       \(r) predict(x, testingResults[r, -anscol])[,1],
                       \(r) predict(x, testingResults[r, -anscol])[,1],
                       \(r) predict(x, testingResults[r, -anscol])[,1],
                       \(r) predict(x, testingResults[r, -anscol])[,1],
                       \(r) predict(x, testingResults[r, -anscol])[,1],
                       \(r) predict(x, testingResults[r, -anscol])[,1],
                       \(r) predict(x, testingResults[r, -anscol])[,1],
                       \(r) attr(predict(x, testingResults[r, -anscol], probability=TRUE), "probabilities")[,"TRUE"],
                       \(r) attr(predict(x, testingResults[r, -anscol], probability=TRUE), "probabilities")[,"TRUE"],
                       \(r) attr(predict(x, testingResults[r, -anscol], probability=TRUE), "probabilities")[,"TRUE"],
                       \(r) attr(predict(x, testingResults[r, -anscol], probability=TRUE), "probabilities")[,"TRUE"])

## test functions on each CV fold
predictvals <- lapply(seq_along(train_functions), \(x) numeric(nrow(testingResults)))
num_seen <- integer(nrow(testingResults))
for(i in seq_len(NUM_CV_FOLD)){
  cat("Fold", i, "/", NUM_CV_FOLD, '\n')
  posvin <- !vapply(posv, \(j) i %in% groupings[j,], logical(1L))
  negvin <- !vapply(negv, \(j) i %in% groupings[j,], logical(1L))
  rowsin <- c(posv[which(posvin)], negv[which(negvin)])
  rowtest <- c(posv[which(!posvin)], negv[which(!negvin)])
  for(j in seq_along(train_functions)){
    #cat('\t-', function_names[j], '\n')
    x <- (train_functions[[j]])(rowsin)
    pred <- (predict_functions[[j]])(rowtest)
    predictvals[[j]][rowtest] <- predictvals[[j]][rowtest] + pred
  }
  num_seen[rowtest] <- num_seen[rowtest] + 1L
}

# store results in statistics object
nsz <- num_seen > 0
for(i in seq_along(function_names)){
  p <- predictvals[[i]]
  p <- p[nsz] / num_seen[nsz]
  EnsembleComplexStatistics[[function_names[i]]] <- vcheckans(p, key[nsz])
  EnsembleComplexStatistics[[function_names[i]]]$Color <- 5L
}

print(vapply(ComplexStatistics, \(x) x$AUROC, numeric(1L)))
print(vapply(EnsembleComplexStatistics, \(x) x$AUROC, numeric(1L)))

RawScores <- subpreds
RawScores$isTP <- key

# This is used in Fig2_FigS1
save(ComplexStatistics, EnsembleComplexStatistics, RawScores,
     file=file.path(outdir, 'ComplexStatistics_fullcomplexholdouts.RData'))

